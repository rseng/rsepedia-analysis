[provide general introduction to the issue and why it is relevant to this repository]

## Context of the issue

[provide more detailed introduction to the issue itself and why it is relevant]

[the remaining entries are only necessary if you are reporting a bug]

## Process to reproduce the issue

[ordered list the process to finding and recreating the issue, example below]

1. User creates TPOT instance
2. User calls TPOT `fit()` function with training data
3. TPOT crashes with a `KeyError` after 5 generations

## Expected result

[describe what you would expect to have resulted from this process]

## Current result

[describe what you currently experience from this process, and thereby explain the bug]

## Possible fix

[not necessary, but suggest fixes or reasons for the bug]

## `name of issue` screenshot

[if relevant, include a screenshot]
Master status: [![Master Build Status - Mac/Linux](https://travis-ci.com/EpistasisLab/tpot.svg?branch=master)](https://travis-ci.com/EpistasisLab/tpot)
[![Master Build Status - Windows](https://ci.appveyor.com/api/projects/status/b7bmpwpkjhifrm7v/branch/master?svg=true)](https://ci.appveyor.com/project/weixuanfu/tpot?branch=master)
[![Master Coverage Status](https://coveralls.io/repos/github/EpistasisLab/tpot/badge.svg?branch=master)](https://coveralls.io/github/EpistasisLab/tpot?branch=master)

Development status: [![Development Build Status - Mac/Linux](https://travis-ci.com/EpistasisLab/tpot.svg?branch=development)](https://travis-ci.com/EpistasisLab/tpot/branches)
[![Development Build Status - Windows](https://ci.appveyor.com/api/projects/status/b7bmpwpkjhifrm7v/branch/development?svg=true)](https://ci.appveyor.com/project/weixuanfu/tpot?branch=development)
[![Development Coverage Status](https://coveralls.io/repos/github/EpistasisLab/tpot/badge.svg?branch=development)](https://coveralls.io/github/EpistasisLab/tpot?branch=development)

Package information: [![Python 3.7](https://img.shields.io/badge/python-3.7-blue.svg)](https://www.python.org/downloads/release/python-370/)
[![License: LGPL v3](https://img.shields.io/badge/license-LGPL%20v3-blue.svg)](http://www.gnu.org/licenses/lgpl-3.0)
[![PyPI version](https://badge.fury.io/py/TPOT.svg)](https://badge.fury.io/py/TPOT)

<p align="center">
<img src="https://raw.githubusercontent.com/EpistasisLab/tpot/master/images/tpot-logo.jpg" width=300 />
</p>

**TPOT** stands for **T**ree-based **P**ipeline **O**ptimization **T**ool. Consider TPOT your **Data Science Assistant**. TPOT is a Python Automated Machine Learning tool that optimizes machine learning pipelines using genetic programming.

![TPOT Demo](https://github.com/EpistasisLab/tpot/blob/master/images/tpot-demo.gif "TPOT Demo")

TPOT will automate the most tedious part of machine learning by intelligently exploring thousands of possible pipelines to find the best one for your data.

![An example Machine Learning pipeline](https://github.com/EpistasisLab/tpot/blob/master/images/tpot-ml-pipeline.png "An example Machine Learning pipeline")

<p align="center"><strong>An example Machine Learning pipeline</strong></p>

Once TPOT is finished searching (or you get tired of waiting), it provides you with the Python code for the best pipeline it found so you can tinker with the pipeline from there.

![An example TPOT pipeline](https://github.com/EpistasisLab/tpot/blob/master/images/tpot-pipeline-example.png "An example TPOT pipeline")

TPOT is built on top of scikit-learn, so all of the code it generates should look familiar... if you're familiar with scikit-learn, anyway.

**TPOT is still under active development** and we encourage you to check back on this repository regularly for updates.

For further information about TPOT, please see the [project documentation](http://epistasislab.github.io/tpot/).

## License

Please see the [repository license](https://github.com/EpistasisLab/tpot/blob/master/LICENSE) for the licensing and usage information for TPOT.

Generally, we have licensed TPOT to make it as widely usable as possible.

## Installation

We maintain the [TPOT installation instructions](http://epistasislab.github.io/tpot/installing/) in the documentation. TPOT requires a working installation of Python.

## Usage

TPOT can be used [on the command line](http://epistasislab.github.io/tpot/using/#tpot-on-the-command-line) or [with Python code](http://epistasislab.github.io/tpot/using/#tpot-with-code).

Click on the corresponding links to find more information on TPOT usage in the documentation.

## Examples

### Classification

Below is a minimal working example with the optical recognition of handwritten digits dataset.

```python
from tpot import TPOTClassifier
from sklearn.datasets import load_digits
from sklearn.model_selection import train_test_split

digits = load_digits()
X_train, X_test, y_train, y_test = train_test_split(digits.data, digits.target,
                                                    train_size=0.75, test_size=0.25, random_state=42)

tpot = TPOTClassifier(generations=5, population_size=50, verbosity=2, random_state=42)
tpot.fit(X_train, y_train)
print(tpot.score(X_test, y_test))
tpot.export('tpot_digits_pipeline.py')
```

Running this code should discover a pipeline that achieves about 98% testing accuracy, and the corresponding Python code should be exported to the `tpot_digits_pipeline.py` file and look similar to the following:

```python
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline, make_union
from sklearn.preprocessing import PolynomialFeatures
from tpot.builtins import StackingEstimator
from tpot.export_utils import set_param_recursive

# NOTE: Make sure that the outcome column is labeled 'target' in the data file
tpot_data = pd.read_csv('PATH/TO/DATA/FILE', sep='COLUMN_SEPARATOR', dtype=np.float64)
features = tpot_data.drop('target', axis=1)
training_features, testing_features, training_target, testing_target = \
            train_test_split(features, tpot_data['target'], random_state=42)

# Average CV score on the training set was: 0.9799428471757372
exported_pipeline = make_pipeline(
    PolynomialFeatures(degree=2, include_bias=False, interaction_only=False),
    StackingEstimator(estimator=LogisticRegression(C=0.1, dual=False, penalty="l1")),
    RandomForestClassifier(bootstrap=True, criterion="entropy", max_features=0.35000000000000003, min_samples_leaf=20, min_samples_split=19, n_estimators=100)
)
# Fix random state for all the steps in exported pipeline
set_param_recursive(exported_pipeline.steps, 'random_state', 42)

exported_pipeline.fit(training_features, training_target)
results = exported_pipeline.predict(testing_features)
```

### Regression

Similarly, TPOT can optimize pipelines for regression problems. Below is a minimal working example with the practice Boston housing prices data set.

```python
from tpot import TPOTRegressor
from sklearn.datasets import load_boston
from sklearn.model_selection import train_test_split

housing = load_boston()
X_train, X_test, y_train, y_test = train_test_split(housing.data, housing.target,
                                                    train_size=0.75, test_size=0.25, random_state=42)

tpot = TPOTRegressor(generations=5, population_size=50, verbosity=2, random_state=42)
tpot.fit(X_train, y_train)
print(tpot.score(X_test, y_test))
tpot.export('tpot_boston_pipeline.py')
```

which should result in a pipeline that achieves about 12.77 mean squared error (MSE), and the Python code in `tpot_boston_pipeline.py` should look similar to:

```python
import numpy as np
import pandas as pd
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import PolynomialFeatures
from tpot.export_utils import set_param_recursive

# NOTE: Make sure that the outcome column is labeled 'target' in the data file
tpot_data = pd.read_csv('PATH/TO/DATA/FILE', sep='COLUMN_SEPARATOR', dtype=np.float64)
features = tpot_data.drop('target', axis=1)
training_features, testing_features, training_target, testing_target = \
            train_test_split(features, tpot_data['target'], random_state=42)

# Average CV score on the training set was: -10.812040755234403
exported_pipeline = make_pipeline(
    PolynomialFeatures(degree=2, include_bias=False, interaction_only=False),
    ExtraTreesRegressor(bootstrap=False, max_features=0.5, min_samples_leaf=2, min_samples_split=3, n_estimators=100)
)
# Fix random state for all the steps in exported pipeline
set_param_recursive(exported_pipeline.steps, 'random_state', 42)

exported_pipeline.fit(training_features, training_target)
results = exported_pipeline.predict(testing_features)
```

Check the documentation for [more examples and tutorials](http://epistasislab.github.io/tpot/examples/).

## Contributing to TPOT

We welcome you to [check the existing issues](https://github.com/EpistasisLab/tpot/issues/) for bugs or enhancements to work on. If you have an idea for an extension to TPOT, please [file a new issue](https://github.com/EpistasisLab/tpot/issues/new) so we can discuss it.

Before submitting any contributions, please review our [contribution guidelines](http://epistasislab.github.io/tpot/contributing/).

## Having problems or have questions about TPOT?

Please [check the existing open and closed issues](https://github.com/EpistasisLab/tpot/issues?utf8=%E2%9C%93&q=is%3Aissue) to see if your issue has already been attended to. If it hasn't, [file a new issue](https://github.com/EpistasisLab/tpot/issues/new) on this repository so we can review your issue.

## Citing TPOT

If you use TPOT in a scientific publication, please consider citing at least one of the following papers:

Trang T. Le, Weixuan Fu and Jason H. Moore (2020). [Scaling tree-based automated machine learning to biomedical big data with a feature set selector](https://academic.oup.com/bioinformatics/article/36/1/250/5511404). *Bioinformatics*.36(1): 250-256.

BibTeX entry:

```bibtex
@article{le2020scaling,
  title={Scaling tree-based automated machine learning to biomedical big data with a feature set selector},
  author={Le, Trang T and Fu, Weixuan and Moore, Jason H},
  journal={Bioinformatics},
  volume={36},
  number={1},
  pages={250--256},
  year={2020},
  publisher={Oxford University Press}
}
```


Randal S. Olson, Ryan J. Urbanowicz, Peter C. Andrews, Nicole A. Lavender, La Creis Kidd, and Jason H. Moore (2016). [Automating biomedical data science through tree-based pipeline optimization](http://link.springer.com/chapter/10.1007/978-3-319-31204-0_9). *Applications of Evolutionary Computation*, pages 123-137.

BibTeX entry:

```bibtex
@inbook{Olson2016EvoBio,
    author={Olson, Randal S. and Urbanowicz, Ryan J. and Andrews, Peter C. and Lavender, Nicole A. and Kidd, La Creis and Moore, Jason H.},
    editor={Squillero, Giovanni and Burelli, Paolo},
    chapter={Automating Biomedical Data Science Through Tree-Based Pipeline Optimization},
    title={Applications of Evolutionary Computation: 19th European Conference, EvoApplications 2016, Porto, Portugal, March 30 -- April 1, 2016, Proceedings, Part I},
    year={2016},
    publisher={Springer International Publishing},
    pages={123--137},
    isbn={978-3-319-31204-0},
    doi={10.1007/978-3-319-31204-0_9},
    url={http://dx.doi.org/10.1007/978-3-319-31204-0_9}
}
```

Randal S. Olson, Nathan Bartley, Ryan J. Urbanowicz, and Jason H. Moore (2016). [Evaluation of a Tree-based Pipeline Optimization Tool for Automating Data Science](http://dl.acm.org/citation.cfm?id=2908918). *Proceedings of GECCO 2016*, pages 485-492.

BibTeX entry:

```bibtex
@inproceedings{OlsonGECCO2016,
    author = {Olson, Randal S. and Bartley, Nathan and Urbanowicz, Ryan J. and Moore, Jason H.},
    title = {Evaluation of a Tree-based Pipeline Optimization Tool for Automating Data Science},
    booktitle = {Proceedings of the Genetic and Evolutionary Computation Conference 2016},
    series = {GECCO '16},
    year = {2016},
    isbn = {978-1-4503-4206-3},
    location = {Denver, Colorado, USA},
    pages = {485--492},
    numpages = {8},
    url = {http://doi.acm.org/10.1145/2908812.2908918},
    doi = {10.1145/2908812.2908918},
    acmid = {2908918},
    publisher = {ACM},
    address = {New York, NY, USA},
}
```

Alternatively, you can cite the repository directly with the following DOI:

[![DOI](https://zenodo.org/badge/20747/rhiever/tpot.svg)](https://zenodo.org/badge/latestdoi/20747/rhiever/tpot)

## Support for TPOT

TPOT was developed in the [Computational Genetics Lab](http://epistasis.org/) at the [University of Pennsylvania](https://www.upenn.edu/) with funding from the [NIH](http://www.nih.gov/) under grant R01 AI117694. We are incredibly grateful for the support of the NIH and the University of Pennsylvania during the development of this project.

The TPOT logo was designed by Todd Newmuis, who generously donated his time to the project.
[please review the [Contribution Guidelines](http://epistasislab.github.io/tpot/contributing/) prior to submitting your pull request. go ahead and delete this line if you've already reviewed said guidelines.]

## What does this PR do?



## Where should the reviewer start?



## How should this PR be tested?



## Any background context you want to provide?



## What are the relevant issues?

[you can link directly to issues by entering # then the number of the issue, for example, #3 links to issue 3]

## Screenshots (if appropriate)



## Questions:

- Do the docs need to be updated?
- Does this PR add new (Python) dependencies?
TPOT was developed in the [Computational Genetics Lab](http://epistasis.org/) at the [University of Pennsylvania](https://www.upenn.edu/) with funding from the [NIH](http://www.nih.gov/) under grant R01 AI117694. We are incredibly grateful for the support of the NIH and the University of Pennsylvania during the development of this project.

The TPOT logo was designed by Todd Newmuis, who generously donated his time to the project.
<center>
<img src="https://raw.githubusercontent.com/EpistasisLab/tpot/master/images/tpot-logo.jpg" width=300 />
</center>

Consider TPOT your **Data Science Assistant**. TPOT is a Python Automated Machine Learning tool that optimizes machine learning pipelines using genetic programming.

<br />

<center>
<img src="https://raw.githubusercontent.com/EpistasisLab/tpot/master/images/tpot-demo.gif" width=800 alt="TPOT Demo" />
</center>

<br />

TPOT will automate the most tedious part of machine learning by intelligently exploring thousands of possible pipelines to find the best one for your data.

<br />

<center>
<img src="https://raw.githubusercontent.com/EpistasisLab/tpot/master/images/tpot-ml-pipeline.png" width=800 alt="An example machine learning pipeline" />

<strong>An example machine learning pipeline</strong>
</center>

<br />

Once TPOT is finished searching (or you get tired of waiting), it provides you with the Python code for the best pipeline it found so you can tinker with the pipeline from there.

<br />

<center>
<img src="https://raw.githubusercontent.com/EpistasisLab/tpot/master/images/tpot-pipeline-example.png" width=800 alt="An example TPOT pipeline" />

<strong>An example TPOT pipeline</strong>
</center>

<br />

TPOT is built on top of scikit-learn, so all of the code it generates should look familiar... if you're familiar with scikit-learn, anyway.

**TPOT is still under active development** and we encourage you to check back on this repository regularly for updates.
# Citing TPOT

If you use TPOT in a scientific publication, please consider citing at least one of the following papers:


Trang T. Le, Weixuan Fu and Jason H. Moore (2020). [Scaling tree-based automated machine learning to biomedical big data with a feature set selector](https://academic.oup.com/bioinformatics/article/36/1/250/5511404). *Bioinformatics*.36(1): 250-256.

BibTeX entry:

```bibtex
@article{le2020scaling,
  title={Scaling tree-based automated machine learning to biomedical big data with a feature set selector},
  author={Le, Trang T and Fu, Weixuan and Moore, Jason H},
  journal={Bioinformatics},
  volume={36},
  number={1},
  pages={250--256},
  year={2020},
  publisher={Oxford University Press}
}
```



Randal S. Olson, Ryan J. Urbanowicz, Peter C. Andrews, Nicole A. Lavender, La Creis Kidd, and Jason H. Moore (2016). [Automating biomedical data science through tree-based pipeline optimization](http://link.springer.com/chapter/10.1007/978-3-319-31204-0_9). *Applications of Evolutionary Computation*, pages 123-137.

BibTeX entry:

```bibtex
@inbook{Olson2016EvoBio,
    author={Olson, Randal S. and Urbanowicz, Ryan J. and Andrews, Peter C. and Lavender, Nicole A. and Kidd, La Creis and Moore, Jason H.},
    editor={Squillero, Giovanni and Burelli, Paolo},
    chapter={Automating Biomedical Data Science Through Tree-Based Pipeline Optimization},
    title={Applications of Evolutionary Computation: 19th European Conference, EvoApplications 2016, Porto, Portugal, March 30 -- April 1, 2016, Proceedings, Part I},
    year={2016},
    publisher={Springer International Publishing},
    pages={123--137},
    isbn={978-3-319-31204-0},
    doi={10.1007/978-3-319-31204-0_9},
    url={http://dx.doi.org/10.1007/978-3-319-31204-0_9}
}
```

Evaluation of a Tree-based Pipeline Optimization Tool for Automating Data Science

Randal S. Olson, Nathan Bartley, Ryan J. Urbanowicz, and Jason H. Moore (2016). [Evaluation of a Tree-based Pipeline Optimization Tool for Automating Data Science](http://dl.acm.org/citation.cfm?id=2908918). *Proceedings of GECCO 2016*, pages 485-492.

BibTeX entry:

```bibtex
@inproceedings{OlsonGECCO2016,
    author = {Olson, Randal S. and Bartley, Nathan and Urbanowicz, Ryan J. and Moore, Jason H.},
    title = {Evaluation of a Tree-based Pipeline Optimization Tool for Automating Data Science},
    booktitle = {Proceedings of the Genetic and Evolutionary Computation Conference 2016},
    series = {GECCO '16},
    year = {2016},
    isbn = {978-1-4503-4206-3},
    location = {Denver, Colorado, USA},
    pages = {485--492},
    numpages = {8},
    url = {http://doi.acm.org/10.1145/2908812.2908918},
    doi = {10.1145/2908812.2908918},
    acmid = {2908918},
    publisher = {ACM},
    address = {New York, NY, USA},
}
```

Alternatively, you can cite the repository directly with the following DOI:

[DOI](https://zenodo.org/badge/latestdoi/20747/rhiever/tpot)
# TPOT API

## Classification

<pre><em>class</em> tpot.<strong style="color:#008AB8">TPOTClassifier</strong>(<em><strong>generations</strong>=100, <strong>population_size</strong>=100,
                          <strong>offspring_size</strong>=None, <strong>mutation_rate</strong>=0.9,
                          <strong>crossover_rate</strong>=0.1,
                          <strong>scoring</strong>='accuracy', <strong>cv</strong>=5,
                          <strong>subsample</strong>=1.0, <strong>n_jobs</strong>=1,
                          <strong>max_time_mins</strong>=None, <strong>max_eval_time_mins</strong>=5,
                          <strong>random_state</strong>=None, <strong>config_dict</strong>=None,
                          <strong>template</strong>=None,
                          <strong>warm_start</strong>=False,
                          <strong>memory</strong>=None,
                          <strong>use_dask</strong>=False,
                          <strong>periodic_checkpoint_folder</strong>=None,
                          <strong>early_stop</strong>=None,
                          <strong>verbosity</strong>=0,
                          <strong>disable_update_check</strong>=False,
                          <strong>log_file</strong>=None
                          </em>)</pre>
<div align="right"><a href="https://github.com/EpistasisLab/tpot/blob/master/tpot/base.py">source</a></div>

Automated machine learning for supervised classification tasks.

The TPOTClassifier performs an intelligent search over machine learning pipelines that can contain supervised classification models,
preprocessors, feature selection techniques, and any other estimator or transformer that follows the [scikit-learn API](http://scikit-learn.org/stable/developers/contributing.html#apis-of-scikit-learn-objects).
The TPOTClassifier will also search over the hyperparameters of all objects in the pipeline.

By default, TPOTClassifier will search over a broad range of supervised classification algorithms, transformers, and their parameters.
However, the algorithms, transformers, and hyperparameters that the TPOTClassifier searches over can be fully customized using the `config_dict` parameter.

Read more in the [User Guide](using/#tpot-with-code).

<table>
<tr>
<td width="20%" style="vertical-align:top; background:#F5F5F5;"><strong>Parameters:</strong></td>
<td width="80%" style="background:white;">
<strong>generations</strong>: int or None optional (default=100)
<blockquote>
Number of iterations to the run pipeline optimization process. It must be a positive number or None. If None, the parameter <em>max_time_mins</em> must be defined as the runtime limit.
<br /><br />
Generally, TPOT will work better when you give it more generations (and therefore time) to optimize the pipeline.
<br /><br />
TPOT will evaluate <em>population_size</em> + <em>generations</em> × <em>offspring_size</em> pipelines in total.
</blockquote>

<strong>population_size</strong>: int, optional (default=100)
<blockquote>
Number of individuals to retain in the genetic programming population every generation. Must be a positive number.
<br /><br />
Generally, TPOT will work better when you give it more individuals with which to optimize the pipeline.
</blockquote>

<strong>offspring_size</strong>: int, optional (default=None)
<blockquote>
Number of offspring to produce in each genetic programming generation. Must be a positive number. By default, the number of offspring is equal to the number of population size.
</blockquote>

<strong>mutation_rate</strong>: float, optional (default=0.9)
<blockquote>
Mutation rate for the genetic programming algorithm in the range [0.0, 1.0]. This parameter tells the GP algorithm how many pipelines to apply random changes to every generation.
<br /><br />
<em>mutation_rate</em> + <em>crossover_rate</em> cannot exceed 1.0.
<br /><br />
We recommend using the default parameter unless you understand how the mutation rate affects GP algorithms.
</blockquote>

<strong>crossover_rate</strong>: float, optional (default=0.1)
<blockquote>
Crossover rate for the genetic programming algorithm in the range [0.0, 1.0]. This parameter tells the genetic programming algorithm how many pipelines to "breed" every generation.
<br /><br />
<em>mutation_rate</em> + <em>crossover_rate</em> cannot exceed 1.0.
<br /><br />
We recommend using the default parameter unless you understand how the crossover rate affects GP algorithms.
</blockquote>

<strong>scoring</strong>: string or callable, optional (default='accuracy')
<blockquote>
Function used to evaluate the quality of a given pipeline for the classification problem. The following built-in scoring functions can be used:
<br /><br/>
'accuracy', 'adjusted_rand_score', 'average_precision', 'balanced_accuracy', 'f1', 'f1_macro', 'f1_micro', 'f1_samples', 'f1_weighted', 'neg_log_loss', 'precision' etc. (suffixes apply as with ‘f1’), 'recall' etc. (suffixes apply as with ‘f1’), ‘jaccard’ etc. (suffixes apply as with ‘f1’), 'roc_auc', ‘roc_auc_ovr’, ‘roc_auc_ovo’, ‘roc_auc_ovr_weighted’, ‘roc_auc_ovo_weighted’
<br /><br/>
If you would like to use a custom scorer, you can pass the callable object/function with signature <em>scorer(estimator, X, y)</em>.
<br /><br/>
See the section on <a href="../using/#scoring-functions">scoring functions</a> for more details.

</blockquote>

<strong>cv</strong>: int, cross-validation generator, or an iterable, optional (default=5)
<blockquote>
Cross-validation strategy used when evaluating pipelines.
<br /><br />
Possible inputs:
<ul>
<li>integer, to specify the number of folds in a StratifiedKFold,</li>
<li>An object to be used as a cross-validation generator, or</li>
<li>An iterable yielding train/test splits.</li>
</blockquote>

<strong>subsample</strong>: float, optional (default=1.0)
<blockquote>
Fraction of training samples that are used during the TPOT optimization process. Must be in the range (0.0, 1.0].
<br /><br />
Setting <em>subsample</em>=0.5 tells TPOT to use a random subsample of half of the training data. This subsample will remain the same during the entire pipeline optimization process.
</blockquote>

<strong>n_jobs</strong>: integer, optional (default=1)
<blockquote>
Number of processes to use in parallel for evaluating pipelines during the TPOT optimization process.
<br /><br />
Setting <em>n_jobs</em>=-1 will use as many cores as available on the computer. For n_jobs below -1, (n_cpus + 1 + n_jobs) are used. Thus for n_jobs = -2, all CPUs but one are used. Beware that using multiple processes on the same machine may cause memory issues for large datasets.
</blockquote>

<strong>max_time_mins</strong>: integer or None, optional (default=None)
<blockquote>
How many minutes TPOT has to optimize the pipeline.
<br /><br />
If not None, this setting will allow TPOT to run until <em>max_time_mins</em> minutes elapsed and then stop. TPOT will stop earlier if <em>generations</em> is set and all generations are already evaluated.
</blockquote>

<strong>max_eval_time_mins</strong>: float, optional (default=5)
<blockquote>
How many minutes TPOT has to evaluate a single pipeline.
<br /><br />
Setting this parameter to higher values will allow TPOT to evaluate more complex pipelines, but will also allow TPOT to run longer. Use this parameter to help prevent TPOT from wasting time on evaluating time-consuming pipelines.
</blockquote>

<strong>random_state</strong>: integer or None, optional (default=None)
<blockquote>
The seed of the pseudo random number generator used in TPOT.
<br /><br />
Use this parameter to make sure that TPOT will give you the same results each time you run it against the same data set with that seed.
</blockquote>

<strong>config_dict</strong>: Python dictionary, string, or None, optional (default=None)
<blockquote>
A configuration dictionary for customizing the operators and parameters that TPOT searches in the optimization process.
<br /><br />
Possible inputs are:
<ul>
<li>Python dictionary, TPOT will use your custom configuration,</li>
<li>string 'TPOT light', TPOT will use a built-in configuration with only fast models and preprocessors, or</li>
<li>string 'TPOT MDR', TPOT will use a built-in configuration specialized for genomic studies, or</li>
<li>string 'TPOT sparse': TPOT will use a configuration dictionary with a one-hot encoder and the operators normally included in TPOT that also support sparse matrices, or</li>
<li>None, TPOT will use the default TPOTClassifier configuration.</li>
</ul>
See the <a href="../using/#built-in-tpot-configurations">built-in configurations</a> section for the list of configurations included with TPOT, and the <a href="../using/#customizing-tpots-operators-and-parameters">custom configuration</a> section for more information and examples of how to create your own TPOT configurations.
</blockquote>

<strong>template</strong>: string (default=None)
<blockquote>
Template of predefined pipeline structure. The option is for specifying a desired structure for the machine learning pipeline evaluated in TPOT.
<br /><br />

So far this option only supports linear pipeline structure. Each step in the pipeline should be a main class of operators (Selector, Transformer, Classifier) or a specific operator (e.g. `SelectPercentile`) defined in TPOT operator configuration. If one step is a main class, TPOT will randomly assign all subclass operators (subclasses of [`SelectorMixin`](https://github.com/scikit-learn/scikit-learn/blob/master/sklearn/feature_selection/base.py#L17), [`TransformerMixin`](https://scikit-learn.org/stable/modules/generated/sklearn.base.TransformerMixin.html), [`ClassifierMixin`](https://scikit-learn.org/stable/modules/generated/sklearn.base.ClassifierMixin.html) in scikit-learn) to that step. Steps in the template are delimited by "-", e.g. "SelectPercentile-Transformer-Classifier". By default value of template is None, TPOT generates tree-based pipeline randomly.

See the <a href="../using/#template-option-in-tpot"> template option in tpot</a> section for more details.
</blockquote>

<strong>warm_start</strong>: boolean, optional (default=False)
<blockquote>
Flag indicating whether the TPOT instance will reuse the population from previous calls to <em>fit()</em>.
<br /><br />
Setting <em>warm_start</em>=True can be useful for running TPOT for a short time on a dataset, checking the results, then resuming the TPOT run from where it left off.
</blockquote>

<strong>memory</strong>: a joblib.Memory object or string, optional (default=None)
<blockquote>
If supplied, pipeline will cache each transformer after calling fit. This feature is used to avoid computing the fit transformers within a pipeline if the parameters and input data are identical with another fitted pipeline during optimization process. More details about memory caching in <a href="http://scikit-learn.org/stable/modules/pipeline.html#caching-transformers-avoid-repeated-computation">scikit-learn documentation</a>
<br /><br />
Possible inputs are:
<ul>
<li>String 'auto': TPOT uses memory caching with a temporary directory and cleans it up upon shutdown, or</li>
<li>Path of a caching directory, TPOT uses memory caching with the provided directory and TPOT does NOT clean the caching directory up upon shutdown, or</li>
<li>Memory object, TPOT uses the instance of joblib.Memory for memory caching and TPOT does NOT clean the caching directory up upon shutdown, or</li>
<li>None, TPOT does not use memory caching.</li>
</ul>
</blockquote>

<strong>use_dask</strong>: boolean, optional (default: False)
<blockquote>
Whether to use Dask-ML's pipeline optimiziations. This avoid re-fitting
the same estimator on the same split of data multiple times. It
will also provide more detailed diagnostics when using Dask's
distributed scheduler.
<br /><br />
See <a href="https://dask-ml.readthedocs.io/en/latest/hyper-parameter-search.html#avoid-repeated-work">avoid repeated work</a> for more details.
</blockquote>

<strong>periodic_checkpoint_folder</strong>: path string, optional (default: None)
<blockquote>
If supplied, a folder in which TPOT will periodically save pipelines in pareto front so far while optimizing.<br /><br />
Currently once per generation but not more often than once per 30 seconds.<br /><br />
Useful in multiple cases:
<ul>
<li>Sudden death before TPOT could save optimized pipeline</li>
<li>Track its progress</li>
<li>Grab pipelines while it's still optimizing</li>
</ul>
</blockquote>

<strong>early_stop</strong>: integer, optional (default: None)
<blockquote>
How many generations TPOT checks whether there is no improvement in optimization process.
<br /><br />
Ends the optimization process if there is no improvement in the given number of generations.
</blockquote>

<strong>verbosity</strong>: integer, optional (default=0)
<blockquote>
How much information TPOT communicates while it's running.
<br /><br />
Possible inputs are:
<ul>
<li>0, TPOT will print nothing,</li>
<li>1, TPOT will print minimal information,</li>
<li>2, TPOT will print more information and provide a progress bar, or</li>
<li>3, TPOT will print everything and provide a progress bar.</li>
</ul>
</blockquote>

<strong>disable_update_check</strong>: boolean, optional (default=False)
<blockquote>
Flag indicating whether the TPOT version checker should be disabled.
<br /><br />
The update checker will tell you when a new version of TPOT has been released.
</blockquote>

<strong>log_file</strong>: file-like class (io.TextIOWrapper or io.StringIO) or string, optional (default: None)
<br /><br />
<blockquote>
Save progress content to a file.
If it is a string for the path and file name of the desired output file,
TPOT will create the file and write log into it.
If it is None, TPOT will output log into sys.stdout
</blockquote>

</td>
</tr>

<tr>
<td width="20%" style="vertical-align:top; background:#F5F5F5;"><strong>Attributes:</strong></td>
<td width="80%" style="background:white;">
<strong>fitted_pipeline_</strong>: scikit-learn Pipeline object
<blockquote>
The best pipeline that TPOT discovered during the pipeline optimization process, fitted on the entire training dataset.
</blockquote>

<strong>pareto_front_fitted_pipelines_</strong>: Python dictionary
<blockquote>
Dictionary containing the all pipelines on the TPOT Pareto front, where the key is the string representation of the pipeline and the value is the corresponding pipeline fitted on the entire training dataset.
<br /><br />
The TPOT Pareto front provides a trade-off between pipeline complexity (i.e., the number of steps in the pipeline) and the predictive performance of the pipeline.
<br /><br />
Note: <em>pareto_front_fitted_pipelines_</em> is only available when <em>verbosity</em>=3.
</blockquote>

<strong>evaluated_individuals_</strong>: Python dictionary
<blockquote>
Dictionary containing all pipelines that were evaluated during the pipeline optimization process, where the key is the string representation of the pipeline and the value is a tuple containing (# of steps in pipeline, accuracy metric for the pipeline).
<br /><br />
This attribute is primarily for internal use, but may be useful for looking at the other pipelines that TPOT evaluated.
</blockquote>
</td>
<tr>
</table>

<strong>Example</strong>

```Python
from tpot import TPOTClassifier
from sklearn.datasets import load_digits
from sklearn.model_selection import train_test_split

digits = load_digits()
X_train, X_test, y_train, y_test = train_test_split(digits.data, digits.target,
                                                    train_size=0.75, test_size=0.25)

tpot = TPOTClassifier(generations=5, population_size=50, verbosity=2)
tpot.fit(X_train, y_train)
print(tpot.score(X_test, y_test))
tpot.export('tpot_digits_pipeline.py')
```

<strong>Functions</strong>

<table width="100%">
<tr>
<td width="25%"><a href="#tpotclassifier-fit">fit</a>(features, classes[, sample_weight, groups])</td>
<td>Run the TPOT optimization process on the given training data.</td>
</tr>

<tr>
<td><a href="#tpotclassifier-predict">predict</a>(features)</td>
<td>Use the optimized pipeline to predict the classes for a feature set.</td>
</tr>

<tr>
<td><a href="#tpotclassifier-predict-proba">predict_proba</a>(features)</td>
<td>Use the optimized pipeline to estimate the class probabilities for a feature set.</td>
</tr>

<tr>
<td><a href="#tpotclassifier-score">score</a>(testing_features, testing_classes)</td>
<td>Returns the optimized pipeline's score on the given testing data using the user-specified scoring function.</td>
</tr>

<tr>
<td><a href="#tpotclassifier-export">export</a>(output_file_name)</td>
<td>Export the optimized pipeline as Python code.</td>
</tr>
</table>


<a name="tpotclassifier-fit"></a>
```Python
fit(features, classes, sample_weight=None, groups=None)
```

<div style="padding-left:5%" width="100%">
Run the TPOT optimization process on the given training data.
<br /><br />
Uses genetic programming to optimize a machine learning pipeline that maximizes the score on the provided features and target. This pipeline optimization procedure uses internal k-fold cross-validaton to avoid overfitting on the provided data. At the end of the pipeline optimization procedure, the best pipeline is then trained on the entire set of provided samples.
<br /><br />
<table width="100%">
<tr>
<td width="20%" style="vertical-align:top; background:#F5F5F5;"><strong>Parameters:</strong></td>
<td width="80%" style="background:white;">
<strong>features</strong>: array-like {n_samples, n_features}
<blockquote>
Feature matrix
<br /><br />
TPOT and all scikit-learn algorithms assume that the features will be numerical and there will be no missing values.
As such, when a feature matrix is provided to TPOT, all missing values will automatically be replaced (i.e., imputed)
using <a href="http://scikit-learn.org/stable/modules/generated/sklearn.impute.SimpleImputer.html">median value imputation</a>.
<br /><br />
If you wish to use a different imputation strategy than median imputation, please make sure to apply imputation to your feature set prior to passing it to TPOT.
</blockquote>

<strong>classes</strong>: array-like {n_samples}
<blockquote>
List of class labels for prediction
</blockquote>

<strong>sample_weight</strong>: array-like {n_samples}, optional
<blockquote>
Per-sample weights. Higher weights indicate more importance. If specified, sample_weight will be passed to any pipeline element whose fit() function accepts a sample_weight argument. By default, using sample_weight does not affect tpot's scoring functions, which determine preferences between pipelines.
</blockquote>

<strong>groups</strong>: array-like, with shape {n_samples, }, optional
<blockquote>
Group labels for the samples used when performing cross-validation.
<br /><br />
This parameter should only be used in conjunction with sklearn's Group cross-validation functions, such as <a href="http://scikit-learn.org/stable/modules/generated/sklearn.model_selection.GroupKFold.html">sklearn.model_selection.GroupKFold</a>.
</blockquote>
</td>
</tr>
<tr>
<td width="20%" style="vertical-align:top; background:#F5F5F5;"><strong>Returns:</strong></td>
<td width="80%" style="background:white;">
<strong>self</strong>: object
<blockquote>
Returns a copy of the fitted TPOT object
</blockquote>
</td>
</tr>
</table>
</div>


<a name="tpotclassifier-predict"></a>
```Python
predict(features)
```

<div style="padding-left:5%" width="100%">
Use the optimized pipeline to predict the classes for a feature set.
<br /><br />
<table width="100%">
<tr>
<td width="20%" style="vertical-align:top; background:#F5F5F5;"><strong>Parameters:</strong></td>
<td width="80%" style="background:white;">
<strong>features</strong>: array-like {n_samples, n_features}
<blockquote>
Feature matrix
</blockquote>
</td>
</tr>
<tr>
<td width="20%" style="vertical-align:top; background:#F5F5F5;"><strong>Returns:</strong></td>
<td width="80%" style="background:white;">
<strong>predictions</strong>: array-like {n_samples}
<blockquote>
Predicted classes for the samples in the feature matrix
</blockquote>
</td>
</tr>
</table>
</div>


<a name="tpotclassifier-predict-proba"></a>
```Python
predict_proba(features)
```

<div style="padding-left:5%" width="100%">
Use the optimized pipeline to estimate the class probabilities for a feature set.
<br /><br />
Note: This function will only work for pipelines whose final classifier supports the <em>predict_proba</em> function. TPOT will raise an error otherwise.
<br /><br />
<table width="100%">
<tr>
<td width="20%" style="vertical-align:top; background:#F5F5F5;"><strong>Parameters:</strong></td>
<td width="80%" style="background:white;">
<strong>features</strong>: array-like {n_samples, n_features}
<blockquote>
Feature matrix
</blockquote>
</td>
</tr>
<tr>
<td width="20%" style="vertical-align:top; background:#F5F5F5;"><strong>Returns:</strong></td>
<td width="80%" style="background:white;">
<strong>predictions</strong>: array-like {n_samples, n_classes}
<blockquote>
The class probabilities of the input samples
</blockquote>
</td>
</tr>
</table>
</div>


<a name="tpotclassifier-score"></a>
```Python
score(testing_features, testing_classes)
```

<div style="padding-left:5%" width="100%">
Returns the optimized pipeline's score on the given testing data using the user-specified scoring function.
<br /><br />
The default scoring function for TPOTClassifier is 'accuracy'.
<br /><br />
<table width="100%">
<tr>
<td width="20%" style="vertical-align:top; background:#F5F5F5;"><strong>Parameters:</strong></td>
<td width="80%" style="background:white;">
<strong>testing_features</strong>: array-like {n_samples, n_features}
<blockquote>
Feature matrix of the testing set
</blockquote>

<strong>testing_classes</strong>: array-like {n_samples}
<blockquote>
List of class labels for prediction in the testing set
</blockquote>
</td>
</tr>
<tr>
<td width="20%" style="vertical-align:top; background:#F5F5F5;"><strong>Returns:</strong></td>
<td width="80%" style="background:white;">
<strong>accuracy_score</strong>: float
<blockquote>
The estimated test set accuracy according to the user-specified scoring function.
</blockquote>
</td>
</tr>
</table>
</div>


<a name="tpotclassifier-export"></a>
```Python
export(output_file_name, data_file_path)
```

<div style="padding-left:5%" width="100%">
Export the optimized pipeline as Python code.
<br /><br />
See the <a href="../using/#tpot-with-code">usage documentation</a> for example usage of the export function.
<br /><br />
<table width="100%">
<tr>
<td width="20%" style="vertical-align:top; background:#F5F5F5;"><strong>Parameters:</strong></td>
<td width="80%" style="background:white;">
<strong>output_file_name</strong>: string
<blockquote>
String containing the path and file name of the desired output file
</blockquote>
<strong>data_file_path</strong>: string
<blockquote>
By default, the path of input dataset is 'PATH/TO/DATA/FILE' by default. If data_file_path is another string, the path will be replaced.
</blockquote>
</tr>
<tr>
<td width="20%" style="vertical-align:top; background:#F5F5F5;"><strong>Returns:</strong></td>
<td width="80%" style="background:white;">
<strong>exported_code_string</strong>: string
<blockquote>
The whole pipeline text as a string should be returned if output_file_name is not specified.
</blockquote>
</td>
</tr>
</table>
</div>




## Regression

<pre><em>class</em> tpot.<strong style="color:#008AB8">TPOTRegressor</strong>(<em><strong>generations</strong>=100, <strong>population_size</strong>=100,
                         <strong>offspring_size</strong>=None, <strong>mutation_rate</strong>=0.9,
                         <strong>crossover_rate</strong>=0.1,
                         <strong>scoring</strong>='neg_mean_squared_error', <strong>cv</strong>=5,
                         <strong>subsample</strong>=1.0, <strong>n_jobs</strong>=1,
                         <strong>max_time_mins</strong>=None, <strong>max_eval_time_mins</strong>=5,
                         <strong>random_state</strong>=None, <strong>config_dict</strong>=None,
                         <strong>template</strong>=None,
                         <strong>warm_start</strong>=False,
                         <strong>memory</strong>=None,
                         <strong>use_dask</strong>=False,
                         <strong>periodic_checkpoint_folder</strong>=None,
                         <strong>early_stop</strong>=None,
                         <strong>verbosity</strong>=0,
                         <strong>disable_update_check</strong>=False</em>)</pre>
<div align="right"><a href="https://github.com/EpistasisLab/tpot/blob/master/tpot/base.py">source</a></div>

Automated machine learning for supervised regression tasks.

The TPOTRegressor performs an intelligent search over machine learning pipelines that can contain supervised regression models,
preprocessors, feature selection techniques, and any other estimator or transformer that follows the [scikit-learn API](http://scikit-learn.org/stable/developers/contributing.html#apis-of-scikit-learn-objects).
The TPOTRegressor will also search over the hyperparameters of all objects in the pipeline.

By default, TPOTRegressor will search over a broad range of supervised regression models, transformers, and their hyperparameters.
However, the models, transformers, and parameters that the TPOTRegressor searches over can be fully customized using the `config_dict` parameter.

Read more in the [User Guide](using/#tpot-with-code).

<table>
<tr>
<td width="20%" style="vertical-align:top; background:#F5F5F5;"><strong>Parameters:</strong></td>
<td width="80%" style="background:white;">
<strong>generations</strong>: int or None, optional (default=100)
<blockquote>
Number of iterations to the run pipeline optimization process. It must be a positive number or None. If None, the parameter <em>max_time_mins</em> must be defined as the runtime limit.
<br /><br />
Generally, TPOT will work better when you give it more generations (and therefore time) to optimize the pipeline.
<br /><br />
TPOT will evaluate <em>population_size</em> + <em>generations</em> × <em>offspring_size</em> pipelines in total.
</blockquote>

<strong>population_size</strong>: int, optional (default=100)
<blockquote>
Number of individuals to retain in the genetic programming population every generation. Must be a positive number.
<br /><br />
Generally, TPOT will work better when you give it more individuals with which to optimize the pipeline.
</blockquote>

<strong>offspring_size</strong>: int, optional (default=None)
<blockquote>
Number of offspring to produce in each genetic programming generation. Must be a positive number. By default, the number of offspring is equal to the number of population size.
</blockquote>

<strong>mutation_rate</strong>: float, optional (default=0.9)
<blockquote>
Mutation rate for the genetic programming algorithm in the range [0.0, 1.0]. This parameter tells the GP algorithm how many pipelines to apply random changes to every generation.
<br /><br />
<em>mutation_rate</em> + <em>crossover_rate</em> cannot exceed 1.0.
<br /><br />
We recommend using the default parameter unless you understand how the mutation rate affects GP algorithms.
</blockquote>

<strong>crossover_rate</strong>: float, optional (default=0.1)
<blockquote>
Crossover rate for the genetic programming algorithm in the range [0.0, 1.0]. This parameter tells the genetic programming algorithm how many pipelines to "breed" every generation.
<br /><br />
<em>mutation_rate</em> + <em>crossover_rate</em> cannot exceed 1.0.
<br /><br />
We recommend using the default parameter unless you understand how the crossover rate affects GP algorithms.
</blockquote>

<strong>scoring</strong>: string or callable, optional (default='neg_mean_squared_error')
<blockquote>
Function used to evaluate the quality of a given pipeline for the regression problem. The following built-in scoring functions can be used:
<br /><br/>
'neg_median_absolute_error', 'neg_mean_absolute_error', 'neg_mean_squared_error', 'r2'
<br /><br/>
Note that we recommend using the <em>neg</em> version of mean squared error and related metrics so TPOT will minimize (instead of maximize) the metric.
<br /><br/>
If you would like to use a custom scorer, you can pass the callable object/function with signature <em>scorer(estimator, X, y)</em>.
<br /><br/>
See the section on <a href="../using/#scoring-functions">scoring functions</a> for more details.
</blockquote>

<strong>cv</strong>: int, cross-validation generator, or an iterable, optional (default=5)
<blockquote>
Cross-validation strategy used when evaluating pipelines.
<br /><br />
Possible inputs:
<ul>
<li>integer, to specify the number of folds in a KFold,</li>
<li>An object to be used as a cross-validation generator, or</li>
<li>An iterable yielding train/test splits.</li>
</ul>
</blockquote>

<strong>subsample</strong>: float, optional (default=1.0)
<blockquote>
Fraction of training samples that are used during the TPOT optimization process. Must be in the range (0.0, 1.0].
<br /><br />
Setting <em>subsample</em>=0.5 tells TPOT to use a random subsample of half of the training data. This subsample will remain the same during the entire pipeline optimization process.
</blockquote>

<strong>n_jobs</strong>: integer, optional (default=1)
<blockquote>
Number of processes to use in parallel for evaluating pipelines during the TPOT optimization process.
<br /><br />
Setting <em>n_jobs</em>=-1 will use as many cores as available on the computer. For n_jobs below -1, (n_cpus + 1 + n_jobs) are used. Thus for n_jobs = -2, all CPUs but one are used. Beware that using multiple processes on the same machine may cause memory issues for large datasets
</blockquote>

<strong>max_time_mins</strong>: integer or None, optional (default=None)
<blockquote>
How many minutes TPOT has to optimize the pipeline.
<br /><br />
If not None, this setting will allow TPOT to run until <em>max_time_mins</em> minutes elapsed and then stop. TPOT will stop earlier if <em>generations</em> is set and all generations are already evaluated.
</blockquote>

<strong>max_eval_time_mins</strong>: float, optional (default=5)
<blockquote>
How many minutes TPOT has to evaluate a single pipeline.
<br /><br />
Setting this parameter to higher values will allow TPOT to evaluate more complex pipelines, but will also allow TPOT to run longer. Use this parameter to help prevent TPOT from wasting time on evaluating time-consuming pipelines.
</blockquote>

<strong>random_state</strong>: integer or None, optional (default=None)
<blockquote>
The seed of the pseudo random number generator used in TPOT.
<br /><br />
Use this parameter to make sure that TPOT will give you the same results each time you run it against the same data set with that seed.
</blockquote>

<strong>config_dict</strong>: Python dictionary, string, or None, optional (default=None)
<blockquote>
A configuration dictionary for customizing the operators and parameters that TPOT searches in the optimization process.
<br /><br />
Possible inputs are:
<ul>
<li>Python dictionary, TPOT will use your custom configuration,</li>
<li>string 'TPOT light', TPOT will use a built-in configuration with only fast models and preprocessors, or</li>
<li>string 'TPOT MDR', TPOT will use a built-in configuration specialized for genomic studies, or</li>
<li>string 'TPOT sparse': TPOT will use a configuration dictionary with a one-hot encoder and the operators normally included in TPOT that also support sparse matrices, or</li>
<li>None, TPOT will use the default TPOTRegressor configuration.</li>
</ul>
See the <a href="../using/#built-in-tpot-configurations">built-in configurations</a> section for the list of configurations included with TPOT, and the <a href="../using/#customizing-tpots-operators-and-parameters">custom configuration</a> section for more information and examples of how to create your own TPOT configurations.
</blockquote>

<strong>template</strong>: string (default=None)
<blockquote>
Template of predefined pipeline structure. The option is for specifying a desired structure for the machine learning pipeline evaluated in TPOT.
<br /><br />

So far this option only supports linear pipeline structure. Each step in the pipeline should be a main class of operators (Selector, Transformer or Regressor) or a specific operator (e.g. `SelectPercentile`) defined in TPOT operator configuration. If one step is a main class, TPOT will randomly assign all subclass operators (subclasses of [`SelectorMixin`](https://github.com/scikit-learn/scikit-learn/blob/master/sklearn/feature_selection/base.py#L17), [`TransformerMixin`](https://scikit-learn.org/stable/modules/generated/sklearn.base.TransformerMixin.html) or [`RegressorMixin`](https://scikit-learn.org/stable/modules/generated/sklearn.base.RegressorMixin.html) in scikit-learn) to that step. Steps in the template are delimited by "-", e.g. "SelectPercentile-Transformer-Regressor". By default value of template is None, TPOT generates tree-based pipeline randomly.

See the <a href="../using/#template-option-in-tpot"> template option in tpot</a> section for more details.
</blockquote>

<strong>warm_start</strong>: boolean, optional (default=False)
<blockquote>
Flag indicating whether the TPOT instance will reuse the population from previous calls to <em>fit()</em>.
<br /><br />
Setting <em>warm_start</em>=True can be useful for running TPOT for a short time on a dataset, checking the results, then resuming the TPOT run from where it left off.
</blockquote>

<strong>memory</strong>: a joblib.Memory object or string, optional (default=None)
<blockquote>
If supplied, pipeline will cache each transformer after calling fit. This feature is used to avoid computing the fit transformers within a pipeline if the parameters and input data are identical with another fitted pipeline during optimization process. More details about memory caching in <a href="http://scikit-learn.org/stable/modules/pipeline.html#caching-transformers-avoid-repeated-computation">scikit-learn documentation</a>
<br /><br />
Possible inputs are:
<ul>
<li>String 'auto': TPOT uses memory caching with a temporary directory and cleans it up upon shutdown, or</li>
<li>Path of a caching directory, TPOT uses memory caching with the provided directory and TPOT does NOT clean the caching directory up upon shutdown, or</li>
<li>Memory object, TPOT uses the instance of joblib.Memory for memory caching and TPOT does NOT clean the caching directory up upon shutdown, or</li>
<li>None, TPOT does not use memory caching.</li>
</ul>
</blockquote>

<strong>use_dask</strong>: boolean, optional (default: False)
<blockquote>
Whether to use Dask-ML's pipeline optimiziations. This avoid re-fitting
the same estimator on the same split of data multiple times. It
will also provide more detailed diagnostics when using Dask's
distributed scheduler.
<br /><br />
See <a href="https://dask-ml.readthedocs.io/en/latest/hyper-parameter-search.html#avoid-repeated-work">avoid repeated work</a> for more details.
</blockquote>

<strong>periodic_checkpoint_folder</strong>: path string, optional (default: None)
<blockquote>
If supplied, a folder in which TPOT will periodically save pipelines in pareto front so far while optimizing.<br /><br />
Currently once per generation but not more often than once per 30 seconds.<br /><br />
Useful in multiple cases:
<ul>
<li>Sudden death before TPOT could save optimized pipeline</li>
<li>Track its progress</li>
<li>Grab pipelines while it's still optimizing</li>
</ul>
</blockquote>

<strong>early_stop</strong>: integer, optional (default: None)
<blockquote>
How many generations TPOT checks whether there is no improvement in optimization process.
<br /><br />
Ends the optimization process if there is no improvement in the given number of generations.
</blockquote>

<strong>verbosity</strong>: integer, optional (default=0)
<blockquote>
How much information TPOT communicates while it's running.
<br /><br />
Possible inputs are:
<ul>
<li>0, TPOT will print nothing,</li>
<li>1, TPOT will print minimal information,</li>
<li>2, TPOT will print more information and provide a progress bar, or</li>
<li>3, TPOT will print everything and provide a progress bar.</li>
</ul>
</blockquote>

<strong>disable_update_check</strong>: boolean, optional (default=False)
<blockquote>
Flag indicating whether the TPOT version checker should be disabled.
<br /><br />
The update checker will tell you when a new version of TPOT has been released.
</blockquote>
</td>
</tr>

<tr>
<td width="20%" style="vertical-align:top; background:#F5F5F5;"><strong>Attributes:</strong></td>
<td width="80%" style="background:white;">
<strong>fitted_pipeline_</strong>: scikit-learn Pipeline object
<blockquote>
The best pipeline that TPOT discovered during the pipeline optimization process, fitted on the entire training dataset.
</blockquote>

<strong>pareto_front_fitted_pipelines_</strong>: Python dictionary
<blockquote>
Dictionary containing the all pipelines on the TPOT Pareto front, where the key is the string representation of the pipeline and the value is the corresponding pipeline fitted on the entire training dataset.
<br /><br />
The TPOT Pareto front provides a trade-off between pipeline complexity (i.e., the number of steps in the pipeline) and the predictive performance of the pipeline.
<br /><br />
Note: <em>_pareto_front_fitted_pipelines</em> is only available when <em>verbosity</em>=3.
</blockquote>

<strong>evaluated_individuals_</strong>: Python dictionary
<blockquote>
Dictionary containing all pipelines that were evaluated during the pipeline optimization process, where the key is the string representation of the pipeline and the value is a tuple containing (# of steps in pipeline, accuracy metric for the pipeline).
<br /><br />
This attribute is primarily for internal use, but may be useful for looking at the other pipelines that TPOT evaluated.
</blockquote>
</td>
<tr>
</table>

<strong>Example</strong>

```Python
from tpot import TPOTRegressor
from sklearn.datasets import load_boston
from sklearn.model_selection import train_test_split

digits = load_boston()
X_train, X_test, y_train, y_test = train_test_split(digits.data, digits.target,
                                                    train_size=0.75, test_size=0.25)

tpot = TPOTRegressor(generations=5, population_size=50, verbosity=2)
tpot.fit(X_train, y_train)
print(tpot.score(X_test, y_test))
tpot.export('tpot_boston_pipeline.py')
```

<strong>Functions</strong>

<table width="100%">
<tr>
<td width="25%"><a href="#tpotregressor-fit">fit</a>(features, target[, sample_weight, groups])</td>
<td>Run the TPOT optimization process on the given training data.</td>
</tr>

<tr>
<td><a href="#tpotregressor-predict">predict</a>(features)</td>
<td>Use the optimized pipeline to predict the target values for a feature set.</td>
</tr>

<tr>
<td><a href="#tpotregressor-score">score</a>(testing_features, testing_target)</td>
<td>Returns the optimized pipeline's score on the given testing data using the user-specified scoring function.</td>
</tr>

<tr>
<td><a href="#tpotregressor-export">export</a>(output_file_name)</td>
<td>Export the optimized pipeline as Python code.</td>
</tr>
</table>


<a name="tpotregressor-fit"></a>
```Python
fit(features, target, sample_weight=None, groups=None)
```

<div style="padding-left:5%" width="100%">
Run the TPOT optimization process on the given training data.
<br /><br />
Uses genetic programming to optimize a machine learning pipeline that maximizes the score on the provided features and target. This pipeline optimization procedure uses internal k-fold cross-validaton to avoid overfitting on the provided data. At the end of the pipeline optimization procedure, the best pipeline is then trained on the entire set of provided samples.
<br /><br />
<table width="100%">
<tr>
<td width="20%" style="vertical-align:top; background:#F5F5F5;"><strong>Parameters:</strong></td>
<td width="80%" style="background:white;">
<strong>features</strong>: array-like {n_samples, n_features}
<blockquote>
Feature matrix
<br /><br />
TPOT and all scikit-learn algorithms assume that the features will be numerical and there will be no missing values.
As such, when a feature matrix is provided to TPOT, all missing values will automatically be replaced (i.e., imputed)
using <a href="http://scikit-learn.org/stable/modules/generated/sklearn.impute.SimpleImputer.html">median value imputation</a>.
<br /><br />
If you wish to use a different imputation strategy than median imputation, please make sure to apply imputation to your feature set prior to passing it to TPOT.
</blockquote>

<strong>target</strong>: array-like {n_samples}
<blockquote>
List of target labels for prediction
</blockquote>

<strong>sample_weight</strong>: array-like {n_samples}, optional
<blockquote>
Per-sample weights. Higher weights indicate more importance. If specified, sample_weight will be passed to any pipeline element whose fit() function accepts a sample_weight argument. By default, using sample_weight does not affect tpot's scoring functions, which determine preferences between pipelines.
</blockquote>

<strong>groups</strong>: array-like, with shape {n_samples, }, optional
<blockquote>
Group labels for the samples used when performing cross-validation.
<br /><br />
This parameter should only be used in conjunction with sklearn's Group cross-validation functions, such as <a href="http://scikit-learn.org/stable/modules/generated/sklearn.model_selection.GroupKFold.html">sklearn.model_selection.GroupKFold</a>.
</blockquote>
</td>
</tr>
<tr>
<td width="20%" style="vertical-align:top; background:#F5F5F5;"><strong>Returns:</strong></td>
<td width="80%" style="background:white;">
<strong>self</strong>: object
<blockquote>
Returns a copy of the fitted TPOT object
</blockquote>
</td>
</tr>
</table>
</div>


<a name="tpotregressor-predict"></a>
```Python
predict(features)
```

<div style="padding-left:5%" width="100%">
Use the optimized pipeline to predict the target values for a feature set.
<br /><br />
<table width="100%">
<tr>
<td width="20%" style="vertical-align:top; background:#F5F5F5;"><strong>Parameters:</strong></td>
<td width="80%" style="background:white;">
<strong>features</strong>: array-like {n_samples, n_features}
<blockquote>
Feature matrix
</blockquote>
</td>
</tr>
<tr>
<td width="20%" style="vertical-align:top; background:#F5F5F5;"><strong>Returns:</strong></td>
<td width="80%" style="background:white;">
<strong>predictions</strong>: array-like {n_samples}
<blockquote>
Predicted target values for the samples in the feature matrix
</blockquote>
</td>
</tr>
</table>
</div>


<a name="tpotregressor-score"></a>
```Python
score(testing_features, testing_target)
```

<div style="padding-left:5%" width="100%">
Returns the optimized pipeline's score on the given testing data using the user-specified scoring function.
<br /><br />
The default scoring function for TPOTRegressor is 'mean_squared_error'.
<br /><br />
<table width="100%">
<tr>
<td width="20%" style="vertical-align:top; background:#F5F5F5;"><strong>Parameters:</strong></td>
<td width="80%" style="background:white;">
<strong>testing_features</strong>: array-like {n_samples, n_features}
<blockquote>
Feature matrix of the testing set
</blockquote>

<strong>testing_target</strong>: array-like {n_samples}
<blockquote>
List of target labels for prediction in the testing set
</blockquote>
</td>
</tr>
<tr>
<td width="20%" style="vertical-align:top; background:#F5F5F5;"><strong>Returns:</strong></td>
<td width="80%" style="background:white;">
<strong>accuracy_score</strong>: float
<blockquote>
The estimated test set accuracy according to the user-specified scoring function.
</blockquote>
</td>
</tr>
</table>
</div>


<a name="tpotregressor-export"></a>
```Python
export(output_file_name)
```

<div style="padding-left:5%" width="100%">
Export the optimized pipeline as Python code.
<br /><br />
See the <a href="../using/#tpot-with-code">usage documentation</a> for example usage of the export function.
<br /><br />
<table width="100%">
<tr>
<td width="20%" style="vertical-align:top; background:#F5F5F5;"><strong>Parameters:</strong></td>
<td width="80%" style="background:white;">
<strong>output_file_name</strong>: string
<blockquote>
String containing the path and file name of the desired output file
</blockquote>
<strong>data_file_path</strong>: string
<blockquote>
By default, the path of input dataset is 'PATH/TO/DATA/FILE' by default. If data_file_path is another string, the path will be replaced.
</blockquote>
</tr>
<tr>
<td width="20%" style="vertical-align:top; background:#F5F5F5;"><strong>Returns:</strong></td>
<td width="80%" style="background:white;">
<strong>exported_code_string</strong>: string
<blockquote>
The whole pipeline text as a string should be returned if output_file_name is not specified.
</blockquote>
</td>
</tr>
</table>
</div>
# Release Notes

## Version 0.11.7

- Fix compatibility issue with scikit-learn 0.24 and xgboost 1.3.0
- Fix a bug causing that TPOT does not work when classifying more than 50 classes
- Add initial support `Resampler` from `imblearn`
- Fix minor bugs


## Version 0.11.6

- Fix a bug causing point mutation function does not work properly with using `template` option
- Add a new built configuration called "TPOT cuML" which TPOT will search over a restricted configuration using the GPU-accelerated estimators in [RAPIDS cuML](https://github.com/rapidsai/cuml) and [DMLC XGBoost](https://github.com/dmlc/xgboost). **This configuration requires an NVIDIA Pascal architecture or better GPU with [compute capability 6.0+](https://developer.nvidia.com/cuda-gpus), and that the library cuML is installed.**
- Add string path support for log/log_file parameter
- Fix a bug in version 0.11.5 causing no update in stdout after each generation
- Fix minor bugs


## Version 0.11.5

- Make `Pytorch` as an optional dependency
- Refine installation documentation

## Version 0.11.4

- Add a new built configuration "TPOT NN" which includes all operators in "Default TPOT" plus additional neural network estimators written in PyTorch (currently `tpot.builtins.PytorchLRClassifier` and `tpot.builtins.PytorchMLPClassifier` for classification tasks only)
- Refine `log_file` parameter's behavior

## Version 0.11.3

- Fix a bug in TPOTRegressor in v0.11.2
- Add `-log` option in command line interface to save process log to a file.

## Version 0.11.2

- Fix `early_stop` parameter does not work properly
- TPOT built-in `OneHotEncoder` can refit to different datasets
- Fix the issue that the attribute `evaluated_individuals_` cannot record correct generation info.
- Add a new parameter `log_file` to output logs to a file instead of `sys.stdout`
- Fix some code quality issues and mistakes in documentations
- Fix minor bugs

## Version 0.11.1

- Fix compatibility issue with scikit-learn v0.22
- `warm_start` now saves both Primitive Sets and evaluated_pipelines_ from previous runs;
- Fix the error that TPOT assign wrong fitness scores to non-evaluated pipelines (interrupted by `max_min_mins` or `KeyboardInterrupt`) ;
- Fix the bug that mutation operator cannot generate new pipeline when template is not default value and `warm_start` is True;
- Fix the bug that `max_time_mins` cannot stop optimization process when search space is limited.  
- Fix a bug in exported codes when the exported pipeline is only 1 estimator
- Fix spelling mistakes in documentations
- Fix some code quality issues

## Version 0.11.0

- **Support for Python 3.4 and below has been officially dropped.** Also support for scikit-learn 0.20 or below has been dropped.
- The support of a metric function with the signature `score_func(y_true, y_pred)` for `scoring parameter` has been dropped.
- Refine `StackingEstimator` for not stacking NaN/Infinity predication probabilities.
- Fix a bug that population doesn't persist by `warm_start=True` when `max_time_mins` is not default value.
- Now the `random_state` parameter in TPOT is used for pipeline evaluation instead of using a fixed random seed of 42 before. The `set_param_recursive` function has been moved to `export_utils.py` and it can be used in exported codes for setting `random_state` recursively in scikit-learn Pipeline. It is used to set `random_state` in `fitted_pipeline_` attribute and exported pipelines.
- TPOT can independently use `generations` and `max_time_mins` to limit the optimization process through using one of the parameters or both.
- `.export()` function will return string of exported pipeline if output filename is not specified.
- Add [`SGDClassifier`](https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.SGDClassifier.html) and [`SGDRegressor`](https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.SGDRegressor.html) into TPOT default configs.
- Documentation has been updated
- Fix minor bugs.

## Version 0.10.2

- **TPOT v0.10.2 is the last version to support Python 2.7 and Python 3.4.**
- Minor updates for fixing compatibility issues with the latest version of scikit-learn (version > 0.21) and xgboost (v0.90)
- Default value of `template` parameter is changed to `None` instead.
- Fix errors in documentation

## Version 0.10.1

- Add `data_file_path` option into `expert` function for replacing `'PATH/TO/DATA/FILE'` to customized dataset path in exported scripts. (Related issue #838)
- Change python version in CI tests to 3.7
- Add CI tests for macOS.

## Version 0.10.0

- Add a new `template` option to specify a desired structure for machine learning pipeline in TPOT. Check [TPOT API](https://epistasislab.github.io/tpot/api/) (it will be updated once it is merge to master branch).
- Add `FeatureSetSelector` operator into TPOT for feature selection based on *priori* export knowledge. Please check our [preprint paper](https://www.biorxiv.org/content/10.1101/502484v1.article-info) for more details (*Note: it was named `DatasetSelector` in 1st version paper but we will rename to FeatureSetSelector in next version of the paper*)
- Refine `n_jobs` parameter to accept value below -1. For n_jobs below -1, (n_cpus + 1 + n_jobs) are used. Thus for n_jobs = -2, all CPUs but one are used.
- Now `memory`  parameter can create memory cache directory if it does not exist.
- Fix minor bugs.

## Version 0.9.6

- Fix a bug causing that `max_time_mins` parameter doesn't work when `use_dask=True` in TPOT 0.9.5
- Now TPOT saves best pareto values best pareto pipeline s in checkpoint folder
- TPOT raises `ImportError` if operators in the TPOT configuration are not available when `verbosity>2`
- Thank @PGijsbers for the suggestions. Now TPOT can save scores of individuals already evaluated in any generation even the evaluation process of that generation is interrupted/stopped. But it is noted that, in this case, TPOT will raise this **warning message**: `WARNING: TPOT may not provide a good pipeline if TPOT is stopped/interrupted in a early generation.`, because the pipelines in early generation, e.g. 1st generation, are evolved/modified very limited times via evolutionary algorithm.
- Fix bugs in configuration of `TPOTRegressor`
- Error fixes in documentation

## Version 0.9.5

- **TPOT now supports integration with Dask for parallelization + smart caching**. Big thanks to the Dask dev team for making this happen!

- TPOT now supports for imputation/sparse matrices into `predict` and `predict_proba` functions.

- `TPOTClassifier` and `TPOTRegressor` now follows scikit-learn estimator API.

- We refined scoring parameter in TPOT API for accepting [`Scorer` object](http://jaquesgrobler.github.io/online-sklearn-build/modules/generated/sklearn.metrics.Scorer.html).

- We refined parameters in VarianceThreshold and FeatureAgglomeration.

- TPOT now supports using memory caching within a Pipeline via an optional `memory` parameter.

- We improved documentation of TPOT.

## Version 0.9

* **TPOT now supports sparse matrices** with a new built-in TPOT configuration, "TPOT sparse". We are using a custom OneHotEncoder implementation that supports missing values and continuous features.

* We have added an "early stopping" option for stopping the optimization process if no improvement is made within a set number of generations. Look up the `early_stop` parameter to access this functionality.

* TPOT now reduces the number of duplicated pipelines between generations, which saves you time during the optimization process.

* TPOT now supports custom scoring functions via the command-line mode.

* We have added a new optional argument, `periodic_checkpoint_folder`, that allows TPOT to periodically save the best pipeline so far to a local folder during optimization process.

* TPOT no longer uses `sklearn.externals.joblib` when `n_jobs=1` to avoid the potential freezing issue [that scikit-learn suffers from](http://scikit-learn.org/stable/faq.html#why-do-i-sometime-get-a-crash-freeze-with-n-jobs-1-under-osx-or-linux).

* We have added `pandas` as a dependency to read input datasets instead of `numpy.recfromcsv`. NumPy's `recfromcsv` function is unable to parse datasets with complex data types.

* Fixed a bug that `DEFAULT` in the parameter(s) of nested estimator raises `KeyError` when exporting pipelines.

* Fixed a bug related to setting `random_state` in nested estimators. The issue would happen with pipeline with `SelectFromModel` (`ExtraTreesClassifier` as nested estimator) or `StackingEstimator` if nested estimator has `random_state` parameter.

* Fixed a bug in the missing value imputation function in TPOT to impute along columns instead rows.

* Refined input checking for sparse matrices in TPOT.

* Refined the TPOT pipeline mutation operator.


## Version 0.8

* **TPOT now detects whether there are missing values in your dataset** and replaces them with the median value of the column.

* TPOT now allows you to set a `group` parameter in the `fit` function so you can use the [GroupKFold](http://scikit-learn.org/stable/modules/generated/sklearn.model_selection.GroupKFold.html) cross-validation strategy.

* TPOT now allows you to set a subsample ratio of the training instance with the `subsample` parameter. For example, setting `subsample`=0.5 tells TPOT to create a fixed subsample of half of the training data for the pipeline optimization process. This parameter can be useful for speeding up the pipeline optimization process, but may give less accurate performance estimates from cross-validation.

* **TPOT now has more [built-in configurations](/using/#built-in-tpot-configurations)**, including TPOT MDR and TPOT light, for both classification and regression problems.

* `TPOTClassifier` and `TPOTRegressor` now expose three useful internal attributes, `fitted_pipeline_`, `pareto_front_fitted_pipelines_`, and `evaluated_individuals_`. These attributes are described in the [API documentation](/api/).

* Oh, **TPOT now has [thorough API documentation](/api/)**. Check it out!

* Fixed a reproducibility issue where setting `random_seed` didn't necessarily result in the same results every time. This bug was present since TPOT v0.7.

* Refined input checking in TPOT.

* Removed Python 2 uncompliant code.


## Version 0.7

* **TPOT now has multiprocessing support.** TPOT allows you to use multiple processes in parallel to accelerate the pipeline optimization process in TPOT with the `n_jobs` parameter.

* TPOT now allows you to **customize the operators and parameters considered during the optimization process**, which can be accomplished with the new `config_dict` parameter. The format of this customized dictionary can be found in the [online documentation](/using/#customizing-tpots-operators-and-parameters), along with a list of [built-in configurations](/using/#built-in-tpot-configurations).

* TPOT now allows you to **specify a time limit for evaluating a single pipeline**  (default limit is 5 minutes) in optimization process with the `max_eval_time_mins` parameter, so TPOT won't spend hours evaluating overly-complex pipelines.

* We tweaked TPOT's underlying evolutionary optimization algorithm to work even better, including using the [mu+lambda algorithm](http://deap.readthedocs.io/en/master/api/algo.html#deap.algorithms.eaMuPlusLambda). This algorithm gives you more control of how many pipelines are generated every iteration with the `offspring_size` parameter.

* Refined the default operators and parameters in TPOT, so TPOT 0.7 should work even better than 0.6.

* TPOT now supports sample weights in the fitness function if some if your samples are more important to classify correctly than others. The sample weights option works the same as in scikit-learn, e.g., `tpot.fit(x_train, y_train, sample_weights=sample_weights)`.

* The default scoring metric in TPOT has been changed from balanced accuracy to accuracy, the same default metric for classification algorithms in scikit-learn. Balanced accuracy can still be used by setting `scoring='balanced_accuracy'` when creating a TPOT instance.


## Version 0.6

* **TPOT now supports regression problems!** We have created two separate `TPOTClassifier` and `TPOTRegressor` classes to support classification and regression problems, respectively. The [command-line interface](/using/#tpot-on-the-command-line) also supports this feature through the `-mode` parameter.

* TPOT now allows you to **specify a time limit** for the optimization process with the `max_time_mins` parameter, so you don't need to guess how long TPOT will take any more to recommend a pipeline to you.

* Added a new operator that performs feature selection using [ExtraTrees](http://scikit-learn.org/stable/modules/ensemble.html#extremely-randomized-trees) feature importance scores.

* **[XGBoost](https://github.com/dmlc/xgboost) has been added as an optional dependency to TPOT.** If you have XGBoost installed, TPOT will automatically detect your installation and use the `XGBoostClassifier` and `XGBoostRegressor` in its pipelines.

* TPOT now offers a verbosity level of 3 ("science mode"), which outputs the entire Pareto front instead of only the current best score. This feature may be useful for users looking to make a trade-off between pipeline complexity and score.

## Version 0.5

* Major refactor: Each operator is defined in a separate class file. Hooray for easier-to-maintain code!
* TPOT now **exports directly to scikit-learn Pipelines** instead of hacky code.
* Internal representation of individuals now uses scikit-learn pipelines.
* Parameters for each operator have been optimized so TPOT spends less time exploring useless parameters.
* We have removed pandas as a dependency and instead use numpy matrices to store the data.
* TPOT now uses **k-fold cross-validation** when evaluating pipelines, with a default k = 3. This k parameter can be tuned when creating a new TPOT instance.
* Improved **scoring function support**: Even though TPOT uses balanced accuracy by default, you can now have TPOT use [any of the scoring functions](http://scikit-learn.org/stable/modules/model_evaluation.html#common-cases-predefined-values) that `cross_val_score` supports.
* Added the scikit-learn [Normalizer](http://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.Normalizer.html) preprocessor.
* [Minor text fixes.](http://knowyourmeme.com/memes/pokemon-go-updates-controversy)

## Version 0.4

In TPOT 0.4, we've made some major changes to the internals of TPOT and added some convenience functions. We've summarized the changes below.

<ul>
<li>Added new sklearn models and preprocessors

<ul>
<li>AdaBoostClassifier</li>
<li>BernoulliNB</li>
<li>ExtraTreesClassifier</li>
<li>GaussianNB</li>
<li>MultinomialNB</li>
<li>LinearSVC</li>
<li>PassiveAggressiveClassifier</li>
<li>GradientBoostingClassifier</li>
<li>RBFSampler</li>
<li>FastICA</li>
<li>FeatureAgglomeration</li>
<li>Nystroem</li>
</ul></li>
<li>Added operator that inserts virtual features for the count of features with values of zero</li>
<li>Reworked parameterization of TPOT operators
<ul>
<li>Reduced parameter search space with information from a scikit-learn benchmark</li>
<li>TPOT no longer generates arbitrary parameter values, but uses a fixed parameter set instead</li>
</ul></li>
<li>Removed XGBoost as a dependency
<ul>
<li>Too many users were having install issues with XGBoost</li>
<li>Replaced with scikit-learn's GradientBoostingClassifier</li>
</ul></li>
<li>Improved descriptiveness of TPOT command line parameter documentation</li>
<li>Removed min/max/avg details during fit() when verbosity &gt; 1

<ul>
<li>Replaced with tqdm progress bar</li>
<li>Added tqdm as a dependency</li>
</ul></li>
<li>Added <code>fit_predict()</code> convenience function</li>
<li>Added <code>get_params()</code> function so TPOT can operate in scikit-learn's <code>cross_val_score</code> & related functions</li>
</ul>

## Version 0.3

* We revised the internal optimization process of TPOT to make it more efficient, in particular in regards to the model parameters that TPOT optimizes over.

## Version 0.2

* TPOT now has the ability to export the optimized pipelines to sklearn code.

* Logistic regression, SVM, and k-nearest neighbors classifiers were added as pipeline operators. Previously, TPOT only included decision tree and random forest classifiers.

* TPOT can now use arbitrary scoring functions for the optimization process.

* TPOT now performs multi-objective Pareto optimization to balance model complexity (i.e., # of pipeline operators) and the score of the pipeline.

## Version 0.1

* First public release of TPOT.

* Optimizes pipelines with decision trees and random forest classifiers as the model, and uses a handful of feature preprocessors.
# Contribution Guide

We welcome you to [check the existing issues](https://github.com/EpistasisLab/tpot/issues/) for bugs or enhancements to work on. If you have an idea for an extension to TPOT, please [file a new issue](https://github.com/EpistasisLab/tpot/issues/new) so we can discuss it.

## Project layout

The latest stable release of TPOT is on the [master branch](https://github.com/EpistasisLab/tpot/tree/master), whereas the latest version of TPOT in development is on the [development branch](https://github.com/EpistasisLab/tpot/tree/development). Make sure you are looking at and working on the correct branch if you're looking to contribute code.

In terms of directory structure:

* All of TPOT's code sources are in the `tpot` directory
* The documentation sources are in the `docs_sources` directory
* Images in the documentation are in the `images` directory
* Tutorials for TPOT are in the `tutorials` directory
* Unit tests for TPOT are in the `tests.py` file

Make sure to familiarize yourself with the project layout before making any major contributions, and especially make sure to send all code changes to the `development` branch.

## How to contribute

The preferred way to contribute to TPOT is to fork the
[main repository](https://github.com/EpistasisLab/tpot/) on
GitHub:

1. Fork the [project repository](https://github.com/EpistasisLab/tpot):
   click on the 'Fork' button near the top of the page. This creates
   a copy of the code under your account on the GitHub server.

2. Clone this copy to your local disk:

          $ git clone git@github.com:YourUsername/tpot.git
          $ cd tpot

3. Create a branch to hold your changes:

          $ git checkout -b my-contribution

4. Make sure your local environment is setup correctly for development. Installation instructions are almost identical to [the user instructions](installing.md) except that TPOT should *not* be installed. If you have TPOT installed on your computer then make sure you are using a virtual environment that does not have TPOT installed. Furthermore, you should make sure you have installed the `nose` package into your development environment so that you can test changes locally.

          $ conda install nose

5. Start making changes on your newly created branch, remembering to never work on the ``master`` branch! Work on this copy on your computer using Git to do the version control.

6. Once some changes are saved locally, you can use your tweaked version of TPOT by navigating to the project's base directory and running TPOT directly from the command line:

          $ python -m tpot.driver

    or by running script that imports and uses the TPOT module with code similar to `from tpot import TPOTClassifier`

7. To check your changes haven't broken any existing tests and to check new tests you've added pass run the following (note, you must have the `nose` package installed within your dev environment for this to work):

          $ nosetests -s -v

8. When you're done editing and local testing, run:

          $ git add modified_files
          $ git commit

   to record your changes in Git, then push them to GitHub with:

          $ git push -u origin my-contribution

Finally, go to the web page of your fork of the TPOT repo, and click 'Pull Request' (PR) to send your changes to the maintainers for review. Make sure that you send your PR to the `development` branch, as the `master` branch is reserved for the latest stable release. This will start the CI server to check all the project's unit tests run and send an email to the maintainers.

(If any of the above seems like magic to you, then look up the
[Git documentation](http://git-scm.com/documentation) on the web.)

## Before submitting your pull request

Before you submit a pull request for your contribution, please work through this checklist to make sure that you have done everything necessary so we can efficiently review and accept your changes.

If your contribution changes TPOT in any way:

* Update the [documentation](https://github.com/EpistasisLab/tpot/tree/master/docs_sources) so all of your changes are reflected there.

* Update the [README](https://github.com/EpistasisLab/tpot/blob/master/README.md) if anything there has changed.

If your contribution involves any code changes:

* Update the [project unit tests](https://github.com/EpistasisLab/tpot/tree/master/tests) to test your code changes.

* Make sure that your code is properly commented with [docstrings](https://www.python.org/dev/peps/pep-0257/) and comments explaining your rationale behind non-obvious coding practices.

* If your code affected any of the pipeline operators, make sure that the corresponding [export functionality](https://github.com/EpistasisLab/tpot/blob/master/tpot/export_utils.py) reflects those changes.

If your contribution requires a new library dependency:

* Double-check that the new dependency is easy to install via `pip` or Anaconda and supports both Python 2 and 3. If the dependency requires a complicated installation, then we most likely won't merge your changes because we want to keep TPOT easy to install.

* Add the required version of the library to [.travis.yml](https://github.com/EpistasisLab/tpot/blob/master/.travis.yml#L7)

* Add a line to pip install the library to [.travis_install.sh](https://github.com/EpistasisLab/tpot/blob/master/ci/.travis_install.sh#L46)

* Add a line to print the version of the library to [.travis_install.sh](https://github.com/EpistasisLab/tpot/blob/master/ci/.travis_install.sh#L63)

* Similarly add a line to print the version of the library to [.travis_test.sh](https://github.com/EpistasisLab/tpot/blob/master/ci/.travis_test.sh#L13)

## After submitting your pull request

After submitting your pull request, [Travis-CI](https://travis-ci.com/) will automatically run unit tests on your changes and make sure that your updated code builds and runs on Python 2 and 3. We also use services that automatically check code quality and test coverage.

Check back shortly after submitting your pull request to make sure that your code passes these checks. If any of the checks come back with a red X, then do your best to address the errors.
Other Automated Machine Learning (AutoML) tools and related projects:

<table>
<tr>
<th width="20%">Name</th>
<th width="15%">Language</th>
<th width="15%">License</th>
<th>Description</th>
</tr>
<tr>
<td><a href="http://www.cs.ubc.ca/labs/beta/Projects/autoweka/">Auto-WEKA</a></td>
<td>Java</td>
<td>GPL-v3</td>
<td>Automated model selection and hyper-parameter tuning for Weka models.</td>
</tr>
<tr>
<td><a href="https://github.com/automl/auto-sklearn">auto-sklearn</a></td>
<td>Python</td>
<td>BSD-3-Clause</td>
<td>An automated machine learning toolkit and a drop-in replacement for a scikit-learn estimator.</td>
</tr>
<tr>
<td><a href="https://github.com/ClimbsRocks/auto_ml">auto_ml</a></td>
<td>Python</td>
<td>MIT</td>
<td>Automated machine learning for analytics & production. Supports manual feature type declarations.</td>
</tr>
<tr>
<td><a href="http://docs.h2o.ai/h2o/latest-stable/h2o-docs/automl.html">H2O AutoML</a></td>
<td>Java with Python, Scala & R APIs and web GUI</td>
<td>Apache 2.0</td>
<td>Automated: data prep, hyperparameter tuning, random grid search and stacked ensembles in a distributed ML platform.</td>
</tr>
<tr>
<td><a href="https://github.com/joeddav/devol">devol</a></td>
<td>Python</td>
<td>MIT</td>
<td>Automated deep neural network design via genetic programming.</td>
</tr>
<tr>
<td><a href="https://github.com/AxeldeRomblay/MLBox">MLBox</a></td>
<td>Python</td>
<td>BSD-3-Clause</td>
<td>Accurate hyper-parameter optimization in high-dimensional space with support for distributed computing.</td>
</tr>
<tr>
<td><a href="https://github.com/RecipeML/Recipe">Recipe</a></td>
<td>C</td>
<td>GPL-v3</td>
<td>Machine-learning pipeline optimization through genetic programming. Uses grammars to define pipeline structure.</td>
</tr>
<tr>
<td><a href="https://github.com/reiinakano/xcessiv">Xcessiv</a></td>
<td>Python</td>
<td>Apache 2.0</td>
<td>A web-based application for quick, scalable, and automated hyper-parameter tuning and stacked ensembling in Python.</td>
</tr>
<tr>
<td><a href="https://github.com/PGijsbers/gama">GAMA</a></td>
<td>Python</td>
<td>Apache 2.0</td>
<td>Machine-learning pipeline optimization through asynchronous evaluation based genetic programming. </td>
</tr>
</table>
# Installation

TPOT is built on top of several existing Python libraries, including:

* [NumPy](http://www.numpy.org/)

* [SciPy](https://www.scipy.org/)

* [scikit-learn](http://www.scikit-learn.org/)

* [DEAP](https://github.com/DEAP/deap)

* [update_checker](https://github.com/bboe/update_checker)

* [tqdm](https://github.com/tqdm/tqdm)

* [stopit](https://github.com/glenfant/stopit)

* [pandas](http://pandas.pydata.org)

* [joblib](https://joblib.readthedocs.io/en/latest/)

* [xgboost](https://xgboost.readthedocs.io/en/latest/)

Most of the necessary Python packages can be installed via the [Anaconda Python distribution](https://www.anaconda.com/products/individual), which we strongly recommend that you use. **Support for Python 3.4 and below has been officially dropped since version 0.11.0.**


You can install TPOT using `pip` or `conda-forge`.

## pip

NumPy, SciPy, scikit-learn, pandas, joblib, and PyTorch can be installed in Anaconda via the command:

```Shell
conda install numpy scipy scikit-learn pandas joblib pytorch
```

DEAP, update_checker, tqdm, stopit and xgboost can be installed with `pip` via the command:

```Shell
pip install deap update_checker tqdm stopit xgboost
```

**Windows users: pip installation may not work on some Windows environments, and it may cause unexpected errors.** If you have issues installing XGBoost, check the [XGBoost installation documentation](http://xgboost.readthedocs.io/en/latest/build.html).

If you plan to use [Dask](http://dask.pydata.org/en/latest/) for parallel training, make sure to install [dask[delay] and dask[dataframe]](https://docs.dask.org/en/latest/install.html) and [dask_ml](https://dask-ml.readthedocs.io/en/latest/install.html). **It is noted that dask-ml>=1.7 requires distributed>=2.4.0 and scikit-learn>=0.23.0.**

```Shell
pip install dask[delayed] dask[dataframe] dask-ml fsspec>=0.3.3 distributed>=2.10.0
```

If you plan to use the [TPOT-MDR configuration](https://arxiv.org/abs/1702.01780), make sure to install [scikit-mdr](https://github.com/EpistasisLab/scikit-mdr) and [scikit-rebate](https://github.com/EpistasisLab/scikit-rebate):

```Shell
pip install scikit-mdr skrebate
```

To enable support for [PyTorch](https://pytorch.org/)-based neural networks (TPOT-NN), you will need to install PyTorch. TPOT-NN will work with either CPU or GPU PyTorch, but we strongly recommend using a GPU version, if possible, as CPU PyTorch models tend to train very slowly.

We recommend following [PyTorch's installation instructions](https://pytorch.org/get-started/locally/) customized for your operating system and Python distribution.

Finally to install TPOT itself, run the following command:

```Shell
pip install tpot
```

## conda-forge

To install tpot and its core dependencies you can use:

```Shell
conda install -c conda-forge tpot
```

To install additional dependencies you can use:

```Shell
conda install -c conda-forge tpot xgboost dask dask-ml scikit-mdr skrebate
```

As mentioned above, we recommend following [PyTorch's installation instructions](https://pytorch.org/get-started/locally/) for installing it to enable support for [PyTorch](https://pytorch.org/)-based neural networks (TPOT-NN).

## Installation for using TPOT-cuML configuration

With "TPOT cuML" configuration (see <a href="../using/#built-in-tpot-configurations">built-in configurations</a>), TPOT will search over a restricted configuration using the GPU-accelerated estimators in [RAPIDS cuML](https://github.com/rapidsai/cuml) and [DMLC XGBoost](https://github.com/dmlc/xgboost). **This configuration requires an NVIDIA Pascal architecture or better GPU with [compute capability 6.0+](https://developer.nvidia.com/cuda-gpus), and that the library cuML is installed.** With this configuration, all model training and predicting will be GPU-accelerated. This configuration is particularly useful for medium-sized and larger datasets on which CPU-based estimators are a common bottleneck, and works for both the `TPOTClassifier` and `TPOTRegressor`.

Please download this conda environment <a href="https://github.com/EpistasisLab/tpot/blob/master/tpot-cuml.yml">yml file</a></td> to install TPOT for using TPOT-cuML configuration.

```
conda env create -f tpot-cuml.yml -n tpot-cuml
conda activate tpot-cuml
```


## Installation problems

Please [file a new issue](https://github.com/EpistasisLab/tpot/issues/new) if you run into installation problems.
# Using TPOT

## What to expect from AutoML software

Automated machine learning (AutoML) takes a higher-level approach to machine learning than most practitioners are used to,
so we've gathered a handful of guidelines on what to expect when running AutoML software such as TPOT.

<h5>AutoML algorithms aren't intended to run for only a few minutes</h5>

Of course, you *can* run TPOT for only a few minutes and it will find a reasonably good pipeline for your dataset.
However, if you don't run TPOT for long enough, it may not find the best possible pipeline for your dataset. It may even not
find any suitable pipeline at all, in which case a `RuntimeError('A pipeline has not yet been optimized. Please call fit() first.')`
will be raised.
Often it is worthwhile to run multiple instances of TPOT in parallel for a long time (hours to days) to allow TPOT to thoroughly search
the pipeline space for your dataset.

<h5>AutoML algorithms can take a long time to finish their search</h5>

AutoML algorithms aren't as simple as fitting one model on the dataset; they are considering multiple machine learning algorithms
(random forests, linear models, SVMs, etc.) in a pipeline with multiple preprocessing steps (missing value imputation, scaling,
PCA, feature selection, etc.), the hyperparameters for all of the models and preprocessing steps, as well as multiple ways
to ensemble or stack the algorithms within the pipeline.

As such, TPOT will take a while to run on larger datasets, but it's important to realize why. With the default TPOT settings
(100 generations with 100 population size), TPOT will evaluate 10,000 pipeline configurations before finishing.
To put this number into context, think about a grid search of 10,000 hyperparameter combinations for a machine learning algorithm
and how long that grid search will take. That is 10,000 model configurations to evaluate with 10-fold cross-validation,
which means that roughly 100,000 models are fit and evaluated on the training data in one grid search.
That's a time-consuming procedure, even for simpler models like decision trees.

Typical TPOT runs will take hours to days to finish (unless it's a small dataset), but you can always interrupt
the run partway through and see the best results so far. TPOT also [provides](/tpot/api/) a `warm_start` parameter that
lets you restart a TPOT run from where it left off.

<h5>AutoML algorithms can recommend different solutions for the same dataset</h5>

If you're working with a reasonably complex dataset or run TPOT for a short amount of time, different TPOT runs
may result in different pipeline recommendations. TPOT's optimization algorithm is stochastic in nature, which means
that it uses randomness (in part) to search the possible pipeline space. When two TPOT runs recommend different
pipelines, this means that the TPOT runs didn't converge due to lack of time *or* that multiple pipelines
perform more-or-less the same on your dataset.

This is actually an advantage over fixed grid search techniques: TPOT is meant to be an assistant that gives
you ideas on how to solve a particular machine learning problem by exploring pipeline configurations that you
might have never considered, then leaves the fine-tuning to more constrained parameter tuning techniques such
as grid search.


## TPOT with code

We've taken care to design the TPOT interface to be as similar as possible to scikit-learn.

TPOT can be imported just like any regular Python module. To import TPOT, type:

```Python
from tpot import TPOTClassifier
```

then create an instance of TPOT as follows:

```Python
pipeline_optimizer = TPOTClassifier()
```

It's also possible to use TPOT for regression problems with the `TPOTRegressor` class. Other than the class name,
a `TPOTRegressor` is used the same way as a `TPOTClassifier`. You can read more about the `TPOTClassifier` and `TPOTRegressor` classes in the [API documentation](/tpot/api/).

Some example code with custom TPOT parameters might look like:

```Python
pipeline_optimizer = TPOTClassifier(generations=5, population_size=20, cv=5,
                                    random_state=42, verbosity=2)
```

Now TPOT is ready to optimize a pipeline for you. You can tell TPOT to optimize a pipeline based on a data set with the `fit` function:

```Python
pipeline_optimizer.fit(X_train, y_train)
```

The `fit` function initializes the genetic programming algorithm to find the highest-scoring pipeline based on average k-fold cross-validation
Then, the pipeline is trained on the entire set of provided samples, and the TPOT instance can be used as a fitted model.

You can then proceed to evaluate the final pipeline on the testing set with the `score` function:

```Python
print(pipeline_optimizer.score(X_test, y_test))
```

Finally, you can tell TPOT to export the corresponding Python code for the optimized pipeline to a text file with the `export` function:

```Python
pipeline_optimizer.export('tpot_exported_pipeline.py')
```

Once this code finishes running, `tpot_exported_pipeline.py` will contain the Python code for the optimized pipeline.

Below is a full example script using TPOT to optimize a pipeline, score it, and export the best pipeline to a file.

```Python
from tpot import TPOTClassifier
from sklearn.datasets import load_digits
from sklearn.model_selection import train_test_split

digits = load_digits()
X_train, X_test, y_train, y_test = train_test_split(digits.data, digits.target,
                                                    train_size=0.75, test_size=0.25)

pipeline_optimizer = TPOTClassifier(generations=5, population_size=20, cv=5,
                                    random_state=42, verbosity=2)
pipeline_optimizer.fit(X_train, y_train)
print(pipeline_optimizer.score(X_test, y_test))
pipeline_optimizer.export('tpot_exported_pipeline.py')
```

Check our [examples](/tpot/examples/) to see TPOT applied to some specific data sets.

## TPOT on the command line

To use TPOT via the command line, enter the following command with a path to the data file:

```Shell
tpot /path_to/data_file.csv
```

An example command-line call to TPOT may look like:

```Shell
tpot data/mnist.csv -is , -target class -o tpot_exported_pipeline.py -g 5 -p 20 -cv 5 -s 42 -v 2
```

TPOT offers several arguments that can be provided at the command line. To see brief descriptions of these arguments,
enter the following command:

```Shell
tpot --help
```

Detailed descriptions of the command-line arguments are below.

<table>
<tr>
<th>Argument</th>
<th>Parameter</th>
<th width="15%">Valid values</th>
<th>Effect</th>
</tr>
<tr>
<td>-is</td>
<td>INPUT_SEPARATOR</td>
<td>Any string</td>
<td>Character used to separate columns in the input file.</td>
</tr>
<tr>
<td>-target</td>
<td>TARGET_NAME</td>
<td>Any string</td>
<td>Name of the target column in the input file.</td>
</tr>
<tr>
<td>-mode</td>
<td>TPOT_MODE</td>
<td>['classification', 'regression']</td>
<td>Whether TPOT is being used for a supervised classification or regression problem.</td>
</tr>
<tr>
<td>-o</td>
<td>OUTPUT_FILE</td>
<td>String path to a file</td>
<td>File to export the code for the final optimized pipeline.</td>
</tr>
<tr>
<td>-g</td>
<td>GENERATIONS</td>
<td>Any positive integer or None</td>
<td>Number of iterations to run the pipeline optimization process. It must be a positive number or None. If None, the parameter max_time_mins must be defined as the runtime limit. Generally, TPOT will work better when you give it more generations (and therefore time) to optimize the pipeline.
<br /><br />
TPOT will evaluate POPULATION_SIZE + GENERATIONS x OFFSPRING_SIZE pipelines in total.</td>
</tr>
<tr>
<td>-p</td>
<td>POPULATION_SIZE</td>
<td>Any positive integer</td>
<td>Number of individuals to retain in the GP population every generation. Generally, TPOT will work better when you give it more individuals (and therefore time) to optimize the pipeline.
<br /><br />
TPOT will evaluate POPULATION_SIZE + GENERATIONS x OFFSPRING_SIZE pipelines in total.</td>
</tr>
<tr>
<td>-os</td>
<td>OFFSPRING_SIZE</td>
<td>Any positive integer</td>
<td>Number of offspring to produce in each GP generation.
<br /><br />
By default, OFFSPRING_SIZE = POPULATION_SIZE.</td>
</tr>
<tr>
<td>-mr</td>
<td>MUTATION_RATE</td>
<td>[0.0, 1.0]</td>
<td>GP mutation rate in the range [0.0, 1.0]. This tells the GP algorithm how many pipelines to apply random changes to every generation.
<br /><br />
We recommend using the default parameter unless you understand how the mutation rate affects GP algorithms.</td>
</tr>
<tr>
<td>-xr</td>
<td>CROSSOVER_RATE</td>
<td>[0.0, 1.0]</td>
<td>GP crossover rate in the range [0.0, 1.0]. This tells the GP algorithm how many pipelines to "breed" every generation.
<br /><br />
We recommend using the default parameter unless you understand how the crossover rate affects GP algorithms.</td>
</tr>
<tr>
<td>-scoring</td>
<td>SCORING_FN</td>
<td>'accuracy', 'adjusted_rand_score', 'average_precision', 'balanced_accuracy',<br />'f1',
'f1_macro', 'f1_micro', 'f1_samples', 'f1_weighted', 'neg_log_loss', 'neg_mean_absolute_error',
'neg_mean_squared_error', 'neg_median_absolute_error', 'precision', 'precision_macro', 'precision_micro',
'precision_samples', 'precision_weighted',<br />'r2', 'recall', 'recall_macro', 'recall_micro', 'recall_samples',
'recall_weighted', 'roc_auc', 'my_module.scorer_name*'</td>
<td>Function used to evaluate the quality of a given pipeline for the problem. By default, accuracy is used for classification and mean squared error (MSE) is used for regression.
<br /><br />
TPOT assumes that any function with "error" or "loss" in the name is meant to be minimized, whereas any other functions will be maximized.
<br /><br />
my_module.scorer_name: You can also specify your own function or a full python path to an existing one.
<br /><br />
See the section on <a href="#scoring-functions">scoring functions</a> for more details.</td>
</tr>
<tr>
<td>-cv</td>
<td>CV</td>
<td>Any integer > 1</td>
<td>Number of folds to evaluate each pipeline over in k-fold cross-validation during the TPOT optimization process.</td>
</tr>
<td>-sub</td>
<td>SUBSAMPLE</td>
<td>(0.0, 1.0]</td>
<td>Subsample ratio of the training instance. Setting it to 0.5 means that TPOT randomly collects half of training samples for pipeline optimization process.</td>
</tr>
<tr>
<td>-njobs</td>
<td>NUM_JOBS</td>
<td>Any positive integer or -1</td>
<td>Number of CPUs for evaluating pipelines in parallel during the TPOT optimization process.
<br /><br />
Assigning this to -1 will use as many cores as available on the computer. For n_jobs below -1, (n_cpus + 1 + n_jobs) are used. Thus for n_jobs = -2, all CPUs but one are used.</td>
</tr>
<tr>
<td>-maxtime</td>
<td>MAX_TIME_MINS</td>
<td>Any positive integer</td>
<td>How many minutes TPOT has to optimize the pipeline.
<br /><br />
How many minutes TPOT has to optimize the pipeline.If not None, this setting will allow TPOT to run until max_time_mins minutes elapsed and then stop. TPOT will stop earlier if generationsis set and all generations are already evaluated.</td>
</tr>
<tr>
<td>-maxeval</td>
<td>MAX_EVAL_MINS</td>
<td>Any positive float</td>
<td>How many minutes TPOT has to evaluate a single pipeline.
<br /><br />
Setting this parameter to higher values will allow TPOT to consider more complex pipelines but will also allow TPOT to run longer.</td>
</tr>
<tr>
<td>-s</td>
<td>RANDOM_STATE</td>
<td>Any positive integer</td>
<td>Random number generator seed for reproducibility.
<br /><br />
Set this seed if you want your TPOT run to be reproducible with the same seed and data set in the future.</td>
</tr>
<tr>
<td>-config</td>
<td>CONFIG_FILE</td>
<td>String or file path</td>
<td>Operators and parameter configurations in TPOT:
<br /><br />
<ul>
<li>Path for configuration file: TPOT will use the path to a configuration file for customizing the operators and parameters that TPOT uses in the optimization process</li>
<li>string 'TPOT light', TPOT will use a built-in configuration with only fast models and preprocessors</li>
<li>string 'TPOT MDR', TPOT will use a built-in configuration specialized for genomic studies</li>
<li>string 'TPOT sparse': TPOT will use a configuration dictionary with a one-hot encoder and the operators normally included in TPOT that also support sparse matrices.</li>
</ul>
See the <a href="../using/#built-in-tpot-configurations">built-in configurations</a> section for the list of configurations included with TPOT, and the <a href="../using/#customizing-tpots-operators-and-parameters">custom configuration</a> section for more information and examples of how to create your own TPOT configurations.
</td>
</tr>
<tr>
<td>-template</td>
<td>TEMPLATE</td>
<td>String</td>
<td>Template of predefined pipeline structure. The option is for specifying a desired structure for the machine learning pipeline evaluated in TPOT. So far this option only supports linear pipeline structure. Each step in the pipeline should be a main class of operators (Selector, Transformer, Classifier or Regressor) or a specific operator (e.g. `SelectPercentile`) defined in TPOT operator configuration. If one step is a main class, TPOT will randomly assign all subclass operators (subclasses of [`SelectorMixin`](https://github.com/scikit-learn/scikit-learn/blob/master/sklearn/feature_selection/base.py#L17), [`TransformerMixin`](https://scikit-learn.org/stable/modules/generated/sklearn.base.TransformerMixin.html), [`ClassifierMixin`](https://scikit-learn.org/stable/modules/generated/sklearn.base.ClassifierMixin.html) or [`RegressorMixin`](https://scikit-learn.org/stable/modules/generated/sklearn.base.RegressorMixin.html) in scikit-learn) to that step. Steps in the template are delimited by "-", e.g. "SelectPercentile-Transformer-Classifier". By default value of template is None, TPOT generates tree-based pipeline randomly.

See the <a href="../using/#template-option-in-tpot"> template option in tpot</a> section for more details.
</td>
</tr>
<tr>
<td>-memory</td>
<td>MEMORY</td>
<td>String or file path</td>
<td>If supplied, pipeline will cache each transformer after calling fit. This feature is used to avoid computing the fit transformers within a pipeline if the parameters and input data are identical with another fitted pipeline during optimization process. Memory caching mode in TPOT:
<br /><br />
<ul>
<li>Path for a caching directory: TPOT uses memory caching with the provided directory and TPOT does NOT clean the caching directory up upon shutdown.</li>
<li>string 'auto': TPOT uses memory caching with a temporary directory and cleans it up upon shutdown.</li>
</ul>
</td>
</tr>
<tr>
<td>-cf</td>
<td>CHECKPOINT_FOLDER</td>
<td>Folder path</td>
<td>
If supplied, a folder you created, in which tpot will periodically save pipelines in pareto front so far while optimizing.
<br /><br />
This is useful in multiple cases:
<ul>
<li>sudden death before tpot could save an optimized pipeline</li>
<li>progress tracking</li>
<li>grabbing a pipeline while tpot is working</li>
</ul>
<br /><br />
Example:
<br />
mkdir my_checkpoints
<br />
-cf ./my_checkpoints
</tr>
<tr>
<td>-es</td>
<td>EARLY_STOP</td>
<td>Any positive integer</td>
<td>
How many generations TPOT checks whether there is no improvement in optimization process.
<br /><br />
End optimization process if there is no improvement in the set number of generations.
</tr>
<tr>
<td>-v</td>
<td>VERBOSITY</td>
<td>{0, 1, 2, 3}</td>
<td>How much information TPOT communicates while it is running.
<br /><br />
0 = none, 1 = minimal, 2 = high, 3 = all.
<br /><br />
A setting of 2 or higher will add a progress bar during the optimization procedure.</td>
</tr>
<tr>
<td>-log</td>
<td>LOG</td>
<td>Folder path</td>
<td>Save progress content to a file.</td>
</tr>
<tr>
<td colspan=3>--no-update-check</td>
<td>Flag indicating whether the TPOT version checker should be disabled.</td>
</tr>
<tr>
<td colspan=3>--version</td>
<td>Show TPOT's version number and exit.</td>
</tr>
<tr>
<td colspan=3>--help</td>
<td>Show TPOT's help documentation and exit.</td>
</tr>
</table>

## Scoring functions

TPOT makes use of `sklearn.model_selection.cross_val_score` for evaluating pipelines, and as such offers the same support for scoring functions. There are two ways to make use of scoring functions with TPOT:

- You can pass in a string to the `scoring` parameter from the list above. Any other strings will cause TPOT to throw an exception.

- You can pass the callable object/function with signature `scorer(estimator, X, y)`, where `estimator` is trained estimator to use for scoring, `X` are features that will be passed to `estimator.predict` and `y` are target values for `X`. To do this, you should implement your own function. See the example below for further explanation.

```Python
from tpot import TPOTClassifier
from sklearn.datasets import load_digits
from sklearn.model_selection import train_test_split
from sklearn.metrics import make_scorer

digits = load_digits()
X_train, X_test, y_train, y_test = train_test_split(digits.data, digits.target,
                                                    train_size=0.75, test_size=0.25)
# Make a custom metric function
def my_custom_accuracy(y_true, y_pred):
    return float(sum(y_pred == y_true)) / len(y_true)

# Make a custom a scorer from the custom metric function
# Note: greater_is_better=False in make_scorer below would mean that the scoring function should be minimized.
my_custom_scorer = make_scorer(my_custom_accuracy, greater_is_better=True)

tpot = TPOTClassifier(generations=5, population_size=20, verbosity=2,
                      scoring=my_custom_scorer)
tpot.fit(X_train, y_train)
print(tpot.score(X_test, y_test))
tpot.export('tpot_digits_pipeline.py')
```

* **my_module.scorer_name**: You can also use a custom `score_func(y_true, y_pred)` or `scorer(estimator, X, y)` function through the command line by adding the argument `-scoring my_module.scorer` to your command-line call. TPOT will import your module and use the custom scoring function from there. TPOT will include your current working directory when importing the module, so you can place it in the same directory where you are going to run TPOT.
Example: `-scoring sklearn.metrics.auc` will use the function auc from sklearn.metrics module.

## Built-in TPOT configurations

TPOT comes with a handful of default operators and parameter configurations that we believe work well for optimizing machine learning pipelines. Below is a list of the current built-in configurations that come with TPOT.

<table>
<tr>
<th align="left">Configuration Name</th>
<th align="left">Description</th>
<th align="left">Operators</th>
</tr>

<tr>
<td>Default TPOT</td>
<td>TPOT will search over a broad range of preprocessors, feature constructors, feature selectors, models, and parameters to find a series of operators that minimize the error of the model predictions. Some of these operators are complex and may take a long time to run, especially on larger datasets.
<br /><br />
<strong>Note: This is the default configuration for TPOT.</strong> To use this configuration, use the default value (None) for the config_dict parameter.</td>
<td align="center"><a href="https://github.com/EpistasisLab/tpot/blob/master/tpot/config/classifier.py">Classification</a>
<br /><br />
<a href="https://github.com/EpistasisLab/tpot/blob/master/tpot/config/regressor.py">Regression</a></td>
</tr>

<tr>
<td>TPOT light</td>
<td>TPOT will search over a restricted range of preprocessors, feature constructors, feature selectors, models, and parameters to find a series of operators that minimize the error of the model predictions. Only simpler and fast-running operators will be used in these pipelines, so TPOT light is useful for finding quick and simple pipelines for a classification or regression problem.
<br /><br />
This configuration works for both the TPOTClassifier and TPOTRegressor.</td>
<td align="center"><a href="https://github.com/EpistasisLab/tpot/blob/master/tpot/config/classifier_light.py">Classification</a>
<br /><br />
<a href="https://github.com/EpistasisLab/tpot/blob/master/tpot/config/regressor_light.py">Regression</a></td>
</tr>

<tr>
<td>TPOT MDR</td>
<td>TPOT will search over a series of feature selectors and <a href="https://en.wikipedia.org/wiki/Multifactor_dimensionality_reduction">Multifactor Dimensionality Reduction</a> models to find a series of operators that maximize prediction accuracy. The TPOT MDR configuration is specialized for <a href="https://en.wikipedia.org/wiki/Genome-wide_association_study">genome-wide association studies (GWAS)</a>, and is described in detail online <a href="https://arxiv.org/abs/1702.01780">here</a>.
<br /><br />
Note that TPOT MDR may be slow to run because the feature selection routines are computationally expensive, especially on large datasets.</td>
<td align="center"><a href="https://github.com/EpistasisLab/tpot/blob/master/tpot/config/classifier_mdr.py">Classification</a>
<br /><br />
<a href="https://github.com/EpistasisLab/tpot/blob/master/tpot/config/regressor_mdr.py">Regression</a></td>
</tr>

<tr>
<td>TPOT sparse</td>
<td>TPOT uses a configuration dictionary with a one-hot encoder and the operators normally included in TPOT that also support sparse matrices.
<br /><br />
This configuration works for both the TPOTClassifier and TPOTRegressor.</td>
<td align="center"><a href="https://github.com/EpistasisLab/tpot/blob/master/tpot/config/classifier_sparse.py">Classification</a>
<br /><br />
<a href="https://github.com/EpistasisLab/tpot/blob/master/tpot/config/regressor_sparse.py">Regression</a></td>
</tr>

<tr>
<td>TPOT NN</td>
<td>TPOT uses the same configuration as "Default TPOT" plus additional neural network estimators written in PyTorch (currently only `tpot.builtins.PytorchLRClassifier` and `tpot.builtins.PytorchMLPClassifier`).
<br /><br />
Currently only classification is supported, but future releases will include regression estimators.</td>
<td align="center"><a href="https://github.com/EpistasisLab/tpot/blob/master/tpot/config/classifier_nn.py">Classification</a></td>
</tr>

<tr>
<td>TPOT cuML</td>
<td>TPOT will search over a restricted configuration using the GPU-accelerated estimators in <a href="https://github.com/rapidsai/cuml">RAPIDS cuML</a> and <a href="https://github.com/dmlc/xgboost">DMLC XGBoost</a>. This configuration requires an NVIDIA Pascal architecture or better GPU with compute capability 6.0+, and that the library cuML is installed. With this configuration, all model training and predicting will be GPU-accelerated.
<br /><br />
This configuration is particularly useful for medium-sized and larger datasets on which CPU-based estimators are a common bottleneck, and works for both the TPOTClassifier and TPOTRegressor.</td>
<td align="center"><a href="https://github.com/EpistasisLab/tpot/blob/master/tpot/config/classifier_cuml.py">Classification</a>
<br /><br />
<a href="https://github.com/EpistasisLab/tpot/blob/master/tpot/config/regressor_cuml.py">Regression</a></td>
</tr>

</table>

To use any of these configurations, simply pass the string name of the configuration to the `config_dict` parameter (or `-config` on the command line). For example, to use the "TPOT light" configuration:

```Python
from tpot import TPOTClassifier
from sklearn.datasets import load_digits
from sklearn.model_selection import train_test_split

digits = load_digits()
X_train, X_test, y_train, y_test = train_test_split(digits.data, digits.target,
                                                    train_size=0.75, test_size=0.25)

tpot = TPOTClassifier(generations=5, population_size=20, verbosity=2,
                      config_dict='TPOT light')
tpot.fit(X_train, y_train)
print(tpot.score(X_test, y_test))
tpot.export('tpot_digits_pipeline.py')

```

## Customizing TPOT's operators and parameters

Beyond the default configurations that come with TPOT, in some cases it is useful to limit the algorithms and parameters that TPOT considers. For that reason, we allow users to provide TPOT with a custom configuration for its operators and parameters.

The custom TPOT configuration must be in nested dictionary format, where the first level key is the path and name of the operator (e.g., `sklearn.naive_bayes.MultinomialNB`) and the second level key is the corresponding parameter name for that operator (e.g., `fit_prior`). The second level key should point to a list of parameter values for that parameter, e.g., `'fit_prior': [True, False]`.

For a simple example, the configuration could be:

```Python
tpot_config = {
    'sklearn.naive_bayes.GaussianNB': {
    },

    'sklearn.naive_bayes.BernoulliNB': {
        'alpha': [1e-3, 1e-2, 1e-1, 1., 10., 100.],
        'fit_prior': [True, False]
    },

    'sklearn.naive_bayes.MultinomialNB': {
        'alpha': [1e-3, 1e-2, 1e-1, 1., 10., 100.],
        'fit_prior': [True, False]
    }
}
```

in which case TPOT would only consider pipelines containing `GaussianNB`, `BernoulliNB`, `MultinomialNB`, and tune those algorithm's parameters in the ranges provided. This dictionary can be passed directly within the code to the `TPOTClassifier`/`TPOTRegressor` `config_dict` parameter, described above. For example:

```Python
from tpot import TPOTClassifier
from sklearn.datasets import load_digits
from sklearn.model_selection import train_test_split

digits = load_digits()
X_train, X_test, y_train, y_test = train_test_split(digits.data, digits.target,
                                                    train_size=0.75, test_size=0.25)

tpot_config = {
    'sklearn.naive_bayes.GaussianNB': {
    },

    'sklearn.naive_bayes.BernoulliNB': {
        'alpha': [1e-3, 1e-2, 1e-1, 1., 10., 100.],
        'fit_prior': [True, False]
    },

    'sklearn.naive_bayes.MultinomialNB': {
        'alpha': [1e-3, 1e-2, 1e-1, 1., 10., 100.],
        'fit_prior': [True, False]
    }
}

tpot = TPOTClassifier(generations=5, population_size=20, verbosity=2,
                      config_dict=tpot_config)
tpot.fit(X_train, y_train)
print(tpot.score(X_test, y_test))
tpot.export('tpot_digits_pipeline.py')
```

Command-line users must create a separate `.py` file with the custom configuration and provide the path to the file to the `tpot` call. For example, if the simple example configuration above is saved in `tpot_classifier_config.py`, that configuration could be used on the command line with the command:

```
tpot data/mnist.csv -is , -target class -config tpot_classifier_config.py -g 5 -p 20 -v 2 -o tpot_exported_pipeline.py
```

When using the command-line interface, the configuration file specified in the `-config` parameter *must* name its custom TPOT configuration `tpot_config`. Otherwise, TPOT will not be able to locate the configuration dictionary.

For more detailed examples of how to customize TPOT's operator configuration, see the default configurations for [classification](https://github.com/EpistasisLab/tpot/blob/master/tpot/config/classifier.py) and [regression](https://github.com/EpistasisLab/tpot/blob/master/tpot/config/regressor.py) in TPOT's source code.

Note that you must have all of the corresponding packages for the operators installed on your computer, otherwise TPOT will not be able to use them. For example, if XGBoost is not installed on your computer, then TPOT will simply not import nor use XGBoost in the pipelines it considers.


## Template option in TPOT

Template option provides a way to specify a desired structure for machine learning pipeline, which may reduce TPOT computation time and potentially provide more interpretable results. Current implementation only supports linear pipelines.

Below is a simple example to use `template` option. The pipelines generated/evaluated in TPOT will follow this structure: 1st step is a feature selector (a subclass of [`SelectorMixin`](https://github.com/scikit-learn/scikit-learn/blob/master/sklearn/feature_selection/base.py#L17)), 2nd step is a feature transformer (a subclass of [`TransformerMixin`](https://scikit-learn.org/stable/modules/generated/sklearn.base.TransformerMixin.html)) and 3rd step is a classifier for classification (a subclass of [`ClassifierMixin`](https://scikit-learn.org/stable/modules/generated/sklearn.base.ClassifierMixin.html)). The last step must be `Classifier` for `TPOTClassifier`'s template but `Regressor` for `TPOTRegressor`. **Note: although `SelectorMixin` is subclass of `TransformerMixin` in scikit-learn, but `Transformer` in this option excludes those subclasses of `SelectorMixin`.**

```Python
tpot_obj = TPOTClassifier(
                template='Selector-Transformer-Classifier'
                )
```

If a specific operator, e.g. `SelectPercentile`, is preferred for usage in the 1st step of the pipeline, the template can be defined like 'SelectPercentile-Transformer-Classifier'.


## FeatureSetSelector in TPOT

`FeatureSetSelector` is a special new operator in TPOT. This operator enables feature selection based on *priori* expert knowledge. For example, in RNA-seq gene expression analysis, this operator can be used to select one or more gene (feature) set(s) based on GO (Gene Ontology) terms or annotated gene sets Molecular Signatures Database ([MSigDB](http://software.broadinstitute.org/gsea/msigdb/index.jsp)) in the 1st step of pipeline via `template` option above, in order to reduce dimensions and TPOT computation time. This operator requires a dataset list in csv format. In this csv file, there are only three columns: 1st column is feature set names, 2nd column is the total number of features in one set and 3rd column is a list of feature names (if input X is pandas.DataFrame) or indexes (if input X is numpy.ndarray) delimited by ";". Below is an example how to use this operator in TPOT.

Please check our [preprint paper](https://www.biorxiv.org/content/10.1101/502484v1.article-info) for more details.

```Python
from tpot import TPOTClassifier
import numpy as np
import pandas as pd
from tpot.config import classifier_config_dict
test_data = pd.read_csv("https://raw.githubusercontent.com/EpistasisLab/tpot/master/tests/tests.csv")
test_X = test_data.drop("class", axis=1)
test_y = test_data['class']

# add FeatureSetSelector into tpot configuration
classifier_config_dict['tpot.builtins.FeatureSetSelector'] = {
    'subset_list': ['https://raw.githubusercontent.com/EpistasisLab/tpot/master/tests/subset_test.csv'],
    'sel_subset': [0,1] # select only one feature set, a list of index of subset in the list above
    #'sel_subset': list(combinations(range(3), 2)) # select two feature sets
}


tpot = TPOTClassifier(generations=5,
                           population_size=50, verbosity=2,
                           template='FeatureSetSelector-Transformer-Classifier',
                           config_dict=classifier_config_dict)
tpot.fit(test_X, test_y)
```

## Pipeline caching in TPOT

With the `memory` parameter, pipelines can cache the results of each transformer after fitting them. This feature is used to avoid repeated computation by transformers within a pipeline if the parameters and input data are identical to another fitted pipeline during optimization process. TPOT allows users to specify a custom directory path or [`joblib.Memory`](https://joblib.readthedocs.io/en/latest/generated/joblib.Memory.html) in case they want to re-use the memory cache in future TPOT runs (or a `warm_start` run).

There are three methods for enabling memory caching in TPOT:

```Python
from tpot import TPOTClassifier
from tempfile import mkdtemp
from joblib import Memory
from shutil import rmtree

# Method 1, auto mode: TPOT uses memory caching with a temporary directory and cleans it up upon shutdown
tpot = TPOTClassifier(memory='auto')

# Method 2, with a custom directory for memory caching
tpot = TPOTClassifier(memory='/to/your/path')

# Method 3, with a Memory object
cachedir = mkdtemp() # Create a temporary folder
memory = Memory(cachedir=cachedir, verbose=0)
tpot = TPOTClassifier(memory=memory)

# Clear the cache directory when you don't need it anymore
rmtree(cachedir)
```

**Note: TPOT does NOT clean up memory caches if users set a custom directory path or Memory object. We recommend that you clean up the memory caches when you don't need it anymore.**

## Crash/freeze issue with n_jobs > 1 under OSX or Linux

Internally, TPOT uses [joblib](http://joblib.readthedocs.io/) to fit estimators in parallel.
This is the same parallelization framework used by scikit-learn. But it may crash/freeze with n_jobs > 1 under OSX or Linux [as scikit-learn does](http://scikit-learn.org/stable/faq.html#why-do-i-sometime-get-a-crash-freeze-with-n-jobs-1-under-osx-or-linux), especially with large datasets.

One solution is to configure Python's `multiprocessing` module to use the `forkserver` start method (instead of the default `fork`) to manage the process pools. You can enable the `forkserver` mode globally for your program by putting the following codes into your main script:

```Python
import multiprocessing

# other imports, custom code, load data, define model...

if __name__ == '__main__':
    multiprocessing.set_start_method('forkserver')

    # call scikit-learn utils or tpot utils with n_jobs > 1 here
```

More information about these start methods can be found in the [multiprocessing documentation](https://docs.python.org/3/library/multiprocessing.html#contexts-and-start-methods).

## Parallel Training with Dask

For large problems or working on Jupyter notebook, we highly recommend that you can distribute the work on a [Dask](http://dask.pydata.org/en/latest/) cluster.
The [dask-examples binder](https://mybinder.org/v2/gh/dask/dask-examples/master?filepath=machine-learning%2Ftpot.ipynb) has a runnable example
with a small dask cluster.

To use your Dask cluster to fit a TPOT model, specify the ``use_dask`` keyword when you create the TPOT estimator. **Note: if `use_dask=True`, TPOT will use as many cores as available on the your Dask cluster. If `n_jobs` is specified, then it will control the chunk size (10*`n_jobs` if it is less then offspring size) of parallel training.**

```python
estimator = TPOTEstimator(use_dask=True, n_jobs=-1)
```

This will use all the workers on your cluster to do the training, and use [Dask-ML's pipeline rewriting](https://dask-ml.readthedocs.io/en/latest/hyper-parameter-search.html#avoid-repeated-work) to avoid re-fitting estimators multiple times on the same set of data.
It will also provide fine-grained diagnostics in the [distributed scheduler UI](https://distributed.readthedocs.io/en/latest/web.html).

Alternatively, Dask implements a joblib backend.
You can instruct TPOT to use the distributed backend during training by specifying a `joblib.parallel_backend`:

```python
import joblib
import distributed.joblib
from dask.distributed import Client

# connect to the cluster
client = Client('schedueler-address')

# create the estimator normally
estimator = TPOTClassifier(n_jobs=-1)

# perform the fit in this context manager
with joblib.parallel_backend("dask"):
    estimator.fit(X, y)
```

See [dask's distributed joblib integration](https://distributed.readthedocs.io/en/latest/joblib.html) for more.

## Neural Networks in TPOT (`tpot.nn`)

Support for neural network models and deep learning is an experimental feature newly added to TPOT. Available neural network architectures are provided by the `tpot.nn` module. Unlike regular `sklearn` estimators, these models need to be written by hand, and must also inherit the appropriate base classes provided by `sklearn` for all of their built-in modules. In other words, they need implement methods like `.fit()`, `fit_transform()`, `get_params()`, etc., as described in detail on [Developing scikit-learn estimators](https://scikit-learn.org/stable/developers/develop.html).

### Telling TPOT to use built-in PyTorch neural network models

Mainly due to the issues described below, TPOT won't use its neural network models unless you explicitly tell it to do so. This is done as follows:

- Use `import tpot.nn` before instantiating any TPOT estimators.

- Use a configuration dictionary that includes one or more `tpot.nn` estimators, either by writing one manually, including one from a file, or by importing the configuration in `tpot/config/classifier_nn.py`. A very simple example that will force TPOT to only use a PyTorch-based logistic regression classifier as its main estimator is as follows:

```python
tpot_config = {
    'tpot.nn.PytorchLRClassifier': {
        'learning_rate': [1e-3, 1e-2, 1e-1, 0.5, 1.]
    }
}
```

- Alternatively, use a template string including `PytorchLRClassifier` or `PytorchMLPClassifier` while loading the TPOT-NN configuration dictionary.

Neural network models are notorious for being extremely sensitive to their initialization parameters, so you may need to heavily adjust `tpot.nn` configuration dictionaries in order to attain good performance on your dataset.

A simple example of using TPOT-NN is shown in [examples](/tpot/examples/).

### Important caveats

- Neural network models (especially when they reach moderately large sizes) take a notoriously large amount of time and computing power to train. You should expect `tpot.nn` neural networks to train several orders of magnitude slower than their `sklearn` alternatives. This can be alleviated somewhat by training the models on computers with CUDA-enabled GPUs.

- TPOT will occasionally learn pipelines that stack several `sklearn` estimators. Mathematically, these can be nearly identical to some deep learning models. For example, by stacking several `sklearn.linear_model.LogisticRegression`s, you end up with a very close approximation of a Multilayer Perceptron; one of the simplest and most well known deep learning architectures. TPOT's genetic programming algorithms generally optimize these 'networks' much faster than PyTorch, which typically uses a more brute-force convex optimization approach.

- The problem of 'black box' model introspection is one of the most substantial criticisms and challenges of deep learning. This problem persists in `tpot.nn`, whereas TPOT's default estimators often are far easier to introspect.
# Overview

The following sections illustrate the usage of TPOT with various datasets, each
belonging to a typical class of machine learning tasks.

| Dataset | Task                    | Task class             | Dataset description | Jupyter notebook                                                                           |
| ------- | ----------------------- | ---------------------- |:-------------------:|:------------------------------------------------------------------------------------------:|
| Iris                  | flower classification   | classification         | [link](https://archive.ics.uci.edu/ml/datasets/iris) | [link](https://github.com/EpistasisLab/tpot/blob/master/tutorials/IRIS.ipynb) |
| Optical Recognition of Handwritten Digits                 | digit recognition       | (image) classification | [link](https://scikit-learn.org/stable/datasets/index.html#digits-dataset) | [link](https://github.com/EpistasisLab/tpot/blob/master/tutorials/Digits.ipynb) |
| Boston                | housing prices modeling | regression             | [link](https://www.cs.toronto.edu/~delve/data/boston/bostonDetail.html) | N/A    |
| Titanic               | survival analysis       | classification         | [link](https://www.kaggle.com/c/titanic/data) | [link](https://github.com/EpistasisLab/tpot/blob/master/tutorials/Titanic_Kaggle.ipynb) |
| Bank Marketing        | subscription prediction | classification         | [link](https://archive.ics.uci.edu/ml/datasets/Bank+Marketing) | [link](https://github.com/EpistasisLab/tpot/blob/master/tutorials/Portuguese%20Bank%20Marketing/Portuguese%20Bank%20Marketing%20Strategy.ipynb) |
| MAGIC Gamma Telescope | event detection         | classification         | [link](https://archive.ics.uci.edu/ml/datasets/MAGIC+Gamma+Telescope) | [link](https://github.com/EpistasisLab/tpot/blob/master/tutorials/MAGIC%20Gamma%20Telescope/MAGIC%20Gamma%20Telescope.ipynb) |
| cuML Classification Example | random classification problem         | classification         | [link](https://scikit-learn.org/stable/modules/generated/sklearn.datasets.make_classification.html) | [link](https://github.com/EpistasisLab/tpot/blob/master/tutorials/cuML_Classification_Example.ipynb) |
| cuML Regression Example | random regression problem         | regression         | [link](https://scikit-learn.org/stable/modules/generated/sklearn.datasets.make_regression.html) | [link](https://github.com/EpistasisLab/tpot/blob/master/tutorials/cuML_Regression_Example.ipynb) |

**Notes:**
- For details on how the `fit()`, `score()` and `export()` methods work, refer to the [usage documentation](/using/).
- Upon re-running the experiments, your resulting pipelines _may_ differ (to some extent) from the ones demonstrated here.

## Iris flower classification

The following code illustrates how TPOT can be employed for performing a simple _classification task_ over the Iris dataset.

```Python
from tpot import TPOTClassifier
from sklearn.datasets import load_iris
from sklearn.model_selection import train_test_split
import numpy as np

iris = load_iris()
X_train, X_test, y_train, y_test = train_test_split(iris.data.astype(np.float64),
    iris.target.astype(np.float64), train_size=0.75, test_size=0.25, random_state=42)

tpot = TPOTClassifier(generations=5, population_size=50, verbosity=2, random_state=42)
tpot.fit(X_train, y_train)
print(tpot.score(X_test, y_test))
tpot.export('tpot_iris_pipeline.py')
```

Running this code should discover a pipeline (exported as `tpot_iris_pipeline.py`) that achieves about 97% test accuracy:

```Python
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsClassifier
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import Normalizer
from tpot.export_utils import set_param_recursive

# NOTE: Make sure that the outcome column is labeled 'target' in the data file
tpot_data = pd.read_csv('PATH/TO/DATA/FILE', sep='COLUMN_SEPARATOR', dtype=np.float64)
features = tpot_data.drop('target', axis=1)
training_features, testing_features, training_target, testing_target = \
            train_test_split(features, tpot_data['target'], random_state=42)

# Average CV score on the training set was: 0.9826086956521738
exported_pipeline = make_pipeline(
    Normalizer(norm="l2"),
    KNeighborsClassifier(n_neighbors=5, p=2, weights="distance")
)
# Fix random state for all the steps in exported pipeline
set_param_recursive(exported_pipeline.steps, 'random_state', 42)

exported_pipeline.fit(training_features, training_target)
results = exported_pipeline.predict(testing_features)
```

## Digits dataset

Below is a minimal working example with the optical recognition of handwritten digits dataset, which is an _image classification problem_.

```Python
from tpot import TPOTClassifier
from sklearn.datasets import load_digits
from sklearn.model_selection import train_test_split

digits = load_digits()
X_train, X_test, y_train, y_test = train_test_split(digits.data, digits.target,
                                                    train_size=0.75, test_size=0.25, random_state=42)

tpot = TPOTClassifier(generations=5, population_size=50, verbosity=2, random_state=42)
tpot.fit(X_train, y_train)
print(tpot.score(X_test, y_test))
tpot.export('tpot_digits_pipeline.py')
```

Running this code should discover a pipeline (exported as `tpot_digits_pipeline.py`) that achieves about 98% test accuracy:

```Python
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline, make_union
from sklearn.preprocessing import PolynomialFeatures
from tpot.builtins import StackingEstimator
from tpot.export_utils import set_param_recursive

# NOTE: Make sure that the outcome column is labeled 'target' in the data file
tpot_data = pd.read_csv('PATH/TO/DATA/FILE', sep='COLUMN_SEPARATOR', dtype=np.float64)
features = tpot_data.drop('target', axis=1)
training_features, testing_features, training_target, testing_target = \
            train_test_split(features, tpot_data['target'], random_state=42)

# Average CV score on the training set was: 0.9799428471757372
exported_pipeline = make_pipeline(
    PolynomialFeatures(degree=2, include_bias=False, interaction_only=False),
    StackingEstimator(estimator=LogisticRegression(C=0.1, dual=False, penalty="l1")),
    RandomForestClassifier(bootstrap=True, criterion="entropy", max_features=0.35000000000000003, min_samples_leaf=20, min_samples_split=19, n_estimators=100)
)
# Fix random state for all the steps in exported pipeline
set_param_recursive(exported_pipeline.steps, 'random_state', 42)

exported_pipeline.fit(training_features, training_target)
results = exported_pipeline.predict(testing_features)
```

## Boston housing prices modeling

The following code illustrates how TPOT can be employed for performing a _regression task_ over the Boston housing prices dataset.

```Python
from tpot import TPOTRegressor
from sklearn.datasets import load_boston
from sklearn.model_selection import train_test_split

housing = load_boston()
X_train, X_test, y_train, y_test = train_test_split(housing.data, housing.target,
                                                    train_size=0.75, test_size=0.25, random_state=42)

tpot = TPOTRegressor(generations=5, population_size=50, verbosity=2, random_state=42)
tpot.fit(X_train, y_train)
print(tpot.score(X_test, y_test))
tpot.export('tpot_boston_pipeline.py')
```

Running this code should discover a pipeline (exported as `tpot_boston_pipeline.py`) that achieves at least 10 mean squared error (MSE) on the test set:

```Python
import numpy as np
import pandas as pd
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import PolynomialFeatures
from tpot.export_utils import set_param_recursive

# NOTE: Make sure that the outcome column is labeled 'target' in the data file
tpot_data = pd.read_csv('PATH/TO/DATA/FILE', sep='COLUMN_SEPARATOR', dtype=np.float64)
features = tpot_data.drop('target', axis=1)
training_features, testing_features, training_target, testing_target = \
            train_test_split(features, tpot_data['target'], random_state=42)

# Average CV score on the training set was: -10.812040755234403
exported_pipeline = make_pipeline(
    PolynomialFeatures(degree=2, include_bias=False, interaction_only=False),
    ExtraTreesRegressor(bootstrap=False, max_features=0.5, min_samples_leaf=2, min_samples_split=3, n_estimators=100)
)
# Fix random state for all the steps in exported pipeline
set_param_recursive(exported_pipeline.steps, 'random_state', 42)

exported_pipeline.fit(training_features, training_target)
results = exported_pipeline.predict(testing_features)
```

## Titanic survival analysis

To see the TPOT applied the Titanic Kaggle dataset, see the Jupyter notebook [here](https://github.com/EpistasisLab/tpot/blob/master/tutorials/Titanic_Kaggle.ipynb). This example shows how to take a messy dataset and preprocess it such that it can be used in scikit-learn and TPOT.

## Portuguese Bank Marketing

The corresponding Jupyter notebook, containing the associated data preprocessing and analysis, can be found [here](https://github.com/EpistasisLab/tpot/blob/master/tutorials/Portuguese%20Bank%20Marketing/Portuguese%20Bank%20Marketing%20Stratergy.ipynb).

## MAGIC Gamma Telescope
The corresponding Jupyter notebook, containing the associated data preprocessing and analysis, can be found [here](https://github.com/EpistasisLab/tpot/blob/master/tutorials/MAGIC%20Gamma%20Telescope/MAGIC%20Gamma%20Telescope.ipynb).

## Neural network classifier using TPOT-NN
By loading the <a href="https://github.com/EpistasisLab/tpot/blob/master/tpot/config/classifier_nn.py">TPOT-NN configuration dictionary</a>, PyTorch estimators will be included for classification. Users can also create their own NN configuration dictionary that includes `tpot.builtins.PytorchLRClassifier` and/or `tpot.builtins.PytorchMLPClassifier`, or they can specify them using a template string, as shown in the following example:

```Python
from tpot import TPOTClassifier
from sklearn.datasets import make_blobs
from sklearn.model_selection import train_test_split

X, y = make_blobs(n_samples=100, centers=2, n_features=3, random_state=42)
X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=0.75, test_size=0.25)

clf = TPOTClassifier(config_dict='TPOT NN', template='Selector-Transformer-PytorchLRClassifier',
                     verbosity=2, population_size=10, generations=10)
clf.fit(X_train, y_train)
print(clf.score(X_test, y_test))
clf.export('tpot_nn_demo_pipeline.py')
```

This example is somewhat trivial, but it should result in nearly 100% classification accuracy.
