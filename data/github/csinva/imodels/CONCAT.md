<p align="center">
	<img align="center" width=60% src="https://csinva.io/imodels/img/imodels_logo.svg?sanitize=True&kill_cache=1"> </img>	 <br/>
	Python package for concise, transparent, and accurate predictive modeling. All sklearn-compatible and easy to use.
</p>
<p align="center">
  <a href="https://csinva.github.io/imodels/">docs</a> •
  <a href="#imodels-overview">imodels overview</a> •
  <a href="#demo-notebooks">demo notebooks</a>
</p>
<p align="center">
  <img src="https://img.shields.io/badge/license-mit-blue.svg">
  <img src="https://img.shields.io/badge/python-3.6--3.9-blue">
  <a href="https://doi.org/10.21105/joss.03192"><img src="https://joss.theoj.org/papers/10.21105/joss.03192/status.svg"></a>
  <a href="https://github.com/csinva/imodels/actions"><img src="https://github.com/csinva/imodels/workflows/tests/badge.svg"></a>
  <!--img src="https://img.shields.io/github/checks-status/csinva/imodels/master"-->
  <img src="https://img.shields.io/pypi/v/imodels?color=orange">
  <img src="https://static.pepy.tech/personalized-badge/imodels?period=total&units=none&left_color=gray&right_color=orange&left_text=downloads&kill_cache=12">
</p>  




## imodels overview

<img align="center" width=100% src="https://csinva.io/imodels/img/anim.gif"> </img>

Modern machine-learning models are increasingly complex, often making them difficult to interpret. This package provides a simple interface for fitting and using state-of-the-art interpretable models, all compatible with scikit-learn. These models can often replace black-box models (e.g. random forests) with simpler models (e.g. rule lists) while improving interpretability and computational efficiency, all without sacrificing predictive accuracy! Simply import a classifier or regressor and use the `fit` and `predict` methods, same as standard scikit-learn models.

```python
from imodels import BoostedRulesClassifier, BayesianRuleListClassifier, GreedyRuleListClassifier, SkopeRulesClassifier # see more models below
from imodels import SLIMRegressor, RuleFitRegressor

model = BoostedRulesClassifier()  # initialize a model
model.fit(X_train, y_train)   # fit model
preds = model.predict(X_test) # discrete predictions: shape is (n_test, 1)
preds_proba = model.predict_proba(X_test) # predicted probabilities: shape is (n_test, n_classes)
print(model) # print the rule-based model

-----------------------------
# the model consists of the following 3 rules
# if X1 > 5: then 80.5% risk
# else if X2 > 5: then 40% risk
# else: 10% risk
```

### Installation
Install with `pip install imodels` (see [here](https://github.com/csinva/imodels/blob/master/docs/troubleshooting.md) for help). 

### Supported models

| Model                       | Reference                                                    | Description                                                  |
| :-------------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ |
| Rulefit rule set            | [🗂️](https://csinva.io/imodels/rule_set/rule_fit.html), [🔗](https://github.com/christophM/rulefit), [📄](http://statweb.stanford.edu/~jhf/ftp/RuleFit.pdf) | Fits a sparse linear model on rules extracted from decision trees |
| Skope rule set              | [🗂️](https://csinva.io/imodels/rule_set/skope_rules.html#imodels.rule_set.skope_rules.SkopeRulesClassifier), [🔗](https://github.com/scikit-learn-contrib/skope-rules) | Extracts rules from gradient-boosted trees, deduplicates them,<br/>then linearly combines them based on their OOB precision |
| Boosted rule set            | [🗂️](https://csinva.io/imodels/rule_set/boosted_rules.html), [🔗](https://github.com/jaimeps/adaboost-implementation), [📄](https://www.sciencedirect.com/science/article/pii/S002200009791504X) | Sequentially fits a set of rules with Adaboost           |
| Slipper rule set            | [🗂️](https://csinva.io/imodels/rule_set/slipper.html), ㅤㅤ [📄](https://www.aaai.org/Papers/AAAI/1999/AAAI99-049.pdf) | Sequentially learns a set of rules with SLIPPER            |
| Bayesian rule set           | [🗂️](https://csinva.io/imodels/rule_set/brs.html#imodels.rule_set.brs.BayesianRuleSetClassifier), [🔗](https://github.com/wangtongada/BOA), [📄](https://www.jmlr.org/papers/volume18/16-003/16-003.pdf) | Finds concise rule set with Bayesian sampling (slow)  |
| Optimal rule list           | [🗂️](https://csinva.io/imodels/rule_list/corels_wrapper.html#imodels.rule_list.corels_wrapper.OptimalRuleListClassifier), [🔗](https://github.com/corels/pycorels), [📄](https://www.jmlr.org/papers/volume18/17-716/17-716.pdf) | Fits rule list using global optimization for sparsity (CORELS) |
| Bayesian rule list          | [🗂️](https://csinva.io/imodels/rule_list/bayesian_rule_list/bayesian_rule_list.html#imodels.rule_list.bayesian_rule_list.bayesian_rule_list.BayesianRuleListClassifier), [🔗](https://github.com/tmadl/sklearn-expertsys), [📄](https://arxiv.org/abs/1602.08610) | Fits compact rule list distribution with Bayesian sampling (slow) |
| Greedy rule list            | [🗂️](https://csinva.io/imodels/rule_list/greedy_rule_list.html), [🔗](https://medium.com/@penggongting/implementing-decision-tree-from-scratch-in-python-c732e7c69aea) | Uses CART to fit a list (only a single path), rather than a tree |
| OneR rule list              | [🗂️](https://csinva.io/imodels/rule_list/one_r.html), ㅤㅤ[📄](https://link.springer.com/article/10.1023/A:1022631118932) | Fits rule list restricted to only one feature              |
| Optimal rule tree           | [🗂️](https://csinva.io/imodels/tree/gosdt/pygosdt.html#imodels.tree.gosdt.pygosdt.OptimalTreeClassifier), [🔗](https://github.com/Jimmy-Lin/GeneralizedOptimalSparseDecisionTrees), [📄](https://arxiv.org/abs/2006.08690) | Fits succinct tree using global optimization for sparsity (GOSDT) |
| Greedy rule tree            | [🗂️](https://csinva.io/imodels/tree/cart_wrapper.html), [🔗](https://scikit-learn.org/stable/modules/tree.html), [📄](https://www.taylorfrancis.com/books/mono/10.1201/9781315139470/classification-regression-trees-leo-breiman-jerome-friedman-richard-olshen-charles-stone) | Greedily fits tree using CART                              |
| C4.5 rule tree        | [🗂️](https://csinva.io/imodels/tree/c45_tree/c45_tree.html#imodels.tree.c45_tree.c45_tree.C45TreeClassifier), [🔗](https://github.com/RaczeQ/scikit-learn-C4.5-tree-classifier), [📄](https://link.springer.com/article/10.1007/BF00993309) | Greedily fits tree using C4.5                           |
| Iterative random<br/>forest | [🗂️](https://csinva.io/imodels/tree/iterative_random_forest/iterative_random_forest.html), [🔗](https://github.com/Yu-Group/iterative-Random-Forest), [📄](https://www.pnas.org/content/115/8/1943) | Repeatedly fit random forest, giving features with<br/>high importance a higher chance of being selected |
| Sparse integer<br/>linear model | [🗂️](https://csinva.io/imodels/algebraic/slim.html), ㅤㅤ[📄](https://link.springer.com/article/10.1007/s10994-015-5528-6) | Sparse linear model with integer coefficients                           |
| Tree sums | [🗂️](https://csinva.io/imodels/tree/figs.html#imodels.tree.figs), ㅤㅤ[📄]() | Sum of small trees with very few total rules (FIGS)                          |
| Hierarchical shrinkage<br/>wrapper | [🗂️](https://csinva.io/imodels/tree/hierarchical_shrinkage.html), ㅤㅤ[📄]() | Use regularization to improve trees |
| Distillation<br/>wrapper | [🗂️](https://csinva.io/imodels/util/distillation.html), ㅤㅤ[📄]() | Train a black-box model,<br/>then distill it into an interpretable model |
| More models                 | ⌛                                                            | (Coming soon!) Lightweight Rule Induction, MLRules, ... |

<p align="center">
Docs <a href="https://csinva.io/imodels/">🗂️</a>, Reference code implementation 🔗, Research paper 📄
</br>
</p>

<details>
<summary>Also see our <a href="https://csinva.io/imodels/util/data_util.html#imodels.util.data_util.get_clean_dataset">util functions</a> for downloading popular tabular datasets (e.g. compas).</summary>
These functions, in conjunction with <a href="https://github.com/csinva/imodels-data">imodels-data</a> and <a href="https://github.com/Yu-Group/imodels-experiments">imodels-experiments</a>, make it simple to download data and run experiments on new models.
</details>

<details>
<summary>Also see our <a href="https://csinva.io/imodels/util/explain_errors.html">simple function</a> for explaining classification errors.</summary>
Fit an interpretable model to explain a previous model's errors (ex. in <a href="https://github.com/csinva/imodels/blob/master/notebooks/error_detection_demo.ipynb">this notebook📓</a>).
</details>

<details>
<summary>Also see our <a href="https://csinva.io/imodels/discretization/index.html">fast and effective discretizers</a> for data preprocessing.</summary>
<table>
<thead>
<tr>
<th>Discretizer</th>
<th>Reference</th>
<th>Description</th>
</tr>
</thead>
<tbody>
<tr>
<td>MDLP</td>
<td><a href="https://csinva.io/imodels/discretization/mdlp.html#imodels.discretization.mdlp.MDLPDiscretizer">🗂️</a>, <a href="https://github.com/navicto/Discretization-MDLPC">🔗</a>, <a href="https://trs.jpl.nasa.gov/handle/2014/35171">📄</a></td>
<td>Discretize using entropy minimization heuristic</td>
</tr>
<tr>
<td>Simple</td>
<td><a href="https://csinva.io/imodels/discretization/simple.html#imodels.discretization.simple.SimpleDiscretizer">🗂️</a>, <a href="https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.KBinsDiscretizer.html">🔗</a></td>
<td>Simple KBins discretization</td>
</tr>
<tr>
<td>Random Forest</td>
<td><a href="https://csinva.io/imodels/discretization/discretizer.html#imodels.discretization.discretizer.RFDiscretizer">🗂️</a></td>
<td>Discretize into bins based on random forest split popularity</td>
</tr>
</tbody>
</table>
</details>



The final form of the above models takes one of the following forms, which aim to be simultaneously simple to understand and highly predictive:

|                           Rule set                           |                        Rule list                        |                        Rule tree                        |                       Algebraic models                       |
| :----------------------------------------------------------: | :-----------------------------------------------------: | :-----------------------------------------------------: | :----------------------------------------------------------: |
| <img src="https://csinva.io/imodels/img/rule_set.jpg" width="100%"> | <img src="https://csinva.io/imodels/img/rule_list.jpg"> | <img src="https://csinva.io/imodels/img/rule_tree.jpg"> | <img src="https://csinva.io/imodels/img/algebraic_models.jpg"> |

Different models and algorithms vary not only in their final form but also in different choices made during modeling. In particular, many models differ in the 3 steps given by the table below.

<details>
<summary>ex. RuleFit and SkopeRules</summary>
RuleFit and SkopeRules differ only in the way they prune rules: RuleFit uses a linear model whereas SkopeRules heuristically deduplicates rules sharing overlap.
</details>

<details>
<summary>ex. Bayesian rule lists and greedy rule lists</summary>
Bayesian rule lists and greedy rule lists differ in how they select rules; bayesian rule lists perform a global optimization over possible rule lists while Greedy rule lists pick splits sequentially to maximize a given criterion.
</details>

<details>
<summary>ex. FPSkope and SkopeRules</summary>
FPSkope and SkopeRules differ only in the way they generate candidate rules: FPSkope uses FPgrowth whereas SkopeRules extracts rules from decision trees.
</details>

See the docs for individual models for futher descriptions.

|                  Rule candidate generation                   |                       Rule selection                       |                Rule postprocessing|
| :----------------------------------------------------------: | :--------------------------------------------------------: | :-------------------------------------------------------: |
| <img src="https://csinva.io/imodels/img/rule_candidates.jpg"> | <img src="https://csinva.io/imodels/img/rule_overfit.jpg"> | <img src="https://csinva.io/imodels/img/rule_pruned.jpg"> |

The code here contains many useful and customizable functions for rule-based learning in the [util folder](https://csinva.io/imodels/util/index.html). This includes functions / classes for rule deduplication, rule screening, and converting between trees, rulesets, and neural networks.

## Demo notebooks

Demos are contained in the [notebooks](notebooks) folder.

<details>
<summary><a href="notebooks/imodels_demo.ipynb">imodels demo</a></summary>
Shows how to fit, predict, and visualize with different interpretable models
</details>

<details>
<summary><a href="https://colab.research.google.com/drive/1WfqvSjegygT7p0gyqiWpRpiwz2ePtiao#scrollTo=bLnLknIuoWtQ">imodels colab demo</a> <a href="https://colab.research.google.com/drive/1WfqvSjegygT7p0gyqiWpRpiwz2ePtiao#scrollTo=bLnLknIuoWtQ"> <img src="https://colab.research.google.com/assets/colab-badge.svg"></a></summary>
Shows how to fit, predict, and visualize with different interpretable models
</details>

<details>
<summary><a href="https://github.com/csinva/iai-clinical-decision-rule/blob/master/notebooks/05_fit_interpretable_models.ipynb">clinical decision rule notebook</a></summary>
Shows an example of using <code>imodels</code> for deriving a clinical decision rule
</details>

<details>
<summary>posthoc analysis</summary>
We also include some demos of posthoc analysis, which occurs after fitting models:
<a href="notebooks/posthoc_analysis.ipynb">posthoc.ipynb</a> shows different simple analyses to interpret a trained model and 
<a href="notebooks/uncertainty_analysis.ipynb">uncertainty.ipynb</a> contains basic code to get uncertainty estimates for a model
</details>



## Support for different tasks

Different models support different machine-learning tasks. Current support for different models is given below (each of these models can be imported directly from imodels (e.g. `from imodels import RuleFitClassifier`):

| Model                       |                    Binary classification                     |                          Regression                          | Notes |
| :-------------------------- | :----------------------------------------------------------: | :----------------------------------------------------------: | --------------------------- |
| Rulefit rule set            | [RuleFitClassifier](https://csinva.io/imodels/rule_set/rule_fit.html#imodels.rule_set.rule_fit.RuleFitClassifier) | [RuleFitRegressor](https://csinva.io/imodels/rule_set/rule_fit.html#imodels.rule_set.rule_fit.RuleFitRegressor) |  |
| Skope rule set              | [SkopeRulesClassifier](https://csinva.io/imodels/rule_set/slipper.html#imodels.rule_set.slipper.SlipperClassifier) |                                                              |  |
| Boosted rule set            | [BoostedRulesClassifier](https://csinva.io/imodels/rule_set/boosted_rules.html#imodels.rule_set.boosted_rules.BoostedRulesClassifier) |                                                              |  |
| SLIPPER rule set            | [SlipperClassifier](https://csinva.io/imodels/rule_set/slipper.html#imodels.rule_set.slipper.SlipperClassifier) |                                                              |  |
| Bayesian rule set           | [BayesianRuleSetClassifier](https://csinva.io/imodels/rule_set/brs.html#imodels.rule_set.brs.BayesianRuleSetClassifier) |                                                              | Fails for large problems |
| Optimal rule list (CORELS)  | [OptimalRuleListClassifier](https://csinva.io/imodels/rule_list/corels_wrapper.html#imodels.rule_list.corels_wrapper.OptimalRuleListClassifier) |                                                              | Requires [corels](https://pypi.org/project/corels/), fails for large problems |
| Bayesian rule list          | [BayesianRuleListClassifier](https://csinva.io/imodels/rule_list/bayesian_rule_list/bayesian_rule_list.html#imodels.rule_list.bayesian_rule_list.bayesian_rule_list.BayesianRuleListClassifier) |                                                              |  |
| Greedy rule list            | [GreedyRuleListClassifier](https://csinva.io/imodels/rule_list/greedy_rule_list.html#imodels.rule_list.greedy_rule_list.GreedyRuleListClassifier) |                                                              |  |
| OneR rule list              | [OneRClassifier](https://csinva.io/imodels/rule_list/one_r.html#imodels.rule_list.one_r.OneRClassifier) |                                                              |  |
| Optimal rule tree (GOSDT)   | [OptimalTreeClassifier](https://csinva.io/imodels/tree/gosdt/pygosdt.html#imodels.tree.gosdt.pygosdt.OptimalTreeClassifier) |                                                              | Requires [gosdt](https://pypi.org/project/gosdt/), fails for large problems |
| Greedy rule tree (CART)     | [GreedyTreeClassifier](https://csinva.io/imodels/tree/cart_wrapper.html#imodels.tree.cart_wrapper.GreedyTreeClassifier) |      [GreedyTreeRegressor](https://csinva.io/imodels/tree/cart_wrapper.html#imodels.tree.cart_wrapper.GreedyTreeRegressor)                                                        |  |
| C4.5 rule tree              | [C45TreeClassifier](https://csinva.io/imodels/tree/c45_tree/c45_tree.html#imodels.tree.c45_tree.c45_tree.C45TreeClassifier) |           |  |
| Iterative random forest     | [IRFClassifier](https://csinva.io/imodels/tree/iterative_random_forest/iterative_random_forest.html#imodels.tree.iterative_random_forest.iterative_random_forest.IRFClassifier)                                                             |                                                              | Requires [irf](https://pypi.org/project/irf/) |
| Sparse integer linear model | [SLIMClassifier](https://csinva.io/imodels/algebraic/slim.html#imodels.algebraic.slim.SLIMClassifier) | [SLIMRegressor](https://csinva.io/imodels/algebraic/slim.html#imodels.algebraic.slim.SLIMRegressor) | Requires extra dependencies for speed |
| Sapling Sums (FIGS) | [FIGSClassifier](https://csinva.io/imodels/tree/figs.html#imodels.tree.figs.FIGSClassifier) | [FIGSRegressor](https://csinva.io/imodels/tree/figs.html#imodels.tree.figs.FIGSRegressor) |                                                              |
| Hierarchical shrinkage | [HSTreeClassifierCV](https://csinva.io/imodels/tree/hierarchical_shrinkage.html#imodels.tree.hierarchical_shrinkage.HSTreeClassifierCV) | [HSTreeRegressorCV](https://csinva.io/imodels/tree/hierarchical_shrinkage.html#imodels.tree.hierarchical_shrinkage.HSTreeRegressorCV) | Wraps any sklearn tree-based model |
| Distillation |  | [DistilledRegressor](https://csinva.io/imodels/docs/util/distillation.html#imodels.util.distillation.DistilledRegressor) | Wraps any sklearn-compatible models |

## References

<details>
<summary>Readings</summary>
<ul>
  <li>Interpretable ML good quick overview: murdoch et al. 2019, <a href="https://arxiv.org/pdf/1901.04592.pdf">pdf</a></li>
	<li>Interpretable ML book: molnar 2019, <a href="https://christophm.github.io/interpretable-ml-book/">pdf</a></li>
	<li>Case for interpretable models rather than post-hoc explanation: rudin 2019, <a href="https://arxiv.org/pdf/1811.10154.pdf">pdf</a></li>
	<li>Review on evaluating interpretability: doshi-velez & kim 2017, <a href="https://arxiv.org/pdf/1702.08608.pdf">pdf</a></li>	
</ul>
</details>

<details>
<summary>Reference implementations (also linked above)</summary>
The code here heavily derives from the wonderful work of previous projects. We seek to to extract out, unify, and maintain key parts of these projects.
<ul>
  <li><a href="https://github.com/corels/pycorels">pycorels</a> - by <a href="https://github.com/fingoldin">@fingoldin</a> and the <a href="https://github.com/corels/corels">original CORELS team</a>
  <li><a href="https://github.com/tmadl/sklearn-expertsys">sklearn-expertsys</a> - by <a href="https://github.com/tmadl">@tmadl</a> and <a href="https://github.com/kenben">@kenben</a> based on original code by <a href="http://lethalletham.com/">Ben Letham</a></li>
  <li><a href="https://github.com/christophM/rulefit">rulefit</a> - by <a href="https://github.com/christophM">@christophM</a></li>
  <li><a href="https://github.com/scikit-learn-contrib/skope-rules">skope-rules</a> - by the <a href="https://github.com/scikit-learn-contrib/skope-rules/blob/master/AUTHORS.rst">skope-rules team</a> (including <a href="https://github.com/ngoix">@ngoix</a>, <a href="https://github.com/floriangardin">@floriangardin</a>, <a href="https://github.com/datajms">@datajms</a>, <a href="">Bibi Ndiaye</a>, <a href="">Ronan Gautier</a>)</li>
  <li><a href="https://github.com/wangtongada/BOA">boa</a> - by <a href="https://github.com/wangtongada">@wangtongada</a></li>	
</ul>
</details>

<details>
<summary>Related packages</summary>
<ul>
  <li><a href="https://github.com/trevorstephens/gplearn/tree/ad57cb18caafdb02cca861aea712f1bf3ed5016e">gplearn</a>: symbolic regression/classification</li>
  <li><a href="https://github.com/MilesCranmer/PySR">pysr</a>: fast symbolic regression</li>
  <li><a href="https://github.com/dswah/pyGAM">pygam</a>: generative additive models</li>
  <li><a href="https://github.com/interpretml/interpret">interpretml</a>: boosting-based gam</li>
  <li><a href="https://github.com/h2oai/h2o-3">h20 ai</a>: gams + glms (and more)</li>
  <li><a href="https://github.com/guillermo-navas-palencia/optbinning">optbinning</a>: data discretization / scoring models</li>	
</ul>
</details>

<details>
<summary>Updates</summary>
<ul>
  <li>For updates, star the repo, <a href="https://github.com/csinva/csinva.github.io">see this related repo</a>, or follow <a href="https://twitter.com/csinva_">@csinva_</a></li>
  <li>Please make sure to give authors of original methods / base implementations appropriate credit!</li>
  <li>Contributing: pull requests <a href="https://github.com/csinva/imodels/blob/master/docs/contributing.md">very welcome</a>!</li>
</ul>
</details>


If it's useful for you, please star/cite the package, and make sure to give authors of original methods / base implementations credit:

```r
@software{
    imodels2021,
    title        = {{imodels: a python package for fitting interpretable models}},
    journal      = {Journal of Open Source Software}
    publisher    = {The Open Journal},
    year         = {2021},
    author       = {Singh, Chandan and Nasseri, Keyan and Tan, Yan Shuo and Tang, Tiffany and Yu, Bin},
    volume       = {6},
    number       = {61},
    pages        = {3192},
    doi          = {10.21105/joss.03192},
    url          = {https://doi.org/10.21105/joss.03192},
}

```
Contributions are very welcome!

We have a queue of open things we are working on [here](https://github.com/csinva/imodels/projects/1). Feel free to open an issue or contact @csinva (cs1@berkeley.edu or https://www.linkedin.com/in/csinva/) if you want to contribute!

Before contributing, it would be good to read the sklearn estimator [contributing guide](https://scikit-learn.org/stable/developers/develop.html) and generally be familiar with sklearn. Also please ensure that any added implementations use the appropriate classes from `imodels.util` (e.g. the `Rule` class). 

[Docs](https://csinva.io/imodels/docs/) are built using [pdoc](https://pdoc3.github.io/pdoc/). Build them by changing to the `docs` directory and then running `build_docs.sh`.

[Tests](tests) are run with [pytest](https://docs.pytest.org/en/stable/) (e.g. run `pytest` in the repo directory) - make sure they pass before pushing code, and that new models pass a reasonable set of tests. Note that you might need to install some additional dependencies in order to get the tests to pass.

Project is also on [pypi](https://pypi.org/project/imodels/). Packaged following [this tutorial](https://realpython.com/pypi-publish-python-package/). Relevant commands:
```bash
pip install twine

python setup.py sdist bdist_wheel

twine check dist/*

twine upload dist/*
```
In case you run into issues with installation, here are some things that could help:

If you don't have permissions to install on your machine, use the --user flag:

`pip install git+https://github.com/csinva/imodels --user`

Note that some models (e.g. the ones below) require extra dependencies:

```python
extra_deps = [
    'cvxpy',  # optionally requires cvxpy for slim
    'corels',  # optinally requires corels for optimalrulelistclassifier
    'gosdt',  # optionally requires gosdt for optimaltreeclassifier
    'irf',  # optionally require irf for iterativeRandomForestClassifier
]
```


To test if everything is successfully installed, just try importing imodels from python:

```
python
> import imodels
```

To run example notebooks and/or develop locally, clone the repo then run

`pip install -e ".[dev]"`
---
title: 'imodels: a python package for fitting interpretable models'
tags:
  - python
  - machine learning
  - interpretability
  - explainability
  - transparency
  - decision rules
authors:
  - name: Chandan Singh^[Equal contribution]
    orcid: 0000-0003-0318-2340
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Keyan Nasseri^[Equal contribution]
    affiliation: 1
  - name: Yan Shuo Tan
    affiliation: 2
  - name: Tiffany Tang
    affiliation: 2
  - name: Bin Yu
    affiliation: "1, 2"
affiliations:
 - name: EECS Department, University of California, Berkeley
   index: 1
 - name: Statistics Department, University of California, Berkeley
   index: 2
date: 27 January 2021
bibliography: docs/paper/references.bib
---

# Summary

`imodels` is a Python package for concise, transparent, and accurate predictive modeling.
It provides users a simple interface for fitting and using state-of-the-art interpretable models, all compatible with scikit-learn [@pedregosa2011scikit].
These models can often replace black-box models while improving interpretability and computational efficiency, all without sacrificing predictive accuracy.
In addition, the package provides a framework for developing custom tools and rule-based models for interpretability.

# Statement of need

Recent advancements in machine learning have led to increasingly complex predictive models, often at the cost of interpretability.
There is often a need for models which are inherently interpretable [@rudin2019stop; @murdoch2019definitions], particularly in high-stakes applications such as medicine, biology, and political science.
In these cases, interpretability can ensure that models behave reasonably, identify when models will make errors, and make the models more trusted by domain experts.
Moreover, interpretable models tend to be much more computationally efficient then larger black-box models.

Despite the development of many methods for fitting interpretable models [@molnar2020interpretable], implementations for such models are often difficult to find, use, and compare to one another.
`imodels` aims to fill this gap by providing a simple unified interface and implementation for many state-of-the-art interpretable modeling techniques.

# Features

Interpretable models can take various forms.
\autoref{fig:models} shows four possible forms a model in the `imodels` package can take.
Each form constrains the final model in order to make it interpretable, but there are different methods for fitting the model which differ in their biases and computational costs.
The `imodels` package contains implementations of various such methods and also useful functions for recombining and extending them.

Rule sets consist of a set of rules which each act independently.
    There are different strategies for deriving a rule set, such as Skope-rules [@skope] or Rulefit [@friedman2008predictive].
Rule lists are composed of a set of rules which act in sequence, and include models such as Bayesian rule lists [@letham2015interpretable] or the oneR algorithm [@holte1993very].
Rule trees are similar to rule lists, but allow branching after rules. This includes models such as CART decision trees [@breiman1984classification].
Algebraic models take a final form of simple algebraic expressions, such as supersparse linear integer models [@ustun2016supersparse].

![Examples of different supported model forms. The bottom of each box shows predictions of the corresponding model as a function of $X_1$ and $X_2$.\label{fig:models}](./docs/img/model_table.png){ width=100% }

# Acknowledgements

The code here heavily derives from the wonderful work of previous projects.
In particular, we build upon the following repos and users: [sklearn-expertsys](https://github.com/tmadl/sklearn-expertsys) - by [Tamas Madl](https://github.com/tmadl) and [Benedict](https://github.com/kenben) based on original code by [Ben Letham](http://lethalletham.com/).
We also based many rule-based models on [skope-rules](https://github.com/scikit-learn-contrib/skope-rules) by the [skope-rules team](https://github.com/scikit-learn-contrib/skope-rules/blob/master/AUTHORS.rst) (including [
Nicolas Goix](https://github.com/ngoix), [Florian Gardin](https://github.com/floriangardin), [Jean-Matthieu Schertzer](https://github.com/datajms), Bibi Ndiaye, and Ronan Gautier). 
We also build upon the [rulefit](https://github.com/christophM/rulefit) repository by [Christoph Molnar](https://github.com/christophM).

# References