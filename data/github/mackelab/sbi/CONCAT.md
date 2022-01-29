# v0.18.0 (next release)

## Breaking changes
- Posteriors saved under `sbi` `v0.17.2` or older can not be loaded under `sbi` 
`v0.18.0` or newer.
- `sample_with` can no longer be passed to `.sample()`. Instead, the user has to rerun
`.build_posterior(sample_with=...)`. (#573)
- the `posterior` no longer has the the method `.sample_conditional()`. Using this 
  feature now requires using the `sampler interface` (see tutorial
  [here](https://www.mackelab.org/sbi/tutorial/07_conditional_distributions/)) (#573)

## Major changes
- bugfix for SNPE-C with mixture density networks (#573)
- new `sampler interface` (#573):
```python
from sbi.inference import SNLE, likelihood_estimator_based_potential

inference = SNLE()  # no more prior needed
likelihood_estimator = inference.append_simulations(theta, x).train()

potential_fn, theta_transform = likelihood_estimator_based_potential(likelihood_estimator, prior, x_o)
posterior = MCMCPosterior(potential_fn, proposal=prior, theta_transform=theta_transform)

samples = posterior.sample((100,))
```

## Minor changes
- pairplot takes `ax` and `fig` (#557)
- bugfix for rejection sampling (#561)


# v0.17.2

## Minor changes
- bug fix for transforms in KDE (#552)

# v0.17.1

## Minor changes
- improve kwarg handling for rejection abc and smcabc
- typo and link fixes (#549, thanks to @pitmonticone)
- tutorial notebook on crafting summary statistics with sbi (#511, thanks to @ybernaerts)
- small fixes and improved documenentation for device handling (#544, thanks to @milagorecki)

# v0.17.0

## Major changes
- New API for specifying sampling methods (#487). Old syntax:
```python
posterior = inference.build_posterior(sample_with_mcmc=True)
```
New syntax:
```python
posterior = inference.build_posterior(sample_with="mcmc")  # or "rejection"
```
- Rejection sampling for likelihood(-ratio)-based posteriors (#487)
- MCMC in unconstrained and z-scored space (#510)
- Prior is now allowed to lie on GPU. The prior has to be on the same device as the one
  passed for training (#519).
- Rejection-ABC and SMC-ABC now return the accepted particles / parameters by default,
  or a KDE fit on those particles (`kde=True`) (#525).
- Fast analytical sampling, evaluation and conditioning for `DirectPosterior` trained
  with MDNs (thanks @jnsbck #458). 

## Minor changes
- `scatter` allowed for diagonal entries in pairplot (#510)
- Changes to default hyperparameters for `SNPE_A` (thanks @famura, #496, #497)
- bugfix for `within_prior` checks (#506)


# v0.16.0

## Major changes
- Implementation of SNPE-A (thanks @famura and @theogruner, #474, #478, #480, #482)
- Option to do inference over iid observations with SNLE and SNRE (#484, #488)

## Minor changes
- Fixed unused argument `num_bins` when using `nsf` as density estimator (#465)
- Fixes to adapt to the new support handling in `torch` `v1.8.0` (#469)
- More scalars for monitoring training progress (thanks @psteinb #471)
- Fixed bug in `minimal.py` (thanks @psteinb, #485)
- Depend on `pyknos` `v0.14.2`


# v0.15.1

- add option to pass `torch.data.DataLoader` kwargs to all inference methods (thanks @narendramukherjee, #445)
- fix bug due to release of `torch` `v1.8.0` (#451)
- expose `leakage_correction` parameters for `log_prob` correction in unnormalized 
  posteriors (thanks @famura, #454)


# v0.15.0

## Major changes
- Active subspaces for sensitivity analysis (#394, [tutorial](https://www.mackelab.org/sbi/tutorial/09_sensitivity_analysis/))
- Method to compute the maximum-a-posteriori estimate from the posterior (#412)

## API changes
- `pairplot()`, `conditional_pairplot()`, and `conditional_corrcoeff()` should now be imported from `sbi.analysis` instead of `sbi.utils` (#394).
- Changed `fig_size` to `figsize` in pairplot (#394).
- moved `user_input_checks` to `sbi.utils` (#430).

## Minor changes
- Depend on new `joblib=1.0.0` and fix progress bar updates for multiprocessing (#421).
- Fix for embedding nets with `SNRE` (thanks @adittmann, #425).
- Is it now optional to pass a prior distribution when using SNPE (#426).
- Support loading of posteriors saved after `sbi v0.15.0` (#427, thanks @psteinb).
- Neural network training can be resumed (#431).
- Allow using NSF to estimate 1D distributions (#438).
- Fix type checks in input checks (thanks @psteinb, #439).
- Bugfix for GPU training with SNRE_A (thanks @glouppe, #442).


# v0.14.3

- Fixup for conditional correlation matrix (thanks @JBeckUniTb, #404)
- z-score data using only the training data (#411)


# v0.14.2

- Small fix for SMC-ABC with semi-automatic summary statistics (#402)


# v0.14.1

- Support for training and sampling on GPU including fixes from `nflows` (#331)
- Bug fix for SNPE with neural spline flow and MCMC (#398)
- Small fix for SMC-ABC particles covariance
- Small fix for rejection-classifier (#396)


# v0.14.0

- New flexible interface API (#378). This is going to be a breaking change for users of 
the flexible interface and you will have to change your code. Old syntax:
```python
from sbi.inference import SNPE, prepare_for_sbi

simulator, prior = prepare_for_sbi(simulator, prior)
inference = SNPE(simulator, prior)

# Simulate, train, and build posterior.
posterior = inference(num_simulation=1000)
```
New syntax:
```python
from sbi.inference import SNPE, prepare_for_sbi, simulate_for_sbi

simulator, prior = prepare_for_sbi(simulator, prior)
inference = SNPE(prior)

theta, x = simulate_for_sbi(simulator, proposal=prior, num_simulations=1000)
density_estimator = inference.append_simulations(theta, x).train()
posterior = inference.build_posterior(density_estimator)  # MCMC kwargs go here.
```
More information can be found here [here](https://www.mackelab.org/sbi/tutorial/02_flexible_interface/).
- Fixed typo in docs for `infer` (thanks @glouppe, #370)
- New `RestrictionEstimator` to learn regions of bad simulation outputs (#390)
- Improvements for and new ABC methods (#395)
    - Linear regression adjustment as in Beaumont et al. 2002 for both MCABC and SMCABC
    - Semi-automatic summary statistics as in Fearnhead & Prangle 2012 for both MCABC and SMCABC
    - Small fixes to perturbation kernel covariance estimation in SMCABC.


# v0.13.2

- Fix bug in SNRE (#363)
- Fix warnings for multi-D x (#361)
- Small improvements to MCMC, verbosity and continuing of chains (#347, #348)


# v0.13.1

- Make logging of vectorized numpy slice sampler slightly less verbose and address NumPy future warning (#347)
- Allow continuation of MCMC chains (#348)


# v0.13.0

- Conditional distributions and correlations for analysing the posterior (#321)
- Moved rarely used arguments from pairplot into kwargs (#321)
- Sampling from conditional posterior (#327)
- Allow inference with multi-dimensional x when appropriate embedding is passed (#335)
- Fixes a bug with clamp_and_warn not overriding num_atoms for SNRE and the warning message itself (#338)
- Compatibility with Pyro 1.4.0 (#339)
- Speed up posterior rejection sampling by introducing batch size (#340, #343)
- Allow vectorized evaluation of numpy potentials (#341)
- Adds vectorized version of numpy slice sampler which allows parallel log prob evaluations across all chains (#344)


# v0.12.2

- Bug fix for zero simulations in later rounds (#318)
- Bug fix for sbi.utils.sbiutils.Standardize; mean and std are now registered in state dict (thanks @plcrodrigues, #325)
- Tutorials on embedding_net and presimulated data (thanks @plcrodrigues, #314, #318)
- FAQ entry for pickling error


# v0.12.1

- Bug fix for broken NSF (#310, thanks @tvwenger).


# v0.12.0

- Add FAQ (#293)
- Fix bug in embedding_net when output dimension does not equal input dimension (#299)
- Expose arguments of functions used to build custom networks (#299)
- Implement non-atomic APT (#301)
- Depend on pyknos 0.12 and nflows 0.12
- Improve documentation (#302, #305, thanks to @agramfort)
- Fix bug for 1D uniform priors (#307).

# v0.11.2

- Fixed pickling of SNRE by moving StandardizeInputs (#291)
- Added check to ensure correct round number when presimulated data is provided
- Subclassed Posterior depending on inference algorithm (#282, #285)
- Pinned pyro to v1.3.1 as a temporary workaround (see #288) 
- Detaching weights for MCMC SIR init immediately to save memory (#292)


# v0.11.1

- Bug fix for log_prob() in SNRE (#280)


# v0.11.0

- Changed the API to do multi-round inference (#273)
- Allow to continue inference (#273)


# v0.10.2

- Added missing type imports (#275)
- Made compatible for Python 3.6 (#275)


# v0.10.1

- Added `mcmc_parameters` to init methods of inference methods (#270)
- Fixed detaching of `log_weights` when using `sir` MCMC init (#270)
- Fixed logging for SMC-ABC


# v0.10.0

- Added option to pass external data (#264)
- Added setters for MCMC parameters (#267)
- Added check for `density_estimator` argument (#263)
- Fixed `NeuralPosterior` pickling error (#265)
- Added code coverage reporting (#269)


# v0.9.0

- Added ABC methods (#250)
- Added multiple chains for MCMC and new init strategy (#247)
- Added options for z-scoring for all inference methods (#256)
- Simplified swapping out neural networks (#256)
- Improved tutorials
- Fixed device keyword argument (#253)
- Removed need for passing x-shapes (#259)


# v0.8.0

- First public version
[![PyPI version](https://badge.fury.io/py/sbi.svg)](https://badge.fury.io/py/sbi)
[![Contributions welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat)](https://github.com/mackelab/sbi/blob/master/CONTRIBUTING.md)
[![Tests](https://github.com/mackelab/sbi/workflows/Tests/badge.svg?branch=main)](https://github.com/mackelab/sbi/actions)
[![codecov](https://codecov.io/gh/mackelab/sbi/branch/main/graph/badge.svg)](https://codecov.io/gh/mackelab/sbi)
[![GitHub license](https://img.shields.io/github/license/mackelab/sbi)](https://github.com/mackelab/sbi/blob/master/LICENSE.txt)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02505/status.svg)](https://doi.org/10.21105/joss.02505)

## sbi: simulation-based inference
[Getting Started](https://www.mackelab.org/sbi/tutorial/00_getting_started/) | [Documentation](https://www.mackelab.org/sbi/)

`sbi` is a PyTorch package for simulation-based inference. Simulation-based inference is  
the process of finding parameters of a simulator from observations.

`sbi` takes a Bayesian approach and returns a full posterior distribution
over the parameters, conditional on the observations. This posterior can be amortized (i.e.
useful for any observation) or focused (i.e. tailored to a particular observation), with different
computational trade-offs.

`sbi` offers a simple interface for one-line posterior inference.

```python
from sbi.inference import infer
# import your simulator, define your prior over the parameters
parameter_posterior = infer(simulator, prior, method='SNPE', num_simulations=100)
```
See below for the available methods of inference, `SNPE`, `SNRE` and `SNLE`.


## Installation

`sbi` requires Python 3.6 or higher. We recommend to use a [`conda`](https://docs.conda.io/en/latest/miniconda.html) virtual
environment ([Miniconda installation instructions](https://docs.conda.io/en/latest/miniconda.html])). If `conda` is installed on the system, an environment for
installing `sbi` can be created as follows:
```commandline
# Create an environment for sbi (indicate Python 3.6 or higher); activate it
$ conda create -n sbi_env python=3.7 && conda activate sbi_env
```

Independent of whether you are using `conda` or not, `sbi` can be installed using `pip`:
```commandline
$ pip install sbi
```

To test the installation, drop into a python prompt and run
```python
from sbi.examples.minimal import simple
posterior = simple()
print(posterior)
```

## Inference Algorithms

The following algorithms are currently available:

#### Sequential Neural Posterior Estimation (SNPE)

* [`SNPE_C`](https://www.mackelab.org/sbi/reference/#sbi.inference.snpe.snpe_c.SNPE_C) or `APT` from Greenberg D, Nonnenmacher M, and Macke J [_Automatic
  Posterior Transformation for likelihood-free
  inference_](https://arxiv.org/abs/1905.07488) (ICML 2019).


#### Sequential Neural Likelihood Estimation (SNLE)
* [`SNLE_A`](https://www.mackelab.org/sbi/reference/#sbi.inference.snle.snle_a.SNLE_A) or just `SNL` from Papamakarios G, Sterrat DC and Murray I [_Sequential
  Neural Likelihood_](https://arxiv.org/abs/1805.07226) (AISTATS 2019).


#### Sequential Neural Ratio Estimation (SNRE)

* [`SNRE_A`](https://www.mackelab.org/sbi/reference/#sbi.inference.snre.snre_a.SNRE_A) or `AALR` from Hermans J, Begy V, and Louppe G. [_Likelihood-free Inference with Amortized Approximate Likelihood Ratios_](https://arxiv.org/abs/1903.04057) (ICML 2020).

* [`SNRE_B`](https://www.mackelab.org/sbi/reference/#sbi.inference.snre.snre_b.SNRE_B) or `SRE` from Durkan C, Murray I, and Papamakarios G. [_On Contrastive Learning for Likelihood-free Inference_](https://arxiv.org/abs/2002.03712) (ICML 2020).


## Feedback and Contributions

We would like to hear how `sbi` is working for your inference problems as well as receive bug reports, pull requests and other feedback (see
[contribute](http://www.mackelab.org/sbi/contribute/)).


## Acknowledgements

`sbi` is the successor (using PyTorch) of the
[`delfi`](https://github.com/mackelab/delfi) package. It was started as a fork of Conor
M. Durkan's `lfi`. `sbi` runs as a community project; development is coordinated at the
[mackelab](https://uni-tuebingen.de/en/research/core-research/cluster-of-excellence-machine-learning/research/research/cluster-research-groups/professorships/machine-learning-in-science/). See also [credits](https://github.com/mackelab/sbi/blob/master/docs/docs/credits.md).


## Support

`sbi` has been supported by the German Federal Ministry of Education and Research (BMBF) through the project ADIMEM, FKZ 01IS18052 A-D). [ADIMEM](https://fit.uni-tuebingen.de/Project/Details?id=9199) is a collaborative project between the groups of Jakob Macke (Uni Tübingen), Philipp Berens (Uni Tübingen), Philipp Hennig (Uni Tübingen) and Marcel Oberlaender (caesar Bonn) which aims to develop inference methods for mechanistic models.


## License

[Affero General Public License v3 (AGPLv3)](https://www.gnu.org/licenses/)


## Citation
If you use `sbi` consider citing the [sbi software paper](https://doi.org/10.21105/joss.02505), in addition to the original research articles describing the specifc sbi-algorithm(s) you are using: 

```
@article{tejero-cantero2020sbi,
  doi = {10.21105/joss.02505},
  url = {https://doi.org/10.21105/joss.02505},
  year = {2020},
  publisher = {The Open Journal},
  volume = {5},
  number = {52},
  pages = {2505},
  author = {Alvaro Tejero-Cantero and Jan Boelts and Michael Deistler and Jan-Matthis Lueckmann and Conor Durkan and Pedro J. Gonçalves and David S. Greenberg and Jakob H. Macke},
  title = {sbi: A toolkit for simulation-based inference},
  journal = {Journal of Open Source Software}
}
```
## User experiences, bugs, and feature requests

If you are using `sbi` to infer the parameters of a simulator, we would be delighted to
know how it worked for you. If it didn't work according to plan, please open up an issue
and tell us more about your use case: the dimensionality of the input parameters and of
the output, as well as the setup you used to run inference (i.e. number of simulations,
number of rounds,...).

To report bugs and suggest features (including better documentation), please equally
head over to [issues on GitHub](https://github.com/mackelab/sbi/issues).


## Code contributions

In general, we use pull requests to make changes to `sbi`.


### Development environment

Clone [the repo](https://github.com/mackelab/sbi) and install all the dependencies using
the `environment.yml` file to create a conda environment: `conda env create -f
environment.yml`. If you already have an `sbi` environment and want to refresh
dependencies, just run `conda env update -f environment.yml --prune`.

Alternatively, you can install via `setup.py` using `pip install -e ".[dev]"` (the dev
flag installs development and testing dependencies).

### Contributing inference algorithms

`sbi` was developed to be extensible and we welcome implementations of additional
inference algorithms. Your new inference algorithm should be a class in
`sbi/inference/your_type_of_algorithm/your_algorithm.py`. The class should have a
`__call__()` function which runs inference and returns a posterior object. The posterior
object itself should have a `.sample()` function following the signature of
`sbi/inference/NeuralPosterior`, allowing to draw samples from the posterior.
Currently, `SNPE`, `SNLE`, and `SNRE` all share the `NeuralPosterior` class in
`sbi/inference/posterior.py`, but future versions of `sbi` will refactor them into
separate classes.

### Style conventions

For docstrings and comments, we use [Google
Style](http://google.github.io/styleguide/pyguide.html#38-comments-and-docstrings).

Code needs to pass through the following tools, which are installed alongside `sbi`:

**[black](https://github.com/psf/black)**: Automatic code formatting for Python. You can
run black manually from the console using `black .` in the top directory of the
repository, which will format all files.

**[isort](https://github.com/timothycrosley/isort)**: Used to consistently order
imports. You can run isort manually from the console using `isort` in the top
directory.


## Online documentation

Most of [the documentation](http://mackelab.org/sbi) is written in markdown ([basic
markdown guide](https://guides.github.com/features/mastering-markdown/)).

You can directly fix mistakes and suggest clearer formulations in markdown files simply
by initiating a PR on through GitHub. Click on [documentation
file](https://github.com/mackelab/sbi/tree/master/docs/docs) and look for the little pencil at top right.
# Documentation

The documentation is available at: <http://mackelab.org/sbi>


## Building the Documentation

You can build the docs locally by running the following command from this subfolder:
```bash
jupyter nbconvert --to markdown ../tutorials/*.ipynb --output-dir docs/tutorial/
jupyter nbconvert --to markdown ../examples/*.ipynb --output-dir docs/examples/
mkdocs serve
```

The docs can be updated to GitHub using:
```bash
jupyter nbconvert --to markdown ../tutorials/*.ipynb --output-dir docs/tutorial/
jupyter nbconvert --to markdown ../examples/*.ipynb --output-dir docs/examples/
mkdocs gh-deploy
```

## Contributing FAQ

Create a new markdown file named `question_XX.md` in the `docs/faq` folder, where `XX` 
is a running index for the questions. The file should start with the question as title 
(i.e. starting with a `#`) and then have the answer below. Additionally, you need to 
add a link to your question in the markdown file `docs/faq.md`.# `sbi`: simulation-based inference

`sbi`: A Python toolbox for simulation-based inference.

![using sbi](static/infer_demo.gif)

Inference can be run in a single
line of code:

```python
posterior = infer(simulator, prior, method='SNPE', num_simulations=1000)
```

- To learn about the general motivation behind simulation-based inference, and the
  inference methods included in `sbi`, read on below.

- For example applications to canonical problems in neuroscience, browse the recent
  research article [Training deep neural density estimators to identify mechanistic models of neural dynamics](https://doi.org/10.7554/eLife.56261).

- If you want to get started using `sbi` on your own problem, jump to
  [installation](install.md) and then check out the [tutorial](tutorial/00_getting_started.md).

## Motivation and approach

Many areas of science and engineering make extensive use of complex, stochastic,
numerical simulations to describe the structure and dynamics of the processes being
investigated.

A key challenge in simulation-based science is constraining these simulation models'
parameters, which are intepretable quantities, with observational data. Bayesian
inference provides a general and powerful framework to invert the simulators, i.e.
describe the parameters which are consistent both with empirical data and prior
knowledge.

In the case of simulators, a key quantity required for statistical inference, the
likelihood of observed data given parameters, $\mathcal{L}(\theta) = p(x_o|\theta)$, is
typically intractable, rendering conventional statistical approaches inapplicable.

`sbi` implements three powerful machine-learning methods that address this problem:

- Sequential Neural Posterior Estimation (SNPE),
- Sequential Neural Likelihood Estimation (SNLE), and
- Sequential Neural Ratio Estimation (SNRE).

Depending on the characteristics of the problem, e.g. the dimensionalities of the
parameter space and the observation space, one of the methods will be more suitable.

![](./static/goal.png)

**Goal: Algorithmically identify mechanistic models which are consistent with data.**

Each of the methods above needs three inputs: A candidate mechanistic model, prior
knowledge or constraints on model parameters, and observational data (or summary statistics
thereof).

The methods then proceed by

1. sampling parameters from the prior followed by simulating synthetic data from
   these parameters,
2. learning the (probabilistic) association between data (or
   data features) and underlying parameters, i.e. to learn statistical inference from
   simulated data. The way in which this association is learned differs between the
   above methods, but all use deep neural networks.
3. This learned neural network is then applied to empirical data to derive the full
   space of parameters consistent with the data and the prior, i.e. the posterior
   distribution. High posterior probability is assigned to parameters which are
   consistent with both the data and the prior, low probability to inconsistent
   parameters. While SNPE directly learns the posterior distribution, SNLE and SNRE need
   an extra MCMC sampling step to construct a posterior.
4. If needed, an initial estimate of the posterior can be used to adaptively generate
   additional informative simulations.


## Publications

See [Cranmer, Brehmer, Louppe (2020)](https://doi.org/10.1073/pnas.1912789117) for a recent
review on simulation-based inference.

The following papers offer additional details on the inference methods included in
`sbi`:


### SNPE

- **Fast ε-free Inference of Simulation Models with Bayesian Conditional Density Estimation**<br> by Papamakarios & Murray (NeurIPS 2016) <br>[[PDF]](https://papers.nips.cc/paper/6084-fast-free-inference-of-simulation-models-with-bayesian-conditional-density-estimation.pdf) [[BibTeX]](https://papers.nips.cc/paper/6084-fast-free-inference-of-simulation-models-with-bayesian-conditional-density-estimation/bibtex)

- **Flexible statistical inference for mechanistic models of neural dynamics** <br> by Lueckmann, Goncalves, Bassetto, Öcal, Nonnenmacher & Macke (NeurIPS 2017) <br>[[PDF]](https://papers.nips.cc/paper/6728-flexible-statistical-inference-for-mechanistic-models-of-neural-dynamics.pdf) [[BibTeX]](https://papers.nips.cc/paper/6728-flexible-statistical-inference-for-mechanistic-models-of-neural-dynamics/bibtex)

- **Automatic posterior transformation for likelihood-free inference**<br>by Greenberg, Nonnenmacher & Macke (ICML 2019) <br>[[PDF]](http://proceedings.mlr.press/v97/greenberg19a/greenberg19a.pdf) [[BibTeX]](data:text/plain;charset=utf-8,%0A%0A%0A%0A%0A%0A%40InProceedings%7Bpmlr-v97-greenberg19a%2C%0A%20%20title%20%3D%20%09%20%7BAutomatic%20Posterior%20Transformation%20for%20Likelihood-Free%20Inference%7D%2C%0A%20%20author%20%3D%20%09%20%7BGreenberg%2C%20David%20and%20Nonnenmacher%2C%20Marcel%20and%20Macke%2C%20Jakob%7D%2C%0A%20%20booktitle%20%3D%20%09%20%7BProceedings%20of%20the%2036th%20International%20Conference%20on%20Machine%20Learning%7D%2C%0A%20%20pages%20%3D%20%09%20%7B2404--2414%7D%2C%0A%20%20year%20%3D%20%09%20%7B2019%7D%2C%0A%20%20editor%20%3D%20%09%20%7BChaudhuri%2C%20Kamalika%20and%20Salakhutdinov%2C%20Ruslan%7D%2C%0A%20%20volume%20%3D%20%09%20%7B97%7D%2C%0A%20%20series%20%3D%20%09%20%7BProceedings%20of%20Machine%20Learning%20Research%7D%2C%0A%20%20address%20%3D%20%09%20%7BLong%20Beach%2C%20California%2C%20USA%7D%2C%0A%20%20month%20%3D%20%09%20%7B09--15%20Jun%7D%2C%0A%20%20publisher%20%3D%20%09%20%7BPMLR%7D%2C%0A%20%20pdf%20%3D%20%09%20%7Bhttp%3A%2F%2Fproceedings.mlr.press%2Fv97%2Fgreenberg19a%2Fgreenberg19a.pdf%7D%2C%0A%20%20url%20%3D%20%09%20%7Bhttp%3A%2F%2Fproceedings.mlr.press%2Fv97%2Fgreenberg19a.html%7D%2C%0A%20%20abstract%20%3D%20%09%20%7BHow%20can%20one%20perform%20Bayesian%20inference%20on%20stochastic%20simulators%20with%20intractable%20likelihoods%3F%20A%20recent%20approach%20is%20to%20learn%20the%20posterior%20from%20adaptively%20proposed%20simulations%20using%20neural%20network-based%20conditional%20density%20estimators.%20However%2C%20existing%20methods%20are%20limited%20to%20a%20narrow%20range%20of%20proposal%20distributions%20or%20require%20importance%20weighting%20that%20can%20limit%20performance%20in%20practice.%20Here%20we%20present%20automatic%20posterior%20transformation%20(APT)%2C%20a%20new%20sequential%20neural%20posterior%20estimation%20method%20for%20simulation-based%20inference.%20APT%20can%20modify%20the%20posterior%20estimate%20using%20arbitrary%2C%20dynamically%20updated%20proposals%2C%20and%20is%20compatible%20with%20powerful%20flow-based%20density%20estimators.%20It%20is%20more%20flexible%2C%20scalable%20and%20efficient%20than%20previous%20simulation-based%20inference%20techniques.%20APT%20can%20operate%20directly%20on%20high-dimensional%20time%20series%20and%20image%20data%2C%20opening%20up%20new%20applications%20for%20likelihood-free%20inference.%7D%0A%7D%0A)

### SNLE

- **Sequential neural likelihood: Fast likelihood-free inference with autoregressive flows**<br>by Papamakarios, Sterratt & Murray (AISTATS 2019) <br>[[PDF]](http://proceedings.mlr.press/v89/papamakarios19a/papamakarios19a.pdf) [[BibTeX]](https://gpapamak.github.io/bibtex/snl.bib)

### SNRE

- **Likelihood-free MCMC with Amortized Approximate Likelihood Ratios**<br>by Hermans, Begy & Louppe (ICML 2020) <br>[[PDF]](http://proceedings.mlr.press/v119/hermans20a/hermans20a.pdf)

- **On Contrastive Learning for Likelihood-free Inference**<br>Durkan, Murray & Papamakarios (ICML 2020) <br>[[PDF]](http://proceedings.mlr.press/v119/durkan20a/durkan20a.pdf).
# Installation

`sbi` requires Python 3.6 or higher. We recommend to use a [`conda`](https://docs.conda.io/en/latest/) virtual
environment ([Miniconda installation instructions](https://docs.conda.io/en/latest/miniconda.html)). If `conda` is installed on the system, an environment for
installing `sbi` can be created as follows:
```commandline
# Create an environment for sbi (indicate Python 3.6 or higher); activate it
$ conda create -n sbi_env python=3.7 && conda activate sbi_env
```

Independent of whether you are using `conda` or not, `sbi` can be installed using `pip`:
```commandline
$ pip install sbi
```

To test the installation, drop into a python prompt and run
```python
from sbi.examples.minimal import simple
posterior = simple()
print(posterior)
```
# Credits

## License

`sbi` is licensed under the [Affero General Public License version 3 (AGPLv3)](https://www.gnu.org/licenses/agpl-3.0.html) and

> Copyright (C) 2020 Álvaro Tejero-Cantero, Jakob H. Macke, Jan-Matthis Lückmann,
> Michael Deistler, Jan F. Bölts.

> Copyright (C) 2020 Conor M. Durkan.

## Support

`sbi` has been supported by the German Federal Ministry of Education and Research (BMBF) through the project ADIMEM, FKZ 01IS18052 A-D). [ADIMEM](https://fit.uni-tuebingen.de/Project/Details?id=9199) is a collaborative project between the groups of Jakob Macke (Uni Tübingen), Philipp Berens (Uni Tübingen), Philipp Hennig (Uni Tübingen) and Marcel Oberlaender (caesar Bonn) which aims to develop inference methods for mechanistic models.

![](static/logo_bmbf.svg)


## Important dependencies and prior art

* `sbi` is the successor to [`delfi`](https://github.com/mackelab/delfi), a Theano-based
  toolbox for sequential neural posterior estimation developed at [mackelab](https://uni-tuebingen.de/en/research/core-research/cluster-of-excellence-machine-learning/research/research/cluster-research-groups/professorships/machine-learning-in-science/). If you were
  using `delfi`, we strongly recommend to move your inference over to `sbi`. Please open
  issues if you find unexpected behaviour or missing features. We will consider these
  bugs and give them priority.

* `sbi` as a PyTorch-based toolbox started as a fork of
  [conormdurkan/lfi](https://github.com/conormdurkan/lfi), by [Conor
  M.Durkan](https://conormdurkan.github.io/).

* `sbi` uses density estimators from
[bayesiains/nflows](https://github.com/bayesiains/nsf) by [Conor
M.Durkan](https://conormdurkan.github.io/), [George
Papamakarios](https://gpapamak.github.io/) and [Artur
Bekasov](https://arturbekasov.github.io/). These are proxied through
[`pyknos`](https://github.com/mackelab/pyknos), a package focused on density estimation.

* `sbi` uses `PyTorch` and tries to align with the interfaces (e.g. for probability
  distributions) adopted by `PyTorch`.

* See [README.md](https://github.com/mackelab/sbi/blob/master/README.md) for a list of
  publications describing the methods implemented in `sbi`.


## Citation
If you use `sbi` consider citing the [corresponding paper](https://doi.org/10.21105/joss.02505):
```
@article{tejero-cantero2020sbi,
  doi = {10.21105/joss.02505},
  url = {https://doi.org/10.21105/joss.02505},
  year = {2020},
  publisher = {The Open Journal},
  volume = {5},
  number = {52},
  pages = {2505},
  author = {Alvaro Tejero-Cantero and Jan Boelts and Michael Deistler and Jan-Matthis Lueckmann and Conor Durkan and Pedro J. Gonçalves and David S. Greenberg and Jakob H. Macke},
  title = {sbi: A toolkit for simulation-based inference},
  journal = {Journal of Open Source Software}
}
```
# API Reference

## Inference

::: sbi.inference.base.infer
    rendering:
      show_root_heading: true

::: sbi.utils.user_input_checks.prepare_for_sbi
    rendering:
      show_root_heading: true
      
::: sbi.inference.base.simulate_for_sbi
    rendering:
      show_root_heading: true

::: sbi.inference.snpe.snpe_a.SNPE_A
    rendering:
      show_root_heading: true
    selection:
      filters: [ "!^_", "^__", "!^__class__" ]
      inherited_members: true

::: sbi.inference.snpe.snpe_c.SNPE_C
    rendering:
      show_root_heading: true
    selection:
      filters: [ "!^_", "^__", "!^__class__" ]
      inherited_members: true

::: sbi.inference.snle.snle_a.SNLE_A
    rendering:
      show_root_heading: true
    selection:
      filters: [ "!^_", "^__", "!^__class__" ]
      inherited_members: true

::: sbi.inference.snre.snre_a.SNRE_A
    rendering:
      show_root_heading: true
    selection:
      filters: [ "!^_", "^__", "!^__class__" ]
      inherited_members: true

::: sbi.inference.snre.snre_b.SNRE_B
    rendering:
      show_root_heading: true
    selection:
      filters: [ "!^_", "^__", "!^__class__" ]
      inherited_members: true

::: sbi.inference.abc.mcabc.MCABC
    rendering:
      show_root_heading: true
    selection:
      filters: [ "!^_", "^__", "!^__class__" ]
      inherited_members: true

::: sbi.inference.abc.smcabc.SMCABC
    rendering:
      show_root_heading: true
    selection:
      filters: [ "!^_", "^__", "!^__class__" ]
      inherited_members: true

## Posteriors

::: sbi.inference.posteriors.direct_posterior.DirectPosterior
    rendering:
      show_root_heading: true
    selection:
      filters: [ "!^_", "^__", "!^__class__" ]
      inherited_members: true
      
::: sbi.inference.posteriors.mcmc_posterior.MCMCPosterior
    rendering:
      show_root_heading: true
    selection:
      filters: [ "!^_", "^__", "!^__class__" ]
      inherited_members: true
      
::: sbi.inference.posteriors.rejection_posterior.RejectionPosterior
    rendering:
      show_root_heading: true
    selection:
      filters: [ "!^_", "^__", "!^__class__" ]
      inherited_members: true

## Models

::: sbi.utils.get_nn_models.posterior_nn
    rendering:
      show_root_heading: true
      show_object_full_path: true

::: sbi.utils.get_nn_models.likelihood_nn
    rendering:
      show_root_heading: true
      show_object_full_path: true

::: sbi.utils.get_nn_models.classifier_nn
    rendering:
      show_root_heading: true
      show_object_full_path: true

## Potentials

::: sbi.inference.potentials.posterior_based_potential.posterior_potential
    rendering:
      show_root_heading: true
    selection:
      filters: [ "!^_", "^__", "!^__class__" ]
      inherited_members: true
      
::: sbi.inference.potentials.likelihood_based_potential.likelihood_potential
    rendering:
      show_root_heading: true
    selection:
      filters: [ "!^_", "^__", "!^__class__" ]
      inherited_members: true
      
::: sbi.inference.potentials.ratio_based_potential.ratio_potential
    rendering:
      show_root_heading: true
    selection:
      filters: [ "!^_", "^__", "!^__class__" ]
      inherited_members: true
  
## Analysis

::: sbi.analysis.plot.pairplot
    rendering:
      show_root_heading: true
      show_object_full_path: true
      
::: sbi.analysis.plot.conditional_pairplot
    rendering:
      show_root_heading: true
      show_object_full_path: true
      
::: sbi.analysis.conditional_density.conditional_corrcoeff
    rendering:
      show_root_heading: true
      show_object_full_path: true
# Frequently asked questions

[Can the algorithms deal with invalid data, e.g. NaN or inf?](faq/question_02.md)

[What should I do when my 'posterior samples are outside of the prior support' in SNPE?](faq/question_01.md)

[When using multiple workers, I get a pickling error. Can I still use multiprocessing?](faq/question_03.md)

[Can I use the GPU for training the density estimator?](faq/question_04.md)

[How should I save and load objects in `sbi`?](faq/question_05.md)

[Can I stop neural network training and resume it later?](faq/question_06.md)
{!CONTRIBUTING.md!}
# Learning summary statistics with a neural net

When doing simulation-based inference, it is very important to use well-chosen summary statistics for describing the data generated by the simulator. Very often, these statistics take into account previous domain knowledge. For instance, in the case of the Hodgkin-Huxley model from this [tutorial](https://www.mackelab.org/sbi/examples/00_HH_simulator/), the summary statistics are defined via the function defined in [here](https://github.com/mackelab/sbi/blob/86d9b07238f5a0176638fecdd5622694d92f2962/examples/HH_helper_functions.py#L159), which takes a 120 ms recording as input (a 12000-dimensional input vector) and outputs a 7-dimensional feature vector containing different statistical descriptors of the recording (e.g., number of spikes, average value, etc.). However, in some occasions, it might be of interest to actually **learn from the data** which summary statistics to use.

`sbi` offers functionality to learn summary statistics from (potentially high-dimensional) simulation outputs with a neural network. In `sbi`, this neural network is referred to as `embedding_net`. If an `embedding_net` is specified, the simulation outputs are passed through the `embedding_net`, whose outputs are then passed to the neural density estimator. The parameters of the `embedding_net` are updated together with the parameters of the neural density estimator.

NB: only `SNPE` and `SNRE` methods can use an `embedding_net` to learn summary statistics from simulation outputs. `SNLE` does not offer this functionality since the simulation outputs $x$ are the outputs of the neural density estimator in `SNLE`.

In the example that follows, we illustrate a situation where the data points generated by the simulator model are high-dimensional (32 by 32 images) and we use a convolutional neural network as summary statistics extractor.

Note, you find the original version of this notebook at [https://github.com/mackelab/sbi/blob/main/tutorials/05_embedding_net.ipynb](https://github.com/mackelab/sbi/blob/main/tutorials/05_embedding_net.ipynb) in the `sbi` repository.

First of all, we import all the packages required for running the tutorial


```python
import matplotlib.pyplot as plt
import numpy as np
import torch
import torch.nn as nn 
import torch.nn.functional as F 
from sbi import utils
from sbi import inference
import numpy as np

# set seed for numpy and torch
seed = 42
np.random.seed(seed)
torch.manual_seed(seed)
```

## The simulator model

The simulator model that we consider has two parameters: $r$ and $\theta$. On each run, it generates 100 two-dimensional points centered around $(r \cos(\theta), r \sin(\theta))$ and perturbed by a Gaussian noise with variance 0.01. Instead of simply outputting the $(x,y)$ coordinates of each data point, the model generates a grayscale image of the scattered points with dimensions 32 by 32. This image is further perturbed by an uniform noise with values betweeen 0 and 0.2. The code below defines such model.


```python
 def simulator_model(parameter, return_points=False):
    """ Simulator model with two-dimensional input parameter and 1024-dimensional output
    
    This simulator serves as a basic example for using a neural net for learning summary features. 
    It has only two input parameters but generates high-dimensional output vectors.
    The data is generated as follows:
        (-) Input:  parameter = [r, theta]
        (1) Generate 100 two-dimensional points centered around (r cos(theta),r sin(theta))  
            and perturbed by a Gaussian noise with variance 0.01
        (2) Create a grayscale image I of the scattered points with dimensions 32 by 32
        (3) Perturb I with an uniform noise with values betweeen 0 and 0.2
        (-) Output: I 
        
    Parameters
    ----------
    parameter : array-like, shape (2)
        The two input parameters of the model, ordered as [r, theta]
    return_points : bool (default: False)
        Whether the simulator should return the coordinates of the simulated data points as well
        
    Returns
    -------
    I: torch tensor, shape (1, 1024)    
        Output flattened image
    (optional) points: array-like, shape (100, 2)
        Coordinates of the 2D simulated data points 
    
    """
    r = parameter[0]
    theta = parameter[1]

    sigma_points = 0.10
    npoints = 100
    points = []
    for _ in range(npoints):
        x = r * np.cos(theta) + sigma_points * np.random.randn()
        y = r * np.sin(theta) + sigma_points * np.random.randn()
        points.append([x, y])
    points = np.array(points)

    nx = 32
    ny = 32
    sigma_image = 0.20
    I = np.zeros((nx, ny))
    for point in points:
        pi = int((point[0] - (-1)) / ((+1) - (-1)) * nx)
        pj = int((point[1] - (-1)) / ((+1) - (-1)) * ny) 
        if (pi < nx) and (pj < ny):   
            I[pi, pj] = 1
    I = I + sigma_image * np.random.rand(nx, ny)    
    I = I.T
    I = I.reshape(1,-1)
    I = torch.tensor(I, dtype=torch.get_default_dtype())

    if return_points:
        return I, points
    else:
        return I
```

The figure below shows an example of the output of the simulator when $r = 0.70$ and $\theta = \pi/4$


```python
# simulate samples
true_parameter = torch.tensor([0.70, np.pi/4])
x_observed, x_points = simulator_model(true_parameter, return_points=True)

# plot the observation
fig, ax = plt.subplots(facecolor='white', figsize=(11.15, 5.61), ncols=2, constrained_layout=True)
circle = plt.Circle((0, 0), 1.0, color='k', ls='--', lw=0.8, fill=False)
ax[0].add_artist(circle)
ax[0].scatter(x_points[:,0], x_points[:,1], s=20)
ax[0].set_xlabel('x')
ax[0].set_ylabel('y')
ax[0].set_xlim(-1, +1)
ax[0].set_xticks([-1, 0.0, +1.0])
ax[0].set_ylim(-1, +1)
ax[0].set_yticks([-1, 0.0, +1.0])
ax[0].set_title(r'original simulated points with $r = 0.70$ and $\theta = \pi/4$')
ax[1].imshow(x_observed.view(32, 32), origin='lower', cmap='gray')
ax[1].set_xticks([]); ax[1].set_yticks([])
ax[1].set_title('noisy observed data (gray image with 32 x 32 pixels)')
```

## Defining an `embedding_net`

An inference procedure applied to the output data from this simulator model determines the posterior distribution of $r$ and $\theta$ given an observation of $x$, which lives in a 1024 dimensional space (32 x 32 = 1024). To avoid working directly on these high-dimensional vectors, one can use a convolutional neural network (CNN) that takes the 32x32 images as input and encodes them into 8-dimensional feature vectors. This CNN is trained along with the neural density estimator of the inference procedure and serves as an automatic summary statistics extractor. 

We define and instantiate the CNN as follows:


```python
class SummaryNet(nn.Module): 
    
    def __init__(self): 
        super().__init__()
        # 2D convolutional layer
        self.conv1 = nn.Conv2d(in_channels=1, out_channels=6, kernel_size=5, padding=2)
        # Maxpool layer that reduces 32x32 image to 4x4
        self.pool = nn.MaxPool2d(kernel_size=8, stride=8)
        # Fully connected layer taking as input the 6 flattened output arrays from the maxpooling layer
        self.fc = nn.Linear(in_features=6*4*4, out_features=8) 
        
    def forward(self, x):
        x = x.view(-1, 1, 32, 32)
        x = self.pool(F.relu(self.conv1(x)))
        x = x.view(-1, 6*4*4)
        x = F.relu(self.fc(x))
        return x

embedding_net = SummaryNet()
```

## The inference procedure

With the `embedding_net` defined and instantiated, we can follow the usual workflow of an inference procedure in `sbi`. The `embedding_net` object appears as an input argument when instantiating the neural density estimator with `utils.posterior_nn`.


```python
# set prior distribution for the parameters 
prior = utils.BoxUniform(low=torch.tensor([0.0, 0.0]), 
                             high=torch.tensor([1.0, 2*np.pi]))                           

# make a SBI-wrapper on the simulator object for compatibility
simulator_wrapper, prior = inference.prepare_for_sbi(simulator_model, prior)

# instantiate the neural density estimator
neural_posterior = utils.posterior_nn(model='maf', 
                                      embedding_net=embedding_net,
                                      hidden_features=10,
                                      num_transforms=2)

# setup the inference procedure with the SNPE-C procedure
inference = inference.SNPE(simulator_wrapper, prior, 
                           density_estimator=neural_posterior, 
                           show_progress_bars=True)

# run the inference procedure on one round and 10000 simulated data points
posterior = inference(num_simulations=10000)
```

## Visualizing the results

We now generate 50000 samples of the posterior distribution of $r$ and $\theta$ when observing an input data point $x$ generated from the `simulator model` with $r = 0.70$ and $\theta = \pi/4$. 


```python
# generate posterior samples
true_parameter = torch.tensor([0.70, np.pi/4])
x_observed = simulator_model(true_parameter)
samples = posterior.set_default_x(x_observed).sample((50000,))
```

The figure below shows the statistics of the generated samples.


```python
# create the figure
fig, ax = utils.pairplot(samples, 
                             points=true_parameter,
                             labels=['r', r'$\theta$'], 
                             limits=[[0, 1], [0, 2*np.pi]],
                             points_colors='r',
                             points_offdiag={'markersize': 6},
                             fig_size=[7.5, 6.4])
```
# What should I do when my 'posterior samples are outside of the prior support' in SNPE?

When working with **multi-round** SNPE, you might have experienced the following
warning: 
```
Only x% posterior samples are within the prior support. It may take a long time to collect the remaining 10000 samples. Consider interrupting (Ctrl-C) and switching to 'sample_with_mcmc=True'.
```

This reason for this issue is described in more detail 
[here](https://arxiv.org/abs/2002.03712) and [here](https://arxiv.org/abs/1905.07488). 
The following fixes are possible:  

- sample with MCMC: `samples = posterior((num_samples,), x=x_o, sample_with_mcmc=True)`.
This will make sampling slower, but samples will not 'leak'.  

- resort to single-round SNPE and (if necessary) increase your simulation budget.  

- if your prior is either Gaussian (torch.distributions.multivariateNormal) or Uniform 
(sbi.utils.BoxUniform), you can avoid leakage by using a mixture density network as 
density estimator. I.e., using the 
[flexible interface](https://www.mackelab.org/sbi/tutorial/03_flexible_interface/), set 
`density_estimator='mdn'`. When running inference, there should be a print statement 
"Using SNPE-C with non-atomic loss"

- use a different algorithm, e.g. SNRE and SNLE. Note, however, that these algorithms
can have different issues and potential pitfalls.  
# When using multiple workers, I get a pickling error. Can I still use multiprocessing?

Yes, but you will have to make a few adjustments to your code. 

Some background: when using `num_workers > 1`, you might experience an error that a 
certain object from your simulator could not be pickled (an example can be found
[here](https://github.com/mackelab/sbi/issues/317)).

This can be fixed by forcing `sbi` to pickle with `dill` instead of the default 
`cloudpickle`. To do so, adjust your code as 
follows:

- Install `dill`:
```
pip install dill
```
- At the very beginning of your python script, set the pickler to `dill`:
```python
from joblib.externals.loky import set_loky_pickler
set_loky_pickler("dill")
```
- Move all imports required by your simulator into the simulator:
```python
# Imports specified outside of the simulator will break dill:
import torch
def my_simulator(parameters):
    return torch.ones(1,10)

# Therefore, move the imports into the simulator:
def my_simulator(parameters):
    import torch
    return torch.ones(1,10)
```

### Alternative: parallelize yourself

You can also write your own code to parallelize simulations with whatever 
multiprocessing framework you prefer. You can then simulate your data outside of `sbi` and pass the simulated data as shown in the 
[flexible interface](https://www.mackelab.org/sbi/tutorial/02_flexible_interface/): 


### Some more background

`sbi` uses `joblib` to parallelize simulations, which in turn uses `pickle` or 
`cloudpickle` to serialize the simulator. Almost all simulators will be picklable with 
`cloudpickle`, but we have experienced issues e.g. with `neuron` simulators, see
[here](https://github.com/mackelab/sbi/issues/317).
# Can I use the GPU for training the density estimator?

TLDR; Yes, by passing `device="cuda"` and by passing a prior that lives on the device
name your passed. But no speed-ups for default density estimators.

Yes. When creating the inference object in the flexible interface, you can pass the
`device` as an argument, e.g.,

```python
inference = SNPE(prior, device="cuda", density_estimator="maf")
```

The device is set to `"cpu"` by default, and it can be set to anything, as long as it
maps to an existing PyTorch CUDA device. `sbi` will take care of copying the `net` and
the training data to and from the `device`. 
Note that the prior must be on the training device already, e.g., when passing `device="cuda:0"`,
make sure to pass a prior object that was created on that device, e.g., 
`prior = torch.distributions.MultivariateNormal(loc=torch.zeros(2, device="cuda:0"), 
                                                covariance_matrix=torch.eye(2, device="cuda:0"))`.

## Performance

Whether or not you reduce your training time when training on a GPU depends on the
problem at hand. We provide a couple of default density estimators for `SNPE`, `SNLE`
and `SNRE`, e.g., a mixture density network (`density_estimator="mdn"`) or a Masked
Autoregressive Flow (`density_estimator="maf"`). For those default density estimators
we do **not** expect a speed up. This is because the underlying neural networks are
quite shallow and not tall, e.g., they do not have many parameters or matrix
operations that profit a lot from being executed on the GPU. 

A speed up through training on the GPU will most likely become visible when you are
using convolutional modules in your neural networks. E.g., when passing an embedding
net for image processing like in this example: [https://github.com/mackelab/sbi/blob/main/tutorials/05_embedding_net.ipynb](https://github.com/mackelab/sbi/blob/main/tutorials/05_embedding_net.ipynb). 

# Can I use a custom prior with sbi?

`sbi` works with torch distributions only so we recommend to use those whenever possible. For example, when you are used to using `scipy.stats` distributions as priors then we recommend using the corresponding `torch.distributions`, most common distributions are implemented there. 

In case you want to use a custom prior that is not in the set of common distributions that's possible as well: You need to write a prior class that mimicks the behaviour of a [`torch.distributions.Distribution`](https://pytorch.org/docs/stable/_modules/torch/distributions/distribution.html#Distribution) class. Then `sbi` will wrap this class to make it a fully functional torch `Distribution`.

Essentially, the class needs two methods:

- `.sample(sample_shape)`, where sample_shape is a shape tuple, e.g., `(n,)`, and returns a batch of n samples, e.g., of shape (n, 2)` for a two dimenional prior.
- `.log_prob(value)` method that returns the "log probs" of parameters under the prior, e.g., for a batches of n parameters with shape `(n, ndims)` it should return a log probs array of shape `(n,)`.

For sbi > 0.17.2 this could look like the following:

```python
class CustomUniformPrior:
    """User defined numpy uniform prior.

    Custom prior with user-defined valid .sample and .log_prob methods.
    """

    def __init__(self, lower: Tensor, upper: Tensor, return_numpy: bool = False):
        self.lower = lower
        self.upper = upper
        self.dist = BoxUniform(lower, upper)
        self.return_numpy = return_numpy

    def sample(self, sample_shape=torch.Size([])):
        samples = self.dist.sample(sample_shape)
        return samples.numpy() if self.return_numpy else samples

    def log_prob(self, values):
        if self.return_numpy:
            values = torch.as_tensor(values)
        log_probs = self.dist.log_prob(values)
        return log_probs.numpy() if self.return_numpy else log_probs
```

Once you have such a class you can wrap into a `Distribution` using the `process_prior` function `sbi` provides:

```python
from sbi.utils import process_prior

custom_prior = CustomUniformPrior(torch.zeros(2), torch.ones(2))
prior, *_ = process_prior(custom_prior)  # Keeping only the first return.
# use this wrapped prior in sbi...
```

In `sbi` it is sometimes required to check the support of the prior, e.g., when the prior support is bounded and one wants to reject samples from the posterior density estimator that lie outside the prior support. In torch `Distributions` this is handled automatically, however, when using a custom prior it is not. Thus, 
if your prior has bounded support (like the one above) it makes sense to pass the bounds to the wrapper function such that `sbi` can pass them to torch `Distributions`:

```python
from sbi.utils import process_prior

custom_prior = CustomUniformPrior(torch.zeros(2), torch.ones(2))
prior = process_prior(custom_prior, 
                      custom_prior_wrapper_kwargs=dict(lower_bound=torch.zeros(2), 
                                                       upper_bound=torch.ones(2)))
# use this wrapped prior in sbi...
```

Note that in `custom_prior_wrapper_kwargs` you can pass additinal arguments for the wrapper, e.g., `validate_args` or `arg_constraints` see the `Distribution` documentation for more details. 

If you are running sbi < 0.17.2 and use `SNLE` the code above will produce a `NotImplementedError` (see [#581](https://github.com/mackelab/sbi/issues/581)). In this case you need to update to a newer version of `sbi` or use `SNPE` instead. 

# Can I stop neural network training and resume it later?

Many clusters have a time limit and `sbi` might exceed this limit. You can circumvent this problem by using the [flexible interface](https://www.mackelab.org/sbi/tutorial/02_flexible_interface/). After simulations are finished, `sbi` trains a neural network. If this process takes too long, you can stop training and resume it later. The syntax is:

```python
inference = SNPE(prior=prior)
inference = inference.append_simulations(theta, x)
inference.train(max_num_epochs=300)  # Pick `max_num_epochs` such that it does not exceed the runtime.

with open("path/to/my/inference.pkl", "wb") as handle:
    dill.dump(inference, handle)

# To resume training:
with open("path/to/my/inference.pkl", "rb") as handle:
    inference_from_disk = dill.load(handle)
inference_from_disk.train(resume_training=True, max_num_epochs=600)  # Run epochs 301 until 600 (or stop early).
posterior = inference_from_disk.build_posterior()
```

Note that the inference object can not be saved with `pickle`. To save it, you will have to install and use [dill](https://pypi.org/project/dill/). Another solution is described [here](https://www.mackelab.org/sbi/faq/question_04/).# Can the algorithms deal with invalid data, e.g. NaN or inf?

Yes. By default, whenever a simulation returns at least one `NaN` or `inf`, it is 
completely excluded from the training data. In other words, the simulation is simply 
discarded.

In cases where a very large fraction of simulations return `NaN` or `inf`, discarding many simulations can be wasteful. There are two options to deal with this: Either, you use the `RestrictionEstimator` to learn regions in parameter space that do not produce `NaN` or `inf`, see [here](https://www.mackelab.org/sbi/tutorial/08_restriction_estimator/). Alternatively, you can manually substitute the 'invalid' values with a reasonable replacement. 
I.e., at the end of your simulation code, you search for invalid entries and replace 
them with a floating point number. Importantly, in order for neural network training 
work well, the floating point number should still be in a reasonable range, i.e. maybe a
 few standard deviations outside of 'good' values.

If you are running **multi-round** SNPE, however, things can go fully wrong if invalid
 data are encountered. In that case, you will get the following warning
```
When invalid simulations are excluded, multi-round SNPE-C can leak into the regions where parameters led to invalid simulations. This can lead to poor results.
```
Hence, if you are running multi-round SNPE and a significant fraction of simulations
returns at least one invalid number, we strongly recommend to manually replace the value
 in your simulation code as described above (or resort to single-round SNPE or use a 
 different method).
# How should I save and load objects in `sbi`?

`NeuralPosterior` objects are picklable.
```python
import pickle

# ... run inference
posterior = inference.build_posterior()

with open("/path/to/my_posterior.pkl", "wb") as handle:
    pickle.dump(posterior, handle)
```

Note: posterior objects that were saved under `sbi v0.17.2` or older can not be loaded under `sbi v0.18.0` or newer.

Note: if you try to load a posterior that was saved under `sbi v0.14.x` or earlier under `sbi v0.15.x` until `sbi v0.17.x`, you have to add:
```python
import sys
from sbi.utils import user_input_checks_utils

sys.modules["sbi.user_input.user_input_checks_utils"] = user_input_checks_utils
```
to your script before loading the posterior.


As of `sbi v0.18.0`, `NeuralInference` objects are also picklable.
```python
import pickle

# ... run inference
posterior = inference.build_posterior()

with open("/path/to/my_inference.pkl", "wb") as handle:
    pickle.dump(inference, handle)
```
However, saving and loading the `inference` object will slightly modify the object (in order to make it serializable). These modifications lead to the following two changes in behaviour:
1) Retraining from scratch is not supported, i.e. `.train(..., retrain_from_scratch=True)` does not work.
2) When the loaded object calls the `.train()` method, it generates a new tensorboard summary writer (instead of appending to the current one).