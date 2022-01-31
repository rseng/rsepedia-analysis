<p align="center">
<img src="images/sysidentpy-logo.svg" width="640" height="320" />
</p>

[![DOI](https://img.shields.io/badge/DOI-10.21105%2Fjoss.02384-%23FF7800)](https://joss.theoj.org/papers/10.21105/joss.02384)
[![PyPI version](https://img.shields.io/pypi/v/sysidentpy?color=%23ff7800)](https://pypi.org/project/sysidentpy/)
[![License](https://img.shields.io/pypi/l/sysidentpy?color=%23FF7800)](https://opensource.org/licenses/BSD-3-Clause)
[![openissues](https://img.shields.io/github/issues/wilsonrljr/sysidentpy?color=%23FF7800)](https://github.com/wilsonrljr/sysidentpy/issues)
[![issuesclosed](https://img.shields.io/github/issues-closed-raw/wilsonrljr/sysidentpy?color=%23FF7800)](https://github.com/wilsonrljr/sysidentpy/issues)
[![downloads](https://img.shields.io/pypi/dm/sysidentpy?color=%23FF7800)](https://pypi.org/project/sysidentpy/)
[![python](https://img.shields.io/pypi/pyversions/sysidentpy?color=%23FF7800)](https://pypi.org/project/sysidentpy/)
[![status](https://img.shields.io/pypi/status/sysidentpy?color=%23FF7800)](https://pypi.org/project/sysidentpy/)
[![discord](https://img.shields.io/discord/711610087700955176?color=%23FF7800&label=discord)](https://discord.gg/7afBSzU4)
[![contributors](https://img.shields.io/github/contributors/wilsonrljr/sysidentpy?color=%23FF7800)](https://github.com/wilsonrljr/sysidentpy/graphs/contributors)
[![forks](https://img.shields.io/github/forks/wilsonrljr/sysidentpy?style=social)](https://github.com/wilsonrljr/sysidentpy/network/members)
[![stars](https://img.shields.io/github/stars/wilsonrljr/sysidentpy?style=social)](https://github.com/wilsonrljr/sysidentpy/stargazers)



SysIdentPy is a Python module for System Identification using **NARMAX** models built on top of **numpy** and is distributed under the 3-Clause BSD license.

# Note
The update **v0.1.7**  has been released with major changes and additional features (Fourier basis function, NAR and NFIR models, possibility to select the lag of the residues for Extended Least Squares algorithm and many more).

There are several API modifications and you will need to change your code to have the new (and upcoming) features.

Check the examples of how to use the new version in the [documentation page](<http://sysidentpy.org/notebooks.html>).
  
For more details, please see the [changelog](<http://sysidentpy.org/changelog/v0.1.7.html>).

# Documentation

- Website: https://sysidentpy.org

# Examples

## SysIdentPy now support NARX Neural Network and General estimators, e.g., sklearn estimators and Catboost.

### Exemples
```python
from torch import nn
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sysidentpy.metrics import mean_squared_error
from sysidentpy.utils.generate_data import get_siso_data


# Generate a dataset of a simulated dynamical system
x_train, x_valid, y_train, y_valid = get_siso_data(n=1000,
                                                   colored_noise=False,
                                                   sigma=0.001,
                                                   train_percentage=80)
```

#### Building Polynomial NARX models with FROLS algorithm

```python
from sysidentpy.model_structure_selection import FROLS
from sysidentpy.basis_function import Polynomial
from sysidentpy.utils.display_results import results
from sysidentpy.utils.plotting import plot_residues_correlation, plot_results
from sysidentpy.residues.residues_correlation import compute_residues_autocorrelation
from sysidentpy.residues.residues_correlation import compute_cross_correlation

basis_function=Polynomial(degree=2)
model = PolynomialNarmax(
  order_selection=True,
  n_info_values=10,
  extended_least_squares=False,
  ylag=2, xlag=2,
  info_criteria='aic',
  estimator='least_squares',
  basis_function=basis_function
)
model.fit(X=x_train, y=y_train)
yhat = model.predict(X=x_valid, y=y_valid)
print(rrse)
r = pd.DataFrame(
	results(
		model.final_model, model.theta, model.err,
		model.n_terms, err_precision=8, dtype='sci'
		),
	columns=['Regressors', 'Parameters', 'ERR'])
print(r)
	
Regressors     Parameters        ERR
0        x1(k-2)     0.9000  0.95556574
1         y(k-1)     0.1999  0.04107943
2  x1(k-1)y(k-1)     0.1000  0.00335113

plot_results(y=y_valid, yhat=yhat, n=1000)
ee = compute_residues_autocorrelation(y_valid, yhat)
plot_residues_correlation(data=ee, title="Residues", ylabel="$e^2$")
x1e = compute_cross_correlation(y_valid, yhat, x2_val)
plot_residues_correlation(data=x1e, title="Residues", ylabel="$x_1e$")
```
![polynomial](./examples/figures/polynomial_narmax.png)

#### NARX Neural Network
```python
from sysidentpy.neural_network import NARXNN

class NARX(nn.Module):
    def __init__(self):
        super().__init__()
        self.lin = nn.Linear(4, 10)
        self.lin2 = nn.Linear(10, 10)
        self.lin3 = nn.Linear(10, 1)
        self.tanh = nn.Tanh()

    def forward(self, xb):
        z = self.lin(xb)
        z = self.tanh(z)
        z = self.lin2(z)
        z = self.tanh(z)
        z = self.lin3(z)
        return z

narx_net = NARXNN(net=NARX(),
                  ylag=2,
                  xlag=2,
                  loss_func='mse_loss',
                  optimizer='Adam',
                  epochs=200,
                  verbose=False,
                  optim_params={'betas': (0.9, 0.999), 'eps': 1e-05} # optional parameters of the optimizer
)

train_dl = narx_net.data_transform(x_train, y_train)
valid_dl = narx_net.data_transform(x_valid, y_valid)
narx_net.fit(train_dl, valid_dl)
yhat = narx_net.predict(x_valid, y_valid)
ee, ex, extras, lam = narx_net.residuals(x_valid, y_valid, yhat)
narx_net.plot_result(y_valid, yhat, ee, ex)
```
![neural](/examples/figures/narx_network.png)

#### Catboost-narx
```python
from sysidentpy.general_estimators import NARX
from catboost import CatBoostRegressor

catboost_narx = NARX(base_estimator=CatBoostRegressor(iterations=300,
                                                      learning_rate=0.1,
                                                      depth=6),
                     xlag=2,
                     ylag=2,
                     fit_params={'verbose': False}
)

catboost_narx.fit(x_train, y_train)
yhat = catboost_narx.predict(x_valid, y_valid)
ee, ex, extras, lam = catboost_narx.residuals(x_valid, y_valid, yhat)
catboost_narx.plot_result(y_valid, yhat, ee, ex)
```
![catboost](/examples/figures/catboost_narx.png)

#### Catboost without NARX configuration

The following is the Catboost performance without the NARX configuration.


```python

def plot_results(yvalid, yhat):
    _, ax = plt.subplots(figsize=(14, 8))
    ax.plot(y_valid[:200], label='Data', marker='o')
    ax.plot(yhat[:200], label='Prediction', marker='*')
    ax.set_xlabel("$n$", fontsize=18)
    ax.set_ylabel("$y[n]$", fontsize=18)
    ax.grid()
    ax.legend(fontsize=18)
    plt.show()

catboost = CatBoostRegressor(iterations=300,
                            learning_rate=0.1,
                            depth=6)
catboost.fit(x_train, y_train, verbose=False)
plot_results(y_valid, catboost.predict(x_valid))
```
![catboost](/examples/figures/catboost.png)

The examples directory has several Jupyter notebooks presenting basic tutorials of how to use the package and some specific applications of sysidentpy. Try it out!

# Requirements

SysIdentPy requires:

- Python (>= 3.6)
- NumPy (>= 1.5.0) for all numerical algorithms
- Matplotlib >= 1.5.2 for static plotting and visualizations
- Pytorch (>=1.7.1) for building feed-forward neural networks

| Platform | Status |
| --------- | -----:|
| Linux | ok |
| Windows | ok |
| macOS | ok |

**SysIdentPy do not to support Python 2.7.**

A few examples require pandas >= 0.18.0. However, it is not required to use sysidentpy.

# Installation

The easiest way to get sysidentpy running is to install it using ``pip``
~~~~~~~~~~~~~~~~~~~~~~
pip install sysidentpy
~~~~~~~~~~~~~~~~~~~~~~

We will make it available at conda repository as soon as possible.

# Changelog

See the [changelog]( <http://sysidentpy.org/changelog/v0.1.6.html>) for a history of notable changes to SysIdentPy.

# Development

We welcome new contributors of all experience levels. The sysidentpy community goals are to be helpful, welcoming, and effective.

*Note*: we use the `pytest` package for testing. The test functions are located in tests subdirectories at each folder inside **SysIdentPy**, which check the validity of the algorithms.

Run the `pytest` in the respective folder to perform all the tests of the corresponding sub-packages.

Currently, we have around 81% of code coverage.

You can install pytest using
~~~~~~~~~~~~~~~~~~~~~~
pip install -U pytest
~~~~~~~~~~~~~~~~~~~~~~

### Example of how to run the tests:

Open a terminal emulator of your choice and go to a subdirectory, e.g,
~~~~~~~~~~~~~~~~~~~~
\sysidentpy\metrics\
~~~~~~~~~~~~~~~~~~~~

Just type `pytest` and you get a result like

~~~~~~~~
========== test session starts ==========

platform linux -- Python 3.7.6, pytest-5.4.2, py-1.8.1, pluggy-0.13.1

rootdir: ~/sysidentpy

plugins: cov-2.8.1

collected 12 items

tests/test_regression.py ............ [100%]

========== 12 passed in 2.45s ==================
~~~~~~~~~~~~~~
You can also see the code coverage using the `pytest-cov` package. First, install `pytest-cov` using
~~~
pip install pytest-cov
~~~
Run the command below in the SysIdentPy root directory, to generate the report.
~~~
pytest --cov=.
~~~

# Important links

- Official source code repo: https://github.com/wilsonrljr/sysidentpy

- Download releases: https://pypi.org/project/sysidentpy/

# Source code

You can check the latest sources with the command::
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
git clone https://github.com/wilsonrljr/sysidentpy.git
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Project History

The project was started by Wilson R. L. Junior, Luan Pascoal and Samir A. M. Martins as a project for System Identification discipline. Samuel joined early in 2019.

The project is actively maintained by Wilson R. L. Junior and looking for contributors.

# Communication

- Discord server: https://discord.gg/8eGE3PQ

  [![discord](https://img.shields.io/discord/711610087700955176?color=%23FF7800&label=discord)](https://discord.gg/7afBSzU4)


- Website: http://sysidentpy.org

# Citation
[![DOI](https://img.shields.io/badge/DOI-10.21105%2Fjoss.02384-%23FF7800)](https://joss.theoj.org/papers/10.21105/joss.02384)

If you use SysIdentPy on your project, please [drop me a line](mailto:wilsonrljr@outlook.com).

If you use SysIdentPy on your scientific publication, we would appreciate citations to the following paper:

- Lacerda et al., (2020). SysIdentPy: A Python package for System Identification using NARMAX models. Journal of Open Source Software, 5(54), 2384, https://doi.org/10.21105/joss.02384

```
@article{Lacerda2020,
  doi = {10.21105/joss.02384},
  url = {https://doi.org/10.21105/joss.02384},
  year = {2020},
  publisher = {The Open Journal},
  volume = {5},
  number = {54},
  pages = {2384},
  author = {Wilson Rocha Lacerda Junior and Luan Pascoal Costa da Andrade and Samuel Carlos Pessoa Oliveira and Samir Angelo Milani Martins},
  title = {SysIdentPy: A Python package for System Identification using NARMAX models},
  journal = {Journal of Open Source Software}
}
```

# Inspiration

The documentation and structure (even this section) is openly inspired by sklearn, einsteinpy, and many others as we used (and keep using) them to learn.
# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, sex characteristics, gender identity and expression,
level of experience, education, socio-economic status, nationality, personal
appearance, race, religion, or sexual identity and orientation.

## Our Standards

Examples of behavior that contributes to creating a positive environment
include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or
 advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic
 address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a
 professional setting

## Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

## Scope

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at wilsonrljr@outlook.com. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see
https://www.contributor-covenant.org/faq
---
title: 'SysIdentPy: A Python package for System Identification using NARMAX models'
tags:
  - Python
  - System Identification
  - NARMAX
  - Dynamical Systems
  - Model Structure Selection
authors:
  - name: Wilson Rocha Lacerda Junior
    orcid: 0000-0002-3263-1152
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: Luan Pascoal da Costa Andrade
    affiliation: 1
  - name: Samuel Carlos Pessoa Oliveira
    affiliation: 1
  - name: Samir Angelo Milani Martins
    affiliation: "1, 2"
affiliations:
 - name: GCoM - Modeling and Control Group at Federal University of São João del-Rei, Brazil
   index: 1
 - name: Department of Electrical Engineering at Federal University of São João del-Rei, Brazil
   index: 2
date: 18 May 2020
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
---

# Summary

The field of System Identification (SI) aims to build mathematical models for static and dynamic behavior from experimental data [@Lju1987]. In particular, nonlinear system identification has become a central issue in the SI community, and from the 1950s onwards many methods have been proposed. In this respect, NARMAX (Nonlinear AutoRegressive Moving Average with eXogenous input) models are among the most well-documented and used model representation of dynamical systems [@Bil2013].

The NARMAX model was proposed by [@BL1981; @LB1985; @CB1989] and can be described as

\begin{equation}
y_k= \mathcal{F}[y_{k-1}, \dotsc, y_{k-n_y},x_{k-d}, x_{k-d-1}, \dotsc, x_{k-d-n_x} + e_{k-1}, \dotsc, e_{k-n_e}] + e_k,
\end{equation}

where $n_y\in \mathbb{N}^*$, $n_x \in \mathbb{N}$, $n_e \in \mathbb{N}$ , are the maximum lags for the system output and input respectively; $x_k \in \mathbb{R}^{n_x}$ is the system input and $y_k \in \mathbb{R}^{n_y}$ is the system output at discrete time $k \in \mathbb{N}^n$; $e_k \in \mathbb{R}^{n_e}$ represents uncertainties and possible noise at discrete time $k$. In this case, $\mathcal{F}$ is some nonlinear function of the input and output regressors and $d$ is a time delay typically set to $d=1$.

Although there are many possible approximations of $\mathcal{F}(\cdot)$ (e.g., Neural Networks, Fuzzy, Wavelet, Radial Basis Function), the power-form Polynomial NARMAX model is the most commonly used [@Bil2013; @KST2020]:

\begin{align}
  y_k = \sum_{i=1}^{p}\Theta_i \times \prod_{j=0}^{n_x}u_{k-j}^{b_i, j}\prod_{l=1}^{n_e}e_{k-l}^{d_i, l}\prod_{m=1}^{n_y}y_{k-m}^{a_i, m}
\label{eq5:narx}
\end{align}
where $p$ is the number of regressors, $\Theta_i$ are the model parameters, and $a_i, m$, $b_i, j$ and $d_i, l \in \mathbb{N}$ are the exponents of the output, input and noise terms, respectively.

The following example is a polynomial NARMAX model where the nonlinearity degree is equal to $2$, identified from experimental data of a DC motor/generator with no prior knowledge of the model form, taken from [@LJAM2017]:
\begin{align}
  y_k =& 1.7813y_{k-1}-0.7962y_{k-2}+0.0339x_{k-1} -0.1597x_{k-1} y_{k-1} +0.0338x_{k-2} + \nonumber \\
  & + 0.1297x_{k-1}y_{k-2} - 0.1396x_{k-2}y_{k-1}+ 0.1086x_{k-2}y_{k-2}+0.0085y_{k-2}^2 + 0.0247e_{k-1}e_{k-2}
  \label{eq5:dcmotor}
\end{align}

The $\Theta$ values are the coefficients of each term of the polynomial equation.

Polynomial basis functions are one of the most used representations of NARMAX models due to several interesting atrributes, such as [@Bil2013; @Agu2004]:

- All polynomial functions are smooth in $\mathbb{R}$.
- The Weierstrass approximation theorem [@Wei1885] states that any continuous real-valued function defined on a closed and bounded space $[a,b]$ can be uniformly approximated using a polynomial on that interval.
- They can describe several nonlinear dynamical systems [@Bil2013], including industrial processes, control systems, structural systems, economic and financial systems, biology, medicine, and social systems [@WMNL2019; @FWHM2003; @GGBW2016; @KGHK2003; @BBWL2018; @CER2001; @Bil2013; @Agu2004; @MA2016].
- Several algorithms have been developed for both structure selection and parameter estimation of polynomial NARMAX models.
- Polynomial NARMAX models can be used both for prediction and inference. The structure of polynomial NARMAX models are easy to interpret and can be related to the underlying system, which is not a trivial task when using, for example, neural or wavelet functions.

Estimating the parameters of NARMAX models is a simple task if the model structure is known *a priori*. However, usually there is no information on what terms one should include in the final model, and selecting the correct terms has to be part of the system identification procedure. Thus, the identification of NARMAX models is twofold: selecting the most significant regressors given a dictionary of candidate terms, which relies on model structure selection algorithms, and estimating their parameters.

# SysIdentPy

`SysIdentPy` is an open-source Python package for system identification using polynomial NARMAX models. The package can handle SISO (Single-Input Single-Output) and MISO (Multiple-Inputs Single-Output) NARMAX model identification and its variants such as NARX, NAR, ARMAX, ARX, and AR models. It provides various tools for both model structure selection and parameter estimation including classical algorithms, e.g., forward regression orthogonal least squares and extended least squares orthogonal forward regression; parameter estimation using ordinary least squares, recursive algorithms and adaptative filters; the Akaike Information Criterion (AIC), Bayesian Information Criterion (BIC), Khinchin's law of iterated logarithm criterion (LILC), and Final Prediction Error (FPE) methods for model order selection [@HK1999]; regression metrics; and residual analysis. The reader is referred to the package documentation for further details.

`SysIdentPy` is designed to be easily expanded and user friendly. Moreover, the package aims to provide useful tools for researchers and students not only in the SI field, but also in correlated areas such as Machine Learning, Statistical Learning and Data Science. Recently, an R package was published [@ayala2020r] with tools to model dynamic systems using NARMAX models. However, to the best of our knowledge, `SysIdentPy` is the first open-source package for system identification using NARMAX models in Python. Moreover, SysIdentPy includes recursive and gradient methods for parameter estimation, e.g., recursive least squares, affine least mean squares, sign-sign least mean squares and many others that are not available in the above-mentioned R package. Also, the user can choose between four methods for model order selection, which is not possible with the mentioned R package.

# Example

The following is an example of how to use `SysIdentPy` to build a NARMAX model from data. For simplicity, the example uses simulated data with $1000$ samples, generated using the method `get_miso_data`:

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sysidentpy.polynomial_basis import PolynomialNarmax
from sysidentpy.metrics import root_relative_squared_error
from sysidentpy.utils.generate_data import get_miso_data

x_train, x_valid, y_train, y_valid = get_miso_data(n=1000,
                                                   colored_noise=False,
                                                   sigma=0.001,
                                                   train_percentage=90)
```

Assuming that there is no information regarding what system generated the data, a dictionary of candidate terms must be created by defining the nonlinearity degree of the polynomial function and the maximum lag of the input and output terms. These parameters are, respectively, `non_degree, ylag, xlag`. The Akaike Information Criterion is chosen for model order selection and the least squares method is used for parameter estimation:

```python
model = PolynomialNarmax(non_degree=2,
                         order_selection=True,
                         ylag=2, xlag=[[1, 2], [1, 2]],
                         info_criteria='aic',
                         estimator='least_squares',
                         )
```

The user can also run a SISO example by replacing `get_miso_data` with `get_siso_data` and the `xlag` values with an integer or a list of integers. If one wants to estimate the parameters using, for example, the recursive least squares algorithm, just set `estimator` to `'recursive_least_squares'`. Replacing the AIC method with BIC, for example, can be done analogously by replacing `'aic'` with `'bic'`.

The `fit` method is used to obtain the model and `predict` to validate the model using new data.
The metric to evaluate is the relative root squared error. To get the root mean square error metric, for example, import it using `from sysidentpy.metrics import root_mean_square_error` and replace the root relative squared error method with it.

```python
model = PolynomialNarmax(non_degree=2,
                         order_selection=True,
                         ylag=2, xlag=[[1, 2], [1, 2]],
                         info_criteria='aic',
                         estimator='least_squares',
                         )

model.fit(x_train, y_train)
yhat = model.predict(x_valid, y_valid)
rrse = root_relative_squared_error(y_valid, yhat)
print(rrse)
```

The `model.results` and `model.residuals` statements return the polynomial model obtained using the `fit` method and plot the results for qualitative analysis.

```python
results = pd.DataFrame(model.results(err_precision=8,
                                     dtype='dec'),
                       columns=['Regressors', 'Parameters', 'ERR'])

print(results)
ee, ex, extras, lam = model.residuals(x_valid, y_valid, yhat)
model.plot_result(y_valid, yhat, ee, ex)
```

The table below and Figure 1 are the ouput of the aforementioned example. Table 1 details the regressors chosen to compose the final model, its respective parameters and the error reduction ratio (ERR), which measure the contribution of each regressor to explain the system output. ERR values can be interpreted as a feature importance metric. Figure 1 depicts the simulation of model prediction and the validation data as well as the autocorrelation of the model residues and the cross-correlation between the input and the residues.

| Regressors     | Parameters | ERR        |
|----------------|------------|------------|
| x2(k-1)        | 0.6000     | 0.90482955 |
| x2(k-2)x1(k-1) | -0.3000    | 0.05072675 |
| y(k-1)^2       | 0.3999     | 0.04410386 |
| x1(k-1)y(k-1)  | 0.1000     | 0.00033239 |

![\label{fig:example}](example1.png)
_Figure 1. Results from modeling a simulated system available with the `SysIdentPy` package. Free run simulation (validation data vs. model prediction), autocorrelation of the residues and cross correlation between residues and the input._

For more information and examples of how to build NARMAX models and its variants using different methods for parameters estimation, model order selection and many more, see the package documentation.

# Future work

Future releases will include new methods for model structure selection of polynomial NARMAX models, new basis functions, multiobjective model structure selection and parameter estimation algorithms, new adaptative filters, frequency domain analysis, and algorithms for using NARMAX models for classification problems.

# References
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps/code to reproduce:
1. Import '...'
2. Use the function '....'
3. See error

**Expected results**
A clear and concise description of what you expected to happen.

**Actual results**

**Environment**
Version of the packages you are using

**Screenshots**
If applicable, add screenshots to help explain your problem.


**Additional context**
Add any other context about the problem here.
---
name: Feature request
about: Suggest an idea for this project
title: ''
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
Contributing
============

SysIdentPy is intended to be a community project, hence all contributions are welcome!

Sugestions and Bug reporting
----------------------------
There exist many possible use cases in System Identification field
and we can not test all scenarios without your help! If you find any
bugs or have suggestions, please report them on 'issue tracker'_ on GitHub.

.. _`issue tracker`: https://github.com/wilsonrljr/sysidentpy/issues


Documentation
-------------

Documentation is as important as the library itself. English is not the primary language of the main authors, so if you find any typo or anything wrong do not hesitate to point out to us.

Development environment
-----------------------

These are some basic steps to help us with code:

1. Install and Setup Git on your computer.
3. `Fork sysidentpy <https://help.github.com/articles/fork-a-repo/>`_.
4. `Clone the fork on your local machine  <https://help.github.com/articles/cloning-a-repository/>`_.
5. Install it in development mode using
   :code:`pip install --editable /path/to/sysidentpy/[dev]`
6. Create a new branch.
7. Make changes following the coding style of the project (or suggesting improvements).
8. Run the tests.
9. Write and/or adapt existing test if needed.
10. Add documentation if needed.
11. Commit.
12. `Push to your fork <https://help.github.com/articles/pushing-to-a-remote/>`_.
13. `Open a pull request! <https://help.github.com/articles/creating-a-pull-request/>`_User Guide
==========

Presenting main functionality
-----------------------------

Example created by Wilson Rocha Lacerda Junior

Here we import the NARMAX model, the metric for model evaluation and the
methods to generate sample data for tests. Also, we import pandas for
specific usage.

.. code:: ipython3

    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from sysidentpy.model_structure_selection import FROLS
    from sysidentpy.basis_function import Polynomial, Fourier
    from sysidentpy.metrics import root_relative_squared_error
    from sysidentpy.utils.generate_data import get_siso_data
    from sysidentpy.utils.display_results import results
    from sysidentpy.utils.plotting import plot_residues_correlation, plot_results
    from sysidentpy.residues.residues_correlation import compute_residues_autocorrelation, compute_cross_correlation


Generating 1 input 1 output sample data
---------------------------------------

The data is generated by simulating the following model:

:math:`y_k = 0.2y_{k-1} + 0.1y_{k-1}x_{k-1} + 0.9x_{k-1} + e_{k}`

If ``colored_noise`` is set to True:

:math:`e_{k} = 0.8\nu_{k-1} + \nu_{k}`

where :math:`x` is a uniformly distributed random variable and
:math:`\nu` is a gaussian distributed variable with :math:`\mu=0` and
:math:`\sigma=0.1`

In the next example we will generate a data with 1000 samples with white
noise and selecting 90% of the data to train the model.

.. code:: ipython3

    x_train, x_valid, y_train, y_valid = get_siso_data(n=1000,
                                                       colored_noise=False,
                                                       sigma=0.001,
                                                       train_percentage=90)

To obtain a NARMAX model we have to choose some values, e.g, the nonlinearity degree (``degree``), the maximum lag for the inputs and output (``xlag`` and ``ylag``).

In addition, you can select the information criteria to be used with the Error Reduction Ratio to select the model order and the method to estimate the model parameters:

-  Information Criteria: aic, bic, lilc, fpe
-  Parameter Estimation: least_squares, total_least_squares,
   recursive_least_squares, least_mean_squares and many other (see the
   docs)

The ``n_terms`` values is optional. It refer to the number of terms to
included in the final model. You can set this value based on the
information criteria (see below) or based on priori information about
the model structure. The default value is ``n_terms=None``, so the
algorithm will choose the minimum value reached by the information
criteria.

To use information criteria you have to set ``order_selection=True``. You
can also select ``n_info_values`` (default = 15).

.. code:: ipython3

    basis_function = Polynomial(degree=2)
    model = FROLS(
        basis_function=basis_function,
        order_selection=True,
        n_info_values=10,
        extended_least_squares=False,
        ylag=2, xlag=2,
        info_criteria='aic',
        estimator='least_squares',
    )

The ``fit`` method executes the Error Reduction Ratio algorithm using Householder reflection to select the model structure.

.. code:: ipython3

    model.fit(X=x_train, y=y_train)




.. parsed-literal::

    <sysidentpy.polynomial_basis.narmax.PolynomialNarmax at 0x7f332d0cc490>



The ``predict`` method is use to generate the predictions. We support ``free run simulation`` (also known as ``infinity steps ahead``), one-step ahead and n-steps ahead prediction.

**Note:** **Free run simulation** means that the ``y`` values used for predictions are the ones predicted in previous iterations. In **one-step ahead** simulation, otherwise, the ``y`` values used are the observed values of the system. In **k-steps ahead**, the ``y`` values are the predicted values but at each *k* iterations the observed values are used.

.. code:: ipython3

    yhat = model.predict(X=x_valid, y=y_valid)

In this example we use the ``root_relative_squared_error`` metric because it is often used in System Identification. More metrics and information about it can be found on documentation.

.. code:: ipython3

    rrse = root_relative_squared_error(y_valid, yhat)
    print(rrse)



.. parsed-literal::

    0.0018758031321337446


``model_object.results`` return the selected model regressors, the estimated parameters and the ERR values. As shown below, the algorithm detect the exact model that was used for simulate the data.

.. code:: ipython3

    r = pd.DataFrame(
    results(
        model.final_model, model.theta, model.err,
        model.n_terms, err_precision=8, dtype='sci'
        ),
    columns=['Regressors', 'Parameters', 'ERR'])
    print(r)


.. parsed-literal::

    Regressors Parameters         ERR
    0        u1(k-2)     0.9001  0.95750813
    1         y(k-1)     0.2000  0.03916822
    2  u1(k-1)y(k-1)     0.1003  0.00332022


In addition, you can access the ``residuals analysis`` and ``plot_result`` methods to take a look at the prediction and residual analysis.

.. code:: ipython3

    plot_results(y=y_valid, yhat = yhat, n=1000)
    ee = compute_residues_autocorrelation(y_valid, yhat)
    plot_residues_correlation(data=ee, title="Residues", ylabel="$e^2$")
    x1e = compute_cross_correlation(y_valid, yhat, x_valid)
    plot_residues_correlation(data=x1e, title="Residues", ylabel="$x_1e$")




.. image:: output_16_0.svg
.. image:: output_16_1.svg
.. image:: output_16_2.svg


In the example above we let the number of terms to compose the final model to be defined as the minimum value of the information criteria. Once you ran the algorithm and choose the best number of parameters, you can turn ``order_selection`` to ``False`` and set the ``n_terms`` value (3 in this example). Here we have a small dataset, but in bigger data this can be critical because running information criteria algorithm is more computational expensive. Since we already know the best number of regressor, we set ``n_terms`` and we get the same result.

However, this is not only critical because computational efficiency. In many situation, the minimum value of the information criteria can lead to over fitting. In some cases, the difference between choosing a model with 30 regressors or 10 is minimal, so you can take the model with 10 terms without loosing accuracy.

In the following we use ``info_values`` to plot the information criteria values. As you can see, the minimum value relies where :math:`xaxis = 5`

.. code:: ipython3

    xaxis = np.arange(1, model.n_info_values + 1)
    plt.plot(xaxis, model.info_values)
    plt.xlabel('n_terms')
    plt.ylabel('Information Criteria')




.. parsed-literal::

    Text(0, 0.5, 'Information Criteria')




.. image:: output_18_1.svg


Important Note:
---------------

Here we are creating random samples with white noise and letting the algorithm choose the number of terms based on the minimum value of information criteria. This is not the best approach in System Identification, but serves as a simple example. The information criteria must be used as an **auxiliary tool** to select ``n_terms``. Plot the information values to help you on that!

If you run the example above several times you might find some cases where the algorithm choose only the first two regressors, or four (depending on the information criteria method selected). This is because the minimum value of information criteria depends on residual variance (affected by noise) and have some limitations in nonlinear scenarios. However, if you check the ERR values (robust to noise) you will see that the ERR is ordering the regressors in the correct way!

We have some examples on ``information_criteria`` notebook!

The ``n_info_values`` limits the number of regressors to apply the information criteria. We choose :math:`n_y = n_x = \ell = 2`, so the candidate regressor is a list of 15 regressors. We can set ``n_info_values = 15`` and see the information values for all regressors. This option can save some amount of computational resources when dealing with multiples inputs and large datasets.

.. code:: ipython3

    basis_function = Polynomial(degree=2)
    model = FROLS(
        basis_function=basis_function,
        order_selection=True,
        n_info_values=15,
        extended_least_squares=False,
        ylag=2, xlag=2,
        info_criteria='aic',
        estimator='least_squares',
        )

    model.fit(X=x_train, y=y_train)

    xaxis = np.arange(1, model.n_info_values + 1)
    plt.plot(xaxis, model.info_values)
    plt.xlabel('n_terms')
    plt.ylabel('Information Criteria')




.. parsed-literal::

    Text(0, 0.5, 'Information Criteria')




.. image:: output_21_1.svg


Now running without executing information criteria methods (setting the ``n_terms``) because we already know the optimal number of regressors

.. code:: ipython3

    basis_function = Polynomial(degree=2)
    model = FROLS(
        basis_function=basis_function,
        # order_selection=True,
        n_terms = 3,
        # n_info_values=15,
        extended_least_squares=False,
        ylag=2, xlag=2,
        info_criteria='aic',
        estimator='least_squares',
    )

    model.fit(X=x_train, y=y_train)
    yhat = model.predict(X=x_valid, y=y_valid)
    rrse = root_relative_squared_error(y_valid, yhat)
    print('rrse: ', rrse)

    r = pd.DataFrame(
    results(
        model.final_model, model.theta, model.err,
        model.n_terms, err_precision=8, dtype='sci'
        ),
    columns=['Regressors', 'Parameters', 'ERR'])
    print(r)


.. parsed-literal::

    rrse:  0.0018758031321337446

           Regressors Parameters         ERR
    0        u1(k-2)     0.9001  0.95750813
    1         y(k-1)     0.2000  0.03916822
    2  u1(k-1)y(k-1)     0.1003  0.00332022


You can access some extra information like the list of all candidate regressors

.. code:: ipython3

    # for now the list is returned as a codification. Here, $0$ is the constant term, $[1001]=y{k-1}, [100n]=y_{k-n}, [200n] = x1_{k-n}, [300n]=x2_{k-n}$ and so on
    model.regressor_code  # list of all possible regressors given non_degree, n_y and n_x values




.. parsed-literal::

    array([[   0,    0],
           [1001,    0],
           [1002,    0],
           [2001,    0],
           [2002,    0],
           [1001, 1001],
           [1002, 1001],
           [2001, 1001],
           [2002, 1001],
           [1002, 1002],
           [2001, 1002],
           [2002, 1002],
           [2001, 2001],
           [2002, 2001],
           [2002, 2002]])



.. code:: ipython3

    print(model.err, '\n\n')  # err values for the selected terms
    print(model.theta)  # estimated parameters for the final model structure


.. parsed-literal::

    [0.95750813 0.03916822 0.00332022 0.         0.         0.
     0.         0.         0.         0.         0.         0.
     0.         0.         0.        ]


    [[0.90008672]
     [0.19998879]
     [0.10026928]]
Code
====

.. toctree::
   :maxdepth: 2
   :caption: Contents:

sysidentpy base
===============
.. automodule:: sysidentpy.narmax_base
   :members:

sysidentpy FROLS
================
.. automodule:: sysidentpy.model_structure_selection.FROLS
   :members:

sysidentpy metamss
==================
.. automodule:: sysidentpy.model_structure_selection.MetaMSS
   :members:

sysidentpy simulation
=====================
.. automodule:: sysidentpy.simulation._simulation
   :members:

sysidentpy basis function
=========================
.. automodule:: sysidentpy.basis_function._basis_function
   :members:

sysidentpy narx_neural_network
==============================
.. automodule:: sysidentpy.narx_neural_network
   :members:

sysidentpy general_estimators
=============================
.. automodule:: sysidentpy.general_estimators.narx
   :members:

sysidentpy bpsogsa
==================
.. automodule:: sysidentpy.metaheuristics.bpsogsa
   :members:

sysidentpy residues
===================
.. automodule:: sysidentpy.residues.residues_correlation
   :members:
   :undoc-members:
   :show-inheritance:
   :private-members: # if you want to document __x attributes
   :special-members:

sysidentpy metrics
==================
.. automodule:: sysidentpy.metrics._regression
   :members:

sysidentpy estimators
=====================
.. automodule:: sysidentpy.parameter_estimation.estimators
   :members:

sysidentpy utils
================
.. automodule:: sysidentpy.utils._check_arrays
   :members:

sysidentpy generate data
========================
.. automodule:: sysidentpy.utils.generate_data
   :members:

sysidentpy plotting
===================
.. automodule:: sysidentpy.utils.plotting
   :members:

sysidentpy display results
==========================
.. automodule:: sysidentpy.utils.display_results
   :members:

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

A brief introduction to NARMAX models.
======================================

Author: Wilson Rocha Lacerda Junior

This is the first in a series of publications explaining a little bit about NARMAX models. I hope the content of these publications will help those who use or would like to use the SysIdentPy library.

**Note**: As I will use the term *Systems Identification* here and there, let me make a brief definition regarding these terms. Systems identification is one of the major areas that deals with the modeling of data-based processes. In this context, the term "system" can be interpreted as any set of operations that process one or more inputs and return one or more outputs. Examples include electrical systems, mechanical systems, biological systems, financial systems, chemical systems … literally anything you can relate to input and output data. The electricity demand is part of a system whose inputs can be, for example, quantity of the population, quantity of water in the reservoirs, season, events. The price of a property is the output of a system whose entries can be the city, per capita income, neighborhood, number of rooms, how old the house is, and many others. You got the idea.

Although there are many things related with Machine Learning, Statistical Learning and other fields,  each field has its particularities.

So, what is a NARMAX model?
---------------------------

You may have noticed the similarity between the acronym NARMAX with the well-known models ARX, ARMAX, etc., which are widely used for forecasting time series. And this resemblance is not by chance. The Autoregressive models with Moving Average and Exogenous Input (ARMAX) and their variations AR, ARX, ARMA (to name just a few) are one of the most used mathematical representations for identifying linear systems.

Let's go back to the model. I said that the *ARX* family of models is commonly used to model linear systems. Linear is the key word here. For nonlinear scenarios we have the *NARMAX* class (*Non-linear Autoregressive Models with Moving Average and Exogenous Input*). As reported by Billings (the creator of NARMAX model) in the book **Nonlinear System Identification: NARMAX Methods in the Time, Frequency, and Spatio-Temporal Domains**,  NARMAX started out as a model name, but soon became a philosophy when it comes to identifying nonlinear systems. Obtaining NARMAX models consists of performing the following steps:

- Dynamical tests and collecting data;
- Choice of mathematical representation;
- Detection of the model structure;
- Estimation of parameters;
- Validation;
- Analysis of the model.

We will cover each of these steps in further publications. The idea of this text is to present an overview of NARMAX models.

NARMAX models **are not**, however, a simple extension of ARMAX models. NARMAX models are able to represent the most different and complex nonlinear systems. Introduced in 1981 by the Electrical Engineer Stephen A. Billings, NARMAX models can be described as:

.. math::

    y_k= F^\ell[y_{k-1}, \dotsc, y_{k-n_y},x_{k-d}, x_{k-d-1}, \dotsc, x_{k-d-n_x} + e_{k-1}, \dotsc, e_{k-n_e}] + e_k

where :math:`n_y\in \mathbb{N}^*`, :math:`n_x \in \mathbb{N}`, :math:`n_e \in \mathbb{N}` , are the maximum lags for the system output and input respectively; :math:`x_k \in \mathbb{R}^{n_x}` is the system input and :math:`y_k \in \mathbb{R}^{n_y}` is the system output at discrete time :math:`k \in \mathbb{N}^n`; :math:`e_k \in \mathbb{R}^{n_e}` stands for uncertainties and possible noise at discrete time :math:`k`. In this case, :math:`\mathcal{F}^\ell` is some nonlinear function of the input and output regressors with nonlinearity degree :math:`\ell \in \mathbb{N}` and :math:`d` is a time delay typically set to :math:`d=1`.

If we do not include noise terms, :math:`e_{k-n_e}`, we have NARX models. If we set :math:`\ell = 1` then we deal with ARMAX models; if :math:`\ell = 1` and we do not include input and noise terms, it turns to AR model (ARX if we include inputs, ARMA if we include noise terms instead); if :math:`\ell>1` and there is no input terms, we have the NARMA. If there is no input or noise terms, we have NAR. There are several variants, but that is sufficient for now.

NARMAX representation
---------------------

There are several nonlinear functions representations to approximate the unknown mapping :math:`\mathrm{f}[\cdot]` in the NARMAX methods, e.g.,

- neural networks;
- fuzzy logic-based models;
- radial basis functions;
- wavelet basis;
- **polynomial basis**;
- generalized additive models;

The remainder of this post contemplates methods related to the power-form polynomial models, which is the most common used representation. Polynomial NARMAX is a mathematical model based on difference equations and relates the current output as a function of past inputs and outputs.

Polynomial NARMAX
-----------------

The polynomial NARMAX model with asymptotically stable equilibrium points can be described as:

.. math::

    y_k =& \sum_{0} + \sum_{i=1}^{p}\Theta_{y}^{i}y_{k-i} + \sum_{j=1}^{q}\Theta_{e}^{j}e_{k-j} + \sum_{m=1}^{r}\Theta_{x}^{m}x_{k-m}\\
    &+ \sum_{i=1}^{p}\sum_{j=1}^{q}\Theta_{ye}^{ij}y_{k-i} e_{k-j} + \sum_{i=1}^{p}\sum_{m=1}^{r}\Theta_{yx}^{im}y_{k-i} x_{k-m} \\
    &+ \sum_{j=1}^{q}\sum_{m=1}^{r}\Theta_{e x}^{jm}e_{k-j} x_{k-m} \\
    &+ \sum_{i=1}^{p}\sum_{j=1}^{q}\sum_{m=1}^{r}\Theta_{y e x}^{ijm}y_{k-i} e_{k-j} x_{k-m} \\
    &+ \sum_{m_1=1}^{r} \sum_{m_2=m_1}^{r}\Theta_{x^2}^{m_1 m_2} x_{k-m_1} x_{k-m_2} \dotsc \\
    &+ \sum_{m_1=1}^{r} \dotsc \sum_{m_l=m_{l-1}}^{r} \Theta_{x^l}^{m_1, \dotsc, m_2} x_{k-m_1} x_{k-m_l}

where :math:`\sum\nolimits_{0}`, :math:`c_{y}^{i}`, :math:`c_{e}^{j}`, :math:`c_{x}^{m}`, :math:`c_{y\e}^{ij}`, :math:`c_{yx}^{im}`, :math:`c_{e x}^{jm}`, :math:`c_{y e x}^{ijm}`, :math:`c_{x^2}^{m_1 m_2} \dotsc c_{x^l}^{m_1, \dotsc, ml}` are constant parameters.

Let's take a look at an example of a NARMAX model for an easy understanding. The following is a NARMAX model of degree~$2$, identified from experimental data of a DC motor/generator with no prior knowledge of the model form. If you want more information about the identification process, I wrote a paper comparing a polynomial NARMAX with a neural NARX model using that data (IN PORTUGUESE: Identificação de um motor/gerador CC por meio de modelos polinomiais autorregressivos e redes neurais artificiais)

.. math::

    y_k =& 1.7813y_{k-1}-0,7962y_{k-2}+0,0339x_{k-1} -0,1597x_{k-1} y_{k-1} +0,0338x_{k-2} \\
    & + 0,1297x_{k-1}y_{k-2} - 0,1396x_{k-2}y_{k-1}+ 0,1086x_{k-2}y_{k-2}+0,0085y_{k-2}^2 + 0.1938e_{k-1}e_{k-2}


But how those terms were selected? How the parameters were estimated? These questions will lead us to model structure selection and parameter estimation topics, but, for now,  let us discuss about those topics in a more simple manner.

First, the "structure" of a model is the set of terms (also called regressors) included in the final model. The parameters are the values multiplying each of theses terms. And looking at the example above we can notice an really important thing regarding polynomial NARMAX models dealt in this text: they have a non-linear structure, but they are linear-in-the-parameters. You will see how this note is important in the post about parameter estimation.

In this respect, consider the case where we have the input and output data of some system. For the sake of simplicity, suppose one input and one output. We have the data, but we do not know which lags to choose for the input or the output. Also, we know nothing about the system non-linearity. So, we have to define some values for maximum lags of the input, output and the noise terms, besides the choice of the :math:`\ell` value. It's worth to notice that many assumptions taken for linear cases are not valid in the nonlinear scenario and therefore select the maximum lags is not straightforward. So, how those values can make the modeling harder?

So we have one input and one output (disregard the noise terms for now). What if we choose the :math:`n_y = n_x = \ell = 2`? With these values, we have the following possibilities for compose the final model:

.. math::

    & constant, y_{k-1}, y_{k-2}, y_{k-1}^2, y_{k-2}^2, x_{k-1}, x_{k-2}, x_{k-1}^2, x_{k-2}^2,y_{k-1}y_{k-2},\\
    & y_{k-1}x_{k-1}, y_{k-1}x_{k-2}, y_{k-2}x_{k-1}, y_{k-2}x_{k-2}, x_{k-1}x_{k-2} .

So we have :math:`15` candidate terms to compose the final model.

Again, we do not know how of those terms are significant to compose the model. One should decide to use all the terms because there are only $15$. This, even in a simple scenario like this, can lead to a very wrong representation of the system that you are trying to modeling. Ok, what if we run a brute force algorithm to test the candidate regressors so we can select only the significant ones? In this case, we have :math:`2^{15} = 32768` possible model structures to be tested. You can think that it is ok, we have computer power for that. But this case is very simple and the system might have lags equal to :math:`10` for input and output. If we define :math:`n_y = n_x = 10` and :math:`\ell=2`, the number of possible models to be tested increases to :math:`2^{231}=3.4508732\times10^{69}`. If the non-linearity is set to :math:`3` then we have :math:`2^{1771} = 1.3308291989700907535925992... \times 10^{533}` candidate models.

Now, think about the case when we have not 1, but 5, 10 or more inputs... and have to include terms for the noise, and maximum lags are higher than 10... and nonlinearity is higher than 3...

And the problem is not solved by only identifying the most significant terms. How do you choose the number of terms to include in the final model. It is not just about check the relevance of each regressor, we have to think about the impact of including 5, 10 or 50 regressors in the model. And do not forget: after selecting the terms, we have to estimate its parameters.

As you can see, to select the most significant terms from a huge dictionary of possible terms is not an easy task. And it is hard not only because the complex combinatorial problem and the uncertainty concerning the model order. Identifying the most significant terms in a nonlinear scenario is very difficult because depends on the type of the non-linearity (sparse singularity or near-singular behavior, memory or dumping effects and many others), dynamical response (spatial-temporal systems, time-dependent), the steady-state response,  frequency of the data, the noise...

Despite all this complexity, NARMAX models are widely used because it is able to represent complex system with simple and transparent models, which terms are selected using robust algorithms for model structure selection. Model structure selection is the core of NARMAX methods and the scientific community is very active on improving classical methods and developing new ones. As I said, I will introduce some of those methods in another post.

I hope this publication served as a brief introduction to NARMAX models. Furthermore, I hope I have sparked your interest in this model class. The link to the other texts will be made available soon, but feel free to contact us if you are interested in collaborating with the SysIdentPy library or if you want to address any questions.
Contributing
============

SysIdentPy is intended to be a community project, hence all contributions are welcome!

Suggestions and Bug reporting
----------------------------
There exist many possible use cases in System Identification field
and we can not test all scenarios without your help! If you find any
bugs or have suggestions, please report them on 'issue tracker'_ on GitHub.

.. _`issue tracker`: https://github.com/wilsonrljr/sysidentpy/issues

Documentation
-------------

Documentation is as important as the library itself. English is not the primary language of the main authors, so if you find any typo or anything wrong do not hesitate to point out to us.

Development environment
-----------------------

These are some basic steps to help us with code:

1. Install and Setup Git on your computer.
3. `Fork sysidentpy <https://help.github.com/articles/fork-a-repo/>`_.
4. `Clone the fork on your local machine  <https://help.github.com/articles/cloning-a-repository/>`_.
5. Install it in development mode using
:code:`pip install --editable /path/to/sysidentpy/[dev]`

6. Create a new branch.
7. Make changes following the coding style of the project (or suggesting improvements).
8. Run the tests.
9. Write and/or adapt existing test if needed.
10. Add documentation if needed.
11. Commit.
12. `Push to your fork <https://help.github.com/articles/pushing-to-a-remote/>`_.
13. `Open a pull request! <https://help.github.com/articles/creating-a-pull-request/>`_
Examples v0.1.6
================

.. toctree::
   :maxdepth: 1
   :caption: V0.1.6 - Basic steps

   /examples-v016/basic_steps.ipynb

.. toctree::
   :maxdepth: 1
   :caption: V0.1.6 - Multiple Inputs usage

   /examples-v016/multiple_inputs_example.ipynb

.. toctree::
   :maxdepth: 1
   :caption: V0.1.6 - Information Criteria - Examples

   /examples-v016/information_criteria_examples.ipynb

.. toctree::
   :maxdepth: 1
   :caption: V0.1.6 - Using Extended Least Squares Algorithm

   /examples-v016/extended_least_squares.ipynb

.. toctree::
   :maxdepth: 1
   :caption: V0.1.6 - Setting specific lags

   /examples-v016/defining_lags.ipynb

.. toctree::
   :maxdepth: 1
   :caption: V0.1.6 - Parameter estimation methods

   /examples-v016/parameter_estimation.ipynb

.. toctree::
   :maxdepth: 1
   :caption: V0.1.6 - MetaMSS Algorithm

   /examples-v016/metamss.ipynb

.. toctree::
   :maxdepth: 1
   :caption: V0.1.6 - Accelerated Orthogonal Least-Squares algorithm

   /examples-v016/aols.ipynb

.. toctree::
   :maxdepth: 1
   :caption: V0.1.6 - F-16 aircraft

   /examples-v016/f_16_benchmark.ipynb

.. toctree::
   :maxdepth: 1
   :caption: V0.1.6 - NARX Neural Network

   /examples-v016/narx_neural_network.ipynb

.. toctree::
   :maxdepth: 1
   :caption: V0.1.6 - General Estimators

   /examples-v016/general_estimators.ipynb
   
.. toctree::
   :maxdepth: 1
   :caption: V0.1.6 - Simulate a predefined model

   /examples-v016/simulating_a_predefined_model.ipynb
   
.. toctree::
   :maxdepth: 1
   :caption: V0.1.6 - System Identification Using Adaptive Filters

   /examples-v016/system_identification_using_adaptative_filters.ipynb
   
.. toctree::
   :maxdepth: 1
   :caption: V0.1.6 - Identification of an electromechanical system

   /examples-v016/identification_of_an_electromechanical_system.ipynb
   
.. toctree::
   :maxdepth: 1
   :caption: V0.1.6 - N-steps-ahead simulation

   /examples-v016/n_steps_ahead_prediction.ipynb
Examples
========

.. toctree::
   :maxdepth: 1
   :caption: Basic steps

   /examples/basic_steps.ipynb

.. toctree::
   :maxdepth: 1
   :caption: Air passengers benchmark

   /examples/air_passenger_benchmark.ipynb

.. toctree::
   :maxdepth: 1
   :caption: Load forecasting benchmark

   /examples/load_forecasting_benchmark.ipynb

.. toctree::
   :maxdepth: 1
   :caption: PV forecasting benchmark

   /examples/PV_forecasting_benchmark.ipynb

.. toctree::
   :maxdepth: 1
   :caption: Multiple Inputs usage

   /examples/multiple_inputs_example.ipynb

.. toctree::
   :maxdepth: 1
   :caption: Information Criteria - Examples

   /examples/information_criteria_examples.ipynb

.. toctree::
   :maxdepth: 1
   :caption: Using Extended Least Squares Algorithm

   /examples/extended_least_squares.ipynb

.. toctree::
   :maxdepth: 1
   :caption: Setting specific lags

   /examples/defining_lags.ipynb

.. toctree::
   :maxdepth: 1
   :caption: Parameter estimation methods

   /examples/parameter_estimation.ipynb

.. toctree::
   :maxdepth: 1
   :caption: MetaMSS Algorithm

   /examples/metamss.ipynb

.. toctree::
   :maxdepth: 1
   :caption: Accelerated Orthogonal Least-Squares algorithm

   /examples/aols.ipynb

.. toctree::
   :maxdepth: 1
   :caption: F-16 aircraft

   /examples/f_16_benchmark.ipynb

.. toctree::
   :maxdepth: 1
   :caption: NARX Neural Network

   /examples/narx_neural_network.ipynb

.. toctree::
   :maxdepth: 1
   :caption: General Estimators

   /examples/general_estimators.ipynb
   
.. toctree::
   :maxdepth: 1
   :caption: Simulate a predefined model

   /examples/simulating_a_predefined_model.ipynb
     
.. toctree::
   :maxdepth: 1
   :caption: Identification of an electromechanical system

   /examples/identification_of_an_electromechanical_system.ipynb
   
.. toctree::
   :maxdepth: 1
   :caption: N-steps-ahead simulation

   /examples/n_steps_ahead_prediction.ipynb
Install Guide
=============

Requirements
------------

SysIdentPy requires:

- Python (>= 3.7)
- NumPy (>= 1.5.0) for all numerical algorithms
- Matplotlib >= 1.5.2 for static plotting and visualizations

==============   ===================
Platform         Status
==============   ===================
Linux            OK
Windows x64      OK
MacOS			 OK
==============   ===================

**SysIdentPy do not to support Python 2.7.**

A few examples require pandas >= 0.18.0. However, it is not required to use sysidentpy.

Installation
------------

The easiest way to get sysidentpy running is to install it using ``pip``   ::

    pip install sysidentpy

We will made it available at conda repository as soon as possible.
.. sys-identpy documentation master file, created by
   sphinx-quickstart on Sun Mar 15 08:37:08 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to SysIdentPy's documentation!
======================================

**SysIdentPy** is a Python module for System Identification using **NARMAX** models built on top of **numpy** and is distributed under the 3-Clause BSD license.

The NARMAX model is described as:

.. math::
	
	 y_k= F[y_{k-1}, \dotsc, y_{k-n_y},x_{k-d}, x_{k-d-1}, \dotsc, x_{k-d-n_x} + e_{k-1}, \dotsc, e_{k-n_e}] + e_k

where :math:`n_y\in \mathbb{N}^*`, :math:`n_x \in \mathbb{N}`, :math:`n_e \in \mathbb{N}`,
are the maximum lags for the system output and input respectively;
:math:`x_k \in \mathbb{R}^{n_x}` is the system input and :math:`y_k \in \mathbb{R}^{n_y}`
is the system output at discrete time :math:`k \in \mathbb{N}^n`;
:math:`e_k \in \mathbb{R}^{n_e}` stands for uncertainties and possible noise
at discrete time :math:`k`. In this case, :math:`\mathcal{F}` is some nonlinear function
of the input and output regressors and :math:`d` is a time delay typically set to :math:`d=1`.

.. note::
	The update **v0.1.7**  has been released with major changes and additional features.
	
	There are several API modifications and you will need to change your code to have the new (and upcoming) features.
	
	Check the examples of how to use the new version in the `documentation page <http://sysidentpy.org/notebooks.html>`__

	For more details, please see the `changelog <http://sysidentpy.org/changelog/v0.1.7.html>`__

.. seealso::
	The examples directory has several Jupyter notebooks presenting basic tutorials of how to use the package and some specific applications of **SysIdentPy**. `Try it out! <http://sysidentpy.org/notebooks.html>`__

.. tip::
	SysIdentPy now support NARX Neural Network and General estimators, e.g., sklearn estimators and Catboost. Check it out!

.. code-block:: python

	import numpy as np
	import pandas as pd
	import matplotlib.pyplot as plt
	from sysidentpy.metrics import mean_squared_error
	from sysidentpy.utils.generate_data import get_siso_data


	# Generate a dataset of a simulated dynamical system
	x_train, x_valid, y_train, y_valid = get_siso_data(
		n=1000,
		colored_noise=False,
        sigma=0.001,
        train_percentage=80
	)


Polynomial NARX
~~~~~~~~~~~~~~~

.. code-block:: python

	from sysidentpy.model_structure_selection import FROLS
	from sysidentpy.basis_function._basis_function import Polynomial
	from sysidentpy.utils.display_results import results
	from sysidentpy.utils.plotting import plot_residues_correlation, plot_results
	from sysidentpy.residues.residues_correlation import compute_residues_autocorrelation, compute_cross_correlation
	
	basis_function = Polynomial(degree=2)
	model = FROLS(
		order_selection=True,
		n_info_values=10,
		extended_least_squares=False,
		ylag=2,
		xlag=2,
		info_criteria='aic',
		estimator='least_squares',
		basis_function=basis_function
	)
	model.fit(X=x_train, y=y_train)
	yhat = model.predict(X=x_valid, y=y_valid)
	rrse = root_relative_squared_error(y_valid, yhat)
	print(rrse)
	r = pd.DataFrame(
		results(
			model.final_model, model.theta, model.err,
			model.n_terms, err_precision=8, dtype='sci'
			),
		columns=['Regressors', 'Parameters', 'ERR'])
	print(r)
	
	Regressors     Parameters        ERR
	0        x1(k-2)     0.9000  0.95556574
	1         y(k-1)     0.1999  0.04107943
	2  x1(k-1)y(k-1)     0.1000  0.00335113

	plot_results(y=y_valid, yhat=yhat, n=1000)
	ee = compute_residues_autocorrelation(y_valid, yhat)
	plot_residues_correlation(data=ee, title="Residues", ylabel="$e^2$")
	x1e = compute_cross_correlation(y_valid, yhat, x2_val)
	plot_residues_correlation(data=x1e, title="Residues", ylabel="$x_1e$")



.. image:: ../../examples/figures/polynomial_narmax.png

NARX Neural Network
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

	from torch import nn
	from sysidentpy.neural_network import NARXNN

	class NARX(nn.Module):
		def __init__(self):
			super().__init__()
			self.lin = nn.Linear(4, 10)
			self.lin2 = nn.Linear(10, 10)
			self.lin3 = nn.Linear(10, 1)
			self.tanh = nn.Tanh()

		def forward(self, xb):
			z = self.lin(xb)
			z = self.tanh(z)
			z = self.lin2(z)
			z = self.tanh(z)
			z = self.lin3(z)
			return z

	narx_net = NARXNN(
		net=NARX(),
		ylag=2,
		xlag=2,
		loss_func='mse_loss',
		optimizer='Adam',
		epochs=200,
		verbose=False,
		optim_params={'betas': (0.9, 0.999), 'eps': 1e-05} # optional parameters of the optimizer
	)

	train_dl = narx_net.data_transform(x_train, y_train)
	valid_dl = narx_net.data_transform(x_valid, y_valid)
	narx_net.fit(train_dl, valid_dl)
	yhat = narx_net.predict(x_valid, y_valid)
	ee, ex, extras, lam = narx_net.residuals(x_valid, y_valid, yhat)
	narx_net.plot_result(y_valid, yhat, ee, ex)


.. image:: ../../examples/figures/narx_network.png

Catboost-narx
~~~~~~~~~~~~~

.. code-block:: python

	from sysidentpy.general_estimators import NARX
	from catboost import CatBoostRegressor

	catboost_narx = NARX(
		base_estimator=CatBoostRegressor(
			iterations=300,
			learning_rate=0.1,
			depth=6),
		xlag=2,
		ylag=2,
		fit_params={'verbose': False}
	)

	catboost_narx.fit(x_train, y_train)
	yhat = catboost_narx.predict(x_valid, y_valid)
	ee, ex, extras, lam = catboost_narx.residuals(x_valid, y_valid, yhat)
	catboost_narx.plot_result(y_valid, yhat, ee, ex)


.. image:: ../../examples/figures/catboost_narx.png

Catboost without NARX configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following is the Catboost performance *without* the NARX configuration.

.. code-block:: python

	def plot_results(yvalid, yhat):
		_, ax = plt.subplots(figsize=(14, 8))
		ax.plot(y_valid[:200], label='Data', marker='o')
		ax.plot(yhat[:200], label='Prediction', marker='*')
		ax.set_xlabel("$n$", fontsize=18)
		ax.set_ylabel("$y[n]$", fontsize=18)
		ax.grid()
		ax.legend(fontsize=18)
		plt.show()

	catboost = CatBoostRegressor(
		iterations=300,
		learning_rate=0.1,
		depth=6
	)

	catboost.fit(x_train, y_train, verbose=False)
	plot_results(y_valid, catboost.predict(x_valid))


.. image:: ../../examples/figures/catboost.png


Changelog
---------

See the `changelog <http://sysidentpy.org/changelog/v0.1.7.html>`__
for a history of notable changes to **SysIdentPy**.


Development
-----------

We welcome new contributors of all experience levels. The **SysIdentPy** community goals are to be helpful, welcoming, and effective.

.. note::
	We use the `pytest` package for testing. The test functions are located in tests subdirectories at each folder inside **SysIdentPy**, which check the validity of the algorithms.

Run the `pytest` in the respective folder to perform all the tests of the corresponding sub-packages.

Currently, we have around 81% of code coverage.

You can install pytest using ::

	pip install -U pytest

Example of how to run the tests:
--------------------------------

Open a terminal emulator of your choice and go to a subdirectory, e.g, ::

	\sysidentpy\metrics\

Just type :code:`pytest` and you get a result like ::


	========== test session starts ==========

	platform linux -- Python 3.7.6, pytest-5.4.2, py-1.8.1, pluggy-0.13.1

	rootdir: ~/sysidentpy

	plugins: cov-2.8.1

	collected 12 items

	tests/test_regression.py ............ [100%]

	========== 12 passed in 2.45s ==================

You can also see the code coverage using the :code:`pytest-cov` package. First, install :code:`pytest-cov` using ::

	pip install pytest-cov

Run the command below in the **SysIdentPy** root directory, to generate the report. ::

	pytest --cov=.


Source code
-----------

You can check the latest sources with the command::

    git clone https://github.com/wilsonrljr/sysidentpy.git

Project History
---------------

The project was started by Wilson R. L. Junior, Luan Pascoal and Samir A. M. Martins as a project for System Identification discipline. Samuel joined early in 2019.

The project is actively maintained by Wilson R. L. Junior and looking for contributors.

Communication
-------------

- Discord server: https://discord.gg/8eGE3PQ
- Website(soon): http://sysidentpy.org

Citation
--------

If you use **SysIdentPy** on your project, please `drop me a line <wilsonrljr@outlook.com>`__.

If you use **SysIdentPy** on your scientific publication, we would appreciate citations to the following paper:

- Lacerda et al., (2020). SysIdentPy: A Python package for System Identification using NARMAX models. Journal of Open Source Software, 5(54), 2384, https://doi.org/10.21105/joss.02384 ::

	@article{Lacerda2020,
	  doi = {10.21105/joss.02384},
	  url = {https://doi.org/10.21105/joss.02384},
	  year = {2020},
	  publisher = {The Open Journal},
	  volume = {5},
	  number = {54},
	  pages = {2384},
	  author = {Wilson Rocha Lacerda Junior and Luan Pascoal Costa da Andrade and Samuel Carlos Pessoa Oliveira and Samir Angelo Milani Martins},
	  title = {SysIdentPy: A Python package for System Identification using NARMAX models},
	  journal = {Journal of Open Source Software}
	}

Inspiration
-----------

The documentation and structure (even this section) is openly inspired by sklearn, einsteinpy, and many others as we used (and keep using) them to learn.

Contents
--------

.. toctree::
    :maxdepth: 1

    installation
    introduction_to_narmax
    user_guide
    dev_guide
    notebooks
    notebooksv016
    changelog/v0.1.8
    code
Changes in SysIdentPy
=====================

v0.1.8
------

CONTRIBUTORS
~~~~~~~~~~~~

- wilsonrljr

CHANGES
~~~~~~~
- The update **v0.1.8**  has been released with additional feature, minor API changes and fixes of the new features added in v0.1.7. 

- MAJOR: Ensemble Basis Functions
    - Now you can use different basis function together. For now we allow to use Fourier combined with Polynomial of different degrees. 

- API change: Add "ensemble" parameter in basis function to combine the features of different basis function.

- Fix: N-steps ahead prediction for model_type="NAR" is working properly now with different forecast horizon.

- DOC: Air passenger benchmark
    - Remove unused code.
    - Use default hyperparameter in SysIdentPy models.

- DOC: Load forecasting benchmark
    - Remove unused code.
    - Use default hyperparameter in SysIdentPy models.

- DOC: PV forecasting benchmark
    - Remove unused code.
    - Use default hyperparameter in SysIdentPy models.

v0.1.7
------

CONTRIBUTORS
~~~~~~~~~~~~

- wilsonrljr

CHANGES
~~~~~~~
- The update **v0.1.7**  has been released with major changes and additional features. There are several API modifications and you will need to change your code to have the new (and upcoming) features. All modifications are meant to make future expansion easier.

- On the user's side, the changes are not that disruptive, but in the background there are many changes that allowed the inclusion of new features and bug fixes that would be complex to solve without the changes. Check the `documentation page <http://sysidentpy.org/notebooks.html>`__

- Many classes were basically rebuild it from scratch, so I suggest to look at the new examples of how to use the new version.

- I will present the main updates below in order to highlight features and usability and then all API changes will be reported.

- MAJOR: NARX models with Fourier basis function `Issue63 <https://github.com/wilsonrljr/sysidentpy/issues/63>`__, `Issue64 <https://github.com/wilsonrljr/sysidentpy/issues/64>`__
    - The user can choose which basis they want by importing it from sysidentpy.basis_function. Check the notebooks with examples of how to use it.
    - Polynomial and Fourier are supported for now. New basis functions will be added in next releases.

- MAJOR: NAR models `Issue58 <https://github.com/wilsonrljr/sysidentpy/issues/58>`__
    - It was already possible to build Polynomial NAR models, but with some hacks. Now the user just need to pass model_type="NAR" to build NAR models.
    - The user doesn't need to pass a vector of zeros as input anymore.
    - Works for any model structure selection algorithm (FROLS, AOLS, MetaMSS)

- Major: NFIR models `Issue59 <https://github.com/wilsonrljr/sysidentpy/issues/59>`__
    - NFIR models are models where the output depends only on past inputs. It was already possible to build Polynomial NFIR models, but with a lot of code on the user's side (much more than NAR, btw). Now the user just need to pass model_type="NFIR" to build NFIR models.
    - Works for any model structure selection algorithm (FROLS, AOLS, MetaMSS)

- Major: Select the order for the residues lags to use in Extended Least Squares - elag
    - The user can select the maximum lag of the residues to be used in the Extended Least Squares algorithm. In previous versions sysidentpy used a predefined subset of residual lags.
    - The degree of the lags follows the degree of the basis function

- Major: Residual analysis methods `Issue60 <https://github.com/wilsonrljr/sysidentpy/issues/60>`__
    - There are now specific functions to calculate the autocorrelation of the residuals and cross-correlation for the analysis of the residuals. In previous versions the calculation was limited to just two inputs, for example, limiting user usability.

- Major: Plotting methods `Issue61 <https://github.com/wilsonrljr/sysidentpy/issues/61>`__
    - The plotting functions are now separated from the models objects, so there are more flexibility regarding what to plot.
    - Residual plots were separated from the forecast plot

- API Change: sysidentpy.polynomial_basis.PolynomialNarmax is deprecated. Use sysidentpy.model_structure_selection.FROLS instead. `Issue64 <https://github.com/wilsonrljr/sysidentpy/issues/62>`__
    - Now the user doesn't need to pass the number of inputs as a parameter.
    - Added the elag parameter for unbiased_estimator. Now the user can define the number of lags of the residues for parameter estimation using the Extended Least Squares algorithm.
    - model_type parameter: now the user can select the model type to be built. The options are "NARMAX", "NAR" and "NFIR". "NARMAX" is the default. If you want to build a NAR model without any "hack", just set model_type="NAR". The same for "NFIR" models.

- API Change: sysidentpy.polynomial_basis.MetaMSS is deprecated. Use sysidentpy.model_structure_selection.MetaMSS instead. `Issue64 <https://github.com/wilsonrljr/sysidentpy/issues/64>`__
    - Now the user doesn't need to pass the number of inputs as a parameter.
    - Added the elag parameter for unbiased_estimator. Now the user can define the number of lags of the residues for parameter estimation using the Extended Least Squares algorithm.

- API Change: sysidentpy.polynomial_basis.AOLS is deprecated. Use sysidentpy.model_structure_selection.AOLS instead. `Issue64 <https://github.com/wilsonrljr/sysidentpy/issues/64>`__

- API Change: sysidentpy.polynomial_basis.SimulatePolynomialNarmax is deprecated. Use sysidentpy.simulation.SimulateNARMAX instead.

- API Change: Introducing sysidentpy.basis_function. Because NARMAX models can be built on different basis function, a new module is added to make easier to implement new basis functions in future updates `Issue64 <https://github.com/wilsonrljr/sysidentpy/issues/64>`__.
    - Each basis function class must have a fit and predict method to be used in training and prediction respectively. 

- API Change: unbiased_estimator method moved to Estimators class.
    - added elag option
    - change the build_information_matrix method to build_output_matrix

- API Change (new): sysidentpy.narmax_base
    - This is the new base for building NARMAX models. The classes have been rewritten to make it easier to expand functionality.

- API Change (new): sysidentpy.narmax_base.GenerateRegressors
    - create_narmax_code: Creates the base coding that allows representation for the NARMAX, NAR, and NFIR models.
    - regressor_space: Creates the encoding representation for the NARMAX, NAR, and NFIR models.

- API Change (new): sysidentpy.narmax_base.ModelInformation
    - _get_index_from_regressor_code: Get the index of the model code representation in regressor space.
    - _list_output_regressor_code: Create a flattened array of output regressors.
    - _list_input_regressor_code: Create a flattened array of input regressors.
    - _get_lag_from_regressor_code: Get the maximum lag from array of regressors.
    - _get_max_lag_from_model_code: the name says it all.
    - _get_max_lag: Get the maximum lag from ylag and xlag.

- API Change (new): sysidentpy.narmax_base.InformationMatrix
    - _create_lagged_X: Create a lagged matrix of inputs without combinations.
    - _create_lagged_y: Create a lagged matrix of the output without combinations.
    - build_output_matrix: Build the information matrix of output values.
    - build_input_matrix: Build the information matrix of input values.
    - build_input_output_matrix: Build the information matrix of input and output values.

- API Change (new): sysidentpy.narmax_base.ModelPrediction
    - predict: base method for prediction. Support infinity_steps ahead, one-step ahead and n-steps ahead prediction and any basis function.
    - _one_step_ahead_prediction: Perform the 1-step-ahead prediction for any basis function.
    - _n_step_ahead_prediction: Perform the n-step-ahead prediction for polynomial basis.
    - _model_prediction: Perform the infinity-step-ahead prediction for polynomial basis.
    - _narmax_predict: wrapper for NARMAX and NAR models.
    - _nfir_predict: wrapper for NFIR models.
    - _basis_function_predict: Perform the infinity-step-ahead prediction for basis functions other than polynomial.
    - basis_function_n_step_prediction: Perform the n-step-ahead prediction for basis functions other than polynomial.

- API Change (new): sysidentpy.model_structure_selection.FROLS `Issue62 <https://github.com/wilsonrljr/sysidentpy/issues/62>`__, `Issue64 <https://github.com/wilsonrljr/sysidentpy/issues/64>`__
    - Based on the old sysidentpy.polynomial_basis.PolynomialNARMAX. The class has been rebuilt with new functions and optimized code.
    - Enforcing keyword-only arguments. This is an effort to promote clear and non-ambiguous use of the library.
    - Add support for new basis functions.
    - The user can choose the residual lags.
    - No need to pass the number of inputs anymore.
    - Improved docstring.
    - Fixed minor grammatical and spelling mistakes.
    - New prediction method.
    - many under the hood changes.
 
- API Change (new): sysidentpy.model_structure_selection.MetaMSS `Issue64 <https://github.com/wilsonrljr/sysidentpy/issues/64>`__
    - Based on the old sysidentpy.polynomial_basis.MetaMSS. The class has been rebuilt with new functions and optimized code.
    - Enforcing keyword-only arguments. This is an effort to promote clear and non-ambiguous use of the library.
    - The user can choose the residual lags.
    - Extended Least Squares support.
    - Add support for new basis functions.
    - No need to pass the number of inputs anymore.
    - Improved docstring.
    - Fixed minor grammatical and spelling mistakes.
    - New prediction method.
    - many under the hood changes.

- API Change (new): sysidentpy.model_structure_selection.AOLS `Issue64 <https://github.com/wilsonrljr/sysidentpy/issues/64>`__
    - Based on the old sysidentpy.polynomial_basis.AOLS. The class has been rebuilt with new functions and optimized code.
    - Enforcing keyword-only arguments. This is an effort to promote clear and non-ambiguous use of the library.
    - Add support for new basis functions.
    - No need to pass the number of inputs anymore.
    - Improved docstring.
    - Change "l" parameter to "L".
    - Fixed minor grammatical and spelling mistakes.
    - New prediction method.
    - many under the hood changes.

- API Change (new): sysidentpy.simulation.SimulateNARMAX
    - Based on the old sysidentpy.polynomial_basis.SimulatePolynomialNarmax. The class has been rebuilt with new functions and optimized code.
    - Fix the Extended Least Squares support.
    - Fix n-steps ahead prediction and 1-step ahead prediction.
    - Enforcing keyword-only arguments. This is an effort to promote clear and non-ambiguous use of the library.
    - The user can choose the residual lags.
    - Improved docstring.
    - Fixed minor grammatical and spelling mistakes.
    - New prediction method.
    - Do not inherit from the structure selection algorithm anymore, only from narmax_base. Avoid circular import and other issues.
    - many under the hood changes.

- API Change (new): sysidentpy.residues
    - compute_residues_autocorrelation: the name says it all.
    - calculate_residues: get the residues from y and yhat.
    - get_unnormalized_e_acf: compute the unnormalized autocorrelation of the residues.
    - compute_cross_correlation: compute cross correlation between two signals.
    - _input_ccf
    - _normalized_correlation: compute the normalized correlation between two signals.

- API Change (new): sysidentpy.utils.plotting
    - plot_results: plot the forecast
    - plot_residues_correlation: the name says it all.

- API Change (new): sysidentpy.utils.display_results
    - results: return the model regressors, estimated parameter and ERR index of the fitted model in a table.

- DOC: Air passenger benchmark
    - Added notebook with Air passenger forecasting benchmark.
    - We compare SysIdentPy against prophet, neuralprophet, autoarima, tbats and many more.

- DOC: Load forecasting benchmark
    - Added notebook with load forecasting benchmark.

- DOC: PV forecasting benchmark
    - Added notebook with PV forecasting benchmark.

- DOC: Presenting main functionality
    - Example rewritten following the new api.
    - Fixed minor grammatical and spelling mistakes.

- DOC: Multiple Inputs usage
    - Example rewritten following the new api
    - Fixed minor grammatical and spelling mistakes.

- DOC: Information Criteria - Examples
    - Example rewritten following the new api.
    - Fixed minor grammatical and spelling mistakes.

- DOC: Important notes and examples of how to use Extended Least Squares
    - Example rewritten following the new api.
    - Fixed minor grammatical and spelling mistakes.

- DOC: Setting specific lags
    - Example rewritten following the new api.
    - Fixed minor grammatical and spelling mistakes.

- DOC: Parameter Estimation
    - Example rewritten following the new api.
    - Fixed minor grammatical and spelling mistakes.

- DOC: Using the Meta-Model Structure Selection (MetaMSS) algorithm for building Polynomial NARX models
    - Example rewritten following the new api.
    - Fixed minor grammatical and spelling mistakes.

- DOC: Using the Accelerated Orthogonal Least-Squares algorithm for building Polynomial NARX models
    - Example rewritten following the new api.
    - Fixed minor grammatical and spelling mistakes.

- DOC: Example: F-16 Ground Vibration Test benchmark
    - Example rewritten following the new api.
    - Fixed minor grammatical and spelling mistakes.

- DOC: Building NARX Neural Network using Sysidentpy
    - Example rewritten following the new api.
    - Fixed minor grammatical and spelling mistakes.

- DOC: Building NARX models using general estimators
    - Example rewritten following the new api.
    - Fixed minor grammatical and spelling mistakes.

- DOC: Simulate a Predefined Model
    - Example rewritten following the new api.
    - Fixed minor grammatical and spelling mistakes.

- DOC: System Identification Using Adaptive Filters
    - Example rewritten following the new api.
    - Fixed minor grammatical and spelling mistakes.

- DOC: Identification of an electromechanical system
    - Example rewritten following the new api.
    - Fixed minor grammatical and spelling mistakes.

- DOC: Example: N-steps-ahead prediction - F-16 Ground Vibration Test benchmark
    - Example rewritten following the new api.
    - Fixed minor grammatical and spelling mistakes.

- DOC: Introduction to NARMAX models
    - Fixed grammatical and spelling mistakes.



v0.1.6
------

CONTRIBUTORS
~~~~~~~~~~~~

- wilsonrljr

CHANGES
~~~~~~~

- MAJOR: Meta-Model Structure Selection Algorithm (Meta-MSS).
    - A new method for build NARMAX models based on metaheuristics. The algorithm uses a Binary hybrid Particle Swarm Optimization and Gravitational Search Algorithm with a new cost function to build parsimonious models.
    
    - New class for the BPSOGSA algorithm. New algorithms can be adapted in the Meta-MSS framework.
	
    - Future updates will add NARX models for classification and multiobjective model structure selection.

- MAJOR: Accelerated Orthogonal Least-Squares algorithm.
    - Added the new class AOLS to build NARX models using the Accelerated Orthogonal Least-Squares algorithm.
    
    - At the best of my knowledge, this is the first time this algorithm is used in the NARMAX framework. The tests I've made are promising, but use it with caution until the results are formalized into a research paper.

- Added notebook with a simple example of how to use MetaMSS and a simple model comparison of the Electromechanical system.

- Added notebook with a simple example of how to use AOLS

- Added ModelInformation class. This class have methods to return model information such as max_lag of a model code.
    - added _list_output_regressor_code
    - added _list_input_regressor_code
    - added _get_lag_from_regressor_code
    - added _get_max_lag_from_model_code

- Minor performance improvement: added the argument "predefined_regressors" in build_information_matrix function on base.py
    to improve the performance of the Simulation method.

- Pytorch is now an optional dependency. Use pip install sysidentpy['full'] 

- Fix code format issues.

- Fixed minor grammatical and spelling mistakes.

- Fix issues related to html on Jupyter notebooks examples on documentation.

- Updated Readme with examples of how to use.

- Improved descriptions and comments in methods.

- metaheuristics.bpsogsa (detailed description on code docstring)
    - added evaluate_objective_function
    - added optimize
    - added generate_random_population
    - added mass_calculation
    - added calculate_gravitational_constant
    - added calculate_acceleration
    - added update_velocity_position

- FIX issue #52


v0.1.5
------

CONTRIBUTORS
~~~~~~~~~~~~

- wilsonrljr

CHANGES
~~~~~~~

- MAJOR: n-steps-ahead prediction.
    - Now you can define the numbers of steps ahead in the predict function.
	- Only for Polynomial models for now. Next update will bring this functionality to Neural NARX and General Estimators.

- MAJOR: Simulating predefined models.
    - Added the new class SimulatePolynomialNarmax to handle the simulation of known model structures.
    - Now you can simulate predefined models by just passing the model structure codification. Check the notebook examples.

- Added 4 new notebooks in the example section.

- Added iterative notebooks. Now you can run the notebooks in Jupyter notebook section of the documentation in Colab.

- Fix code format issues.

- Added new tests for SimulatePolynomialNarmax and generate_data.

- Started changes related to numpy 1.19.4 update. There are still some Deprecation warnings that will be fixed in next update.

- Fix issues related to html on Jupyter notebooks examples on documentation.

- Updated Readme with examples of how to use.



v0.1.4
------

CONTRIBUTORS
~~~~~~~~~~~~

- wilsonrljr

CHANGES
~~~~~~~

- MAJOR: Introducing NARX Neural Network in SysIdentPy.
    - Now you can build NARX Neural Network on SysIdentPy.
    - This feature is built on top of Pytorch. See the docs for more details and examples of how to use.

- MAJOR: Introducing general estimators in SysIdentPy.
    - Now you are able to use any estimator that have Fit/Predict methods (estimators from Sklearn and Catboost, for example) and build NARX models based on those estimators.
    - We use the core functions of SysIdentPy and keep the Fit/Predict approach from those estimators to keep the process easy to use.
    - More estimators are coming soon like XGboost.

- Added notebooks to show how to build NARX neural Network.

- Added notebooks to show how to build NARX models using general estimators.

- Changed the default parameters of the plot_results function.

- NOTE: We will keeping improving the Polynomial NARX models (new model structure selection algorithms and multiobjective identification
is on our roadmap). These recent modifications will allow us to introduce new NARX models like PWARX models very soon.

- New template for the documentation site.

- Fix issues related to html on Jupyter notebooks examples on documentation.

- Updated Readme with examples of how to use.


v0.1.3
------

CONTRIBUTORS
~~~~~~~~~~~~

- wilsonrljr
- renard162

CHANGES
~~~~~~~

- Fixed a bug concerning the xlag and ylag in multiple input scenarios.
- Refactored predict function. Improved performance up to 87% depending on the number of regressors.
- You can set lags with different size for each input.
- Added a new function to get the max value of xlag and ylag. Work with int, list, nested lists.
- Fixed tests for information criteria.
- Added SysIdentPy logo.
- Refactored code of all classes following PEP 8 guidelines to improve readability.
- Added Citation information on Readme.
- Changes on information Criteria tests.
- Added workflow to run the tests when merge branch into master.
- Added new site domain.
- Updated docs.

Changes in SysIdentPy
=====================

v0.1.3
------

CONTRIBUTORS
~~~~~~~~~~~~

- wilsonrljr
- renard162

CHANGES
~~~~~~~

- Fixed a bug concerning the xlag and ylag in multiple input scenarios.
- Refactored predict function. Improved performance up to 87% depending on the number of regressors.
- You can set lags with different size for each input.
- Added a new function to get the max value of xlag and ylag. Work with int, list, nested lists.
- Fixed tests for information criteria.
- Added SysIdentPy logo.
- Refactored code of all classes following PEP 8 guidelines to improve readability.
- Added Citation information on Readme.
- Changes on information Criteria tests.
- Added workflow to run the tests when merge branch into master.
- Added new site domain.
- Updated docs.

v0.1.4
------

CONTRIBUTORS
~~~~~~~~~~~~

- wilsonrljr

CHANGES
~~~~~~~

- MAJOR: Introducing NARX Neural Network in SysIdentPy.
    - Now you can build NARX Neural Network on SysIdentPy.
    - This feature is built on top of Pytorch. See the docs for more details and examples of how to use.

- MAJOR: Introducing general estimators in SysIdentPy.
    - Now you are able to use any estimator that have Fit/Predict methods (estimators from Sklearn and Catboost, for example) and build NARX models based on those estimators.
    - We use the core functions of SysIdentPy and keep the Fit/Predict approach from those estimators to keep the process easy to use.
    - More estimators are coming soon like XGboost.

- Added notebooks to show how to build NARX neural Network.

- Added notebooks to show how to build NARX models using general estimators.

- Changed the default parameters of the plot_results function.

- NOTE: We will keeping improving the Polynomial NARX models (new model structure selection algorithms and multiobjective identification
is on our roadmap). These recent modifications will allow us to introduce new NARX models like PWARX models very soon.

- New template for the documentation site.

- Fix issues related to html on Jupyter notebooks examples on documentation.

- Updated Readme with examples of how to use.


v0.1.5
------
CONTRIBUTORS
CHANGES
Changes in SysIdentPy
=====================

v0.1.9
------

CONTRIBUTORS
~~~~~~~~~~~~

- wilsonrljr

CHANGES
~~~~~~~
- The update **v0.1.9**  has been released with additional feature, minor API changes and fixes of the new features added in v0.1.7. 

- DOC: PV forecasting benchmark
    - FIX AOLS prediction. The example was using the meta_mss model in prediction, so the results for AOLS were wrong.


v0.1.8
------

CONTRIBUTORS
~~~~~~~~~~~~

- wilsonrljr

CHANGES
~~~~~~~
- The update **v0.1.8**  has been released with additional feature, minor API changes and fixes of the new features added in v0.1.7. 

- MAJOR: Ensemble Basis Functions
    - Now you can use different basis function together. For now we allow to use Fourier combined with Polynomial of different degrees. 

- API change: Add "ensemble" parameter in basis function to combine the features of different basis function.

- Fix: N-steps ahead prediction for model_type="NAR" is working properly now with different forecast horizon.

- DOC: Air passenger benchmark
    - Remove unused code.
    - Use default hyperparameter in SysIdentPy models.

- DOC: Load forecasting benchmark
    - Remove unused code.
    - Use default hyperparameter in SysIdentPy models.

- DOC: PV forecasting benchmark
    - Remove unused code.
    - Use default hyperparameter in SysIdentPy models.


v0.1.7
------

CONTRIBUTORS
~~~~~~~~~~~~

- wilsonrljr

CHANGES
~~~~~~~
- The update **v0.1.7**  has been released with major changes and additional features. There are several API modifications and you will need to change your code to have the new (and upcoming) features. All modifications are meant to make future expansion easier.

- On the user's side, the changes are not that disruptive, but in the background there are many changes that allowed the inclusion of new features and bug fixes that would be complex to solve without the changes. Check the `documentation page <http://sysidentpy.org/notebooks.html>`__

- Many classes were basically rebuild it from scratch, so I suggest to look at the new examples of how to use the new version.

- I will present the main updates below in order to highlight features and usability and then all API changes will be reported.

- MAJOR: NARX models with Fourier basis function `Issue63 <https://github.com/wilsonrljr/sysidentpy/issues/63>`__, `Issue64 <https://github.com/wilsonrljr/sysidentpy/issues/64>`__
    - The user can choose which basis they want by importing it from sysidentpy.basis_function. Check the notebooks with examples of how to use it.
    - Polynomial and Fourier are supported for now. New basis functions will be added in next releases.

- MAJOR: NAR models `Issue58 <https://github.com/wilsonrljr/sysidentpy/issues/58>`__
    - It was already possible to build Polynomial NAR models, but with some hacks. Now the user just need to pass model_type="NAR" to build NAR models.
    - The user doesn't need to pass a vector of zeros as input anymore.
    - Works for any model structure selection algorithm (FROLS, AOLS, MetaMSS)

- Major: NFIR models `Issue59 <https://github.com/wilsonrljr/sysidentpy/issues/59>`__
    - NFIR models are models where the output depends only on past inputs. It was already possible to build Polynomial NFIR models, but with a lot of code on the user's side (much more than NAR, btw). Now the user just need to pass model_type="NFIR" to build NFIR models.
    - Works for any model structure selection algorithm (FROLS, AOLS, MetaMSS)

- Major: Select the order for the residues lags to use in Extended Least Squares - elag
    - The user can select the maximum lag of the residues to be used in the Extended Least Squares algorithm. In previous versions sysidentpy used a predefined subset of residual lags.
    - The degree of the lags follows the degree of the basis function

- Major: Residual analysis methods `Issue60 <https://github.com/wilsonrljr/sysidentpy/issues/60>`__
    - There are now specific functions to calculate the autocorrelation of the residuals and cross-correlation for the analysis of the residuals. In previous versions the calculation was limited to just two inputs, for example, limiting user usability.

- Major: Plotting methods `Issue61 <https://github.com/wilsonrljr/sysidentpy/issues/61>`__
    - The plotting functions are now separated from the models objects, so there are more flexibility regarding what to plot.
    - Residual plots were separated from the forecast plot

- API Change: sysidentpy.polynomial_basis.PolynomialNarmax is deprecated. Use sysidentpy.model_structure_selection.FROLS instead. `Issue64 <https://github.com/wilsonrljr/sysidentpy/issues/62>`__
    - Now the user doesn't need to pass the number of inputs as a parameter.
    - Added the elag parameter for unbiased_estimator. Now the user can define the number of lags of the residues for parameter estimation using the Extended Least Squares algorithm.
    - model_type parameter: now the user can select the model type to be built. The options are "NARMAX", "NAR" and "NFIR". "NARMAX" is the default. If you want to build a NAR model without any "hack", just set model_type="NAR". The same for "NFIR" models.

- API Change: sysidentpy.polynomial_basis.MetaMSS is deprecated. Use sysidentpy.model_structure_selection.MetaMSS instead. `Issue64 <https://github.com/wilsonrljr/sysidentpy/issues/64>`__
    - Now the user doesn't need to pass the number of inputs as a parameter.
    - Added the elag parameter for unbiased_estimator. Now the user can define the number of lags of the residues for parameter estimation using the Extended Least Squares algorithm.

- API Change: sysidentpy.polynomial_basis.AOLS is deprecated. Use sysidentpy.model_structure_selection.AOLS instead. `Issue64 <https://github.com/wilsonrljr/sysidentpy/issues/64>`__

- API Change: sysidentpy.polynomial_basis.SimulatePolynomialNarmax is deprecated. Use sysidentpy.simulation.SimulateNARMAX instead.

- API Change: Introducing sysidentpy.basis_function. Because NARMAX models can be built on different basis function, a new module is added to make easier to implement new basis functions in future updates `Issue64 <https://github.com/wilsonrljr/sysidentpy/issues/64>`__.
    - Each basis function class must have a fit and predict method to be used in training and prediction respectively. 

- API Change: unbiased_estimator method moved to Estimators class.
    - added elag option
    - change the build_information_matrix method to build_output_matrix

- API Change (new): sysidentpy.narmax_base
    - This is the new base for building NARMAX models. The classes have been rewritten to make it easier to expand functionality.

- API Change (new): sysidentpy.narmax_base.GenerateRegressors
    - create_narmax_code: Creates the base coding that allows representation for the NARMAX, NAR, and NFIR models.
    - regressor_space: Creates the encoding representation for the NARMAX, NAR, and NFIR models.

- API Change (new): sysidentpy.narmax_base.ModelInformation
    - _get_index_from_regressor_code: Get the index of the model code representation in regressor space.
    - _list_output_regressor_code: Create a flattened array of output regressors.
    - _list_input_regressor_code: Create a flattened array of input regressors.
    - _get_lag_from_regressor_code: Get the maximum lag from array of regressors.
    - _get_max_lag_from_model_code: the name says it all.
    - _get_max_lag: Get the maximum lag from ylag and xlag.

- API Change (new): sysidentpy.narmax_base.InformationMatrix
    - _create_lagged_X: Create a lagged matrix of inputs without combinations.
    - _create_lagged_y: Create a lagged matrix of the output without combinations.
    - build_output_matrix: Build the information matrix of output values.
    - build_input_matrix: Build the information matrix of input values.
    - build_input_output_matrix: Build the information matrix of input and output values.

- API Change (new): sysidentpy.narmax_base.ModelPrediction
    - predict: base method for prediction. Support infinity_steps ahead, one-step ahead and n-steps ahead prediction and any basis function.
    - _one_step_ahead_prediction: Perform the 1-step-ahead prediction for any basis function.
    - _n_step_ahead_prediction: Perform the n-step-ahead prediction for polynomial basis.
    - _model_prediction: Perform the infinity-step-ahead prediction for polynomial basis.
    - _narmax_predict: wrapper for NARMAX and NAR models.
    - _nfir_predict: wrapper for NFIR models.
    - _basis_function_predict: Perform the infinity-step-ahead prediction for basis functions other than polynomial.
    - basis_function_n_step_prediction: Perform the n-step-ahead prediction for basis functions other than polynomial.

- API Change (new): sysidentpy.model_structure_selection.FROLS `Issue62 <https://github.com/wilsonrljr/sysidentpy/issues/62>`__, `Issue64 <https://github.com/wilsonrljr/sysidentpy/issues/64>`__
    - Based on the old sysidentpy.polynomial_basis.PolynomialNARMAX. The class has been rebuilt with new functions and optimized code.
    - Enforcing keyword-only arguments. This is an effort to promote clear and non-ambiguous use of the library.
    - Add support for new basis functions.
    - The user can choose the residual lags.
    - No need to pass the number of inputs anymore.
    - Improved docstring.
    - Fixed minor grammatical and spelling mistakes.
    - New prediction method.
    - many under the hood changes.
 
- API Change (new): sysidentpy.model_structure_selection.MetaMSS `Issue64 <https://github.com/wilsonrljr/sysidentpy/issues/64>`__
    - Based on the old sysidentpy.polynomial_basis.MetaMSS. The class has been rebuilt with new functions and optimized code.
    - Enforcing keyword-only arguments. This is an effort to promote clear and non-ambiguous use of the library.
    - The user can choose the residual lags.
    - Extended Least Squares support.
    - Add support for new basis functions.
    - No need to pass the number of inputs anymore.
    - Improved docstring.
    - Fixed minor grammatical and spelling mistakes.
    - New prediction method.
    - many under the hood changes.

- API Change (new): sysidentpy.model_structure_selection.AOLS `Issue64 <https://github.com/wilsonrljr/sysidentpy/issues/64>`__
    - Based on the old sysidentpy.polynomial_basis.AOLS. The class has been rebuilt with new functions and optimized code.
    - Enforcing keyword-only arguments. This is an effort to promote clear and non-ambiguous use of the library.
    - Add support for new basis functions.
    - No need to pass the number of inputs anymore.
    - Improved docstring.
    - Change "l" parameter to "L".
    - Fixed minor grammatical and spelling mistakes.
    - New prediction method.
    - many under the hood changes.

- API Change (new): sysidentpy.simulation.SimulateNARMAX
    - Based on the old sysidentpy.polynomial_basis.SimulatePolynomialNarmax. The class has been rebuilt with new functions and optimized code.
    - Fix the Extended Least Squares support.
    - Fix n-steps ahead prediction and 1-step ahead prediction.
    - Enforcing keyword-only arguments. This is an effort to promote clear and non-ambiguous use of the library.
    - The user can choose the residual lags.
    - Improved docstring.
    - Fixed minor grammatical and spelling mistakes.
    - New prediction method.
    - Do not inherit from the structure selection algorithm anymore, only from narmax_base. Avoid circular import and other issues.
    - many under the hood changes.

- API Change (new): sysidentpy.residues
    - compute_residues_autocorrelation: the name says it all.
    - calculate_residues: get the residues from y and yhat.
    - get_unnormalized_e_acf: compute the unnormalized autocorrelation of the residues.
    - compute_cross_correlation: compute cross correlation between two signals.
    - _input_ccf
    - _normalized_correlation: compute the normalized correlation between two signals.

- API Change (new): sysidentpy.utils.plotting
    - plot_results: plot the forecast
    - plot_residues_correlation: the name says it all.

- API Change (new): sysidentpy.utils.display_results
    - results: return the model regressors, estimated parameter and ERR index of the fitted model in a table.

- DOC: Air passenger benchmark
    - Added notebook with Air passenger forecasting benchmark.
    - We compare SysIdentPy against prophet, neuralprophet, autoarima, tbats and many more.

- DOC: Load forecasting benchmark
    - Added notebook with load forecasting benchmark.

- DOC: PV forecasting benchmark
    - Added notebook with PV forecasting benchmark.

- DOC: Presenting main functionality
    - Example rewritten following the new api.
    - Fixed minor grammatical and spelling mistakes.

- DOC: Multiple Inputs usage
    - Example rewritten following the new api
    - Fixed minor grammatical and spelling mistakes.

- DOC: Information Criteria - Examples
    - Example rewritten following the new api.
    - Fixed minor grammatical and spelling mistakes.

- DOC: Important notes and examples of how to use Extended Least Squares
    - Example rewritten following the new api.
    - Fixed minor grammatical and spelling mistakes.

- DOC: Setting specific lags
    - Example rewritten following the new api.
    - Fixed minor grammatical and spelling mistakes.

- DOC: Parameter Estimation
    - Example rewritten following the new api.
    - Fixed minor grammatical and spelling mistakes.

- DOC: Using the Meta-Model Structure Selection (MetaMSS) algorithm for building Polynomial NARX models
    - Example rewritten following the new api.
    - Fixed minor grammatical and spelling mistakes.

- DOC: Using the Accelerated Orthogonal Least-Squares algorithm for building Polynomial NARX models
    - Example rewritten following the new api.
    - Fixed minor grammatical and spelling mistakes.

- DOC: Example: F-16 Ground Vibration Test benchmark
    - Example rewritten following the new api.
    - Fixed minor grammatical and spelling mistakes.

- DOC: Building NARX Neural Network using Sysidentpy
    - Example rewritten following the new api.
    - Fixed minor grammatical and spelling mistakes.

- DOC: Building NARX models using general estimators
    - Example rewritten following the new api.
    - Fixed minor grammatical and spelling mistakes.

- DOC: Simulate a Predefined Model
    - Example rewritten following the new api.
    - Fixed minor grammatical and spelling mistakes.

- DOC: System Identification Using Adaptive Filters
    - Example rewritten following the new api.
    - Fixed minor grammatical and spelling mistakes.

- DOC: Identification of an electromechanical system
    - Example rewritten following the new api.
    - Fixed minor grammatical and spelling mistakes.

- DOC: Example: N-steps-ahead prediction - F-16 Ground Vibration Test benchmark
    - Example rewritten following the new api.
    - Fixed minor grammatical and spelling mistakes.

- DOC: Introduction to NARMAX models
    - Fixed grammatical and spelling mistakes.



v0.1.6
------

CONTRIBUTORS
~~~~~~~~~~~~

- wilsonrljr

CHANGES
~~~~~~~

- MAJOR: Meta-Model Structure Selection Algorithm (Meta-MSS).
    - A new method for build NARMAX models based on metaheuristics. The algorithm uses a Binary hybrid Particle Swarm Optimization and Gravitational Search Algorithm with a new cost function to build parsimonious models.
    
    - New class for the BPSOGSA algorithm. New algorithms can be adapted in the Meta-MSS framework.
	
    - Future updates will add NARX models for classification and multiobjective model structure selection.

- MAJOR: Accelerated Orthogonal Least-Squares algorithm.
    - Added the new class AOLS to build NARX models using the Accelerated Orthogonal Least-Squares algorithm.
    
    - At the best of my knowledge, this is the first time this algorithm is used in the NARMAX framework. The tests I've made are promising, but use it with caution until the results are formalized into a research paper.

- Added notebook with a simple example of how to use MetaMSS and a simple model comparison of the Electromechanical system.

- Added notebook with a simple example of how to use AOLS

- Added ModelInformation class. This class have methods to return model information such as max_lag of a model code.
    - added _list_output_regressor_code
    - added _list_input_regressor_code
    - added _get_lag_from_regressor_code
    - added _get_max_lag_from_model_code

- Minor performance improvement: added the argument "predefined_regressors" in build_information_matrix function on base.py
    to improve the performance of the Simulation method.

- Pytorch is now an optional dependency. Use pip install sysidentpy['full'] 

- Fix code format issues.

- Fixed minor grammatical and spelling mistakes.

- Fix issues related to html on Jupyter notebooks examples on documentation.

- Updated Readme with examples of how to use.

- Improved descriptions and comments in methods.

- metaheuristics.bpsogsa (detailed description on code docstring)
    - added evaluate_objective_function
    - added optimize
    - added generate_random_population
    - added mass_calculation
    - added calculate_gravitational_constant
    - added calculate_acceleration
    - added update_velocity_position

- FIX issue #52


v0.1.5
------

CONTRIBUTORS
~~~~~~~~~~~~

- wilsonrljr

CHANGES
~~~~~~~

- MAJOR: n-steps-ahead prediction.
    - Now you can define the numbers of steps ahead in the predict function.
	- Only for Polynomial models for now. Next update will bring this functionality to Neural NARX and General Estimators.

- MAJOR: Simulating predefined models.
    - Added the new class SimulatePolynomialNarmax to handle the simulation of known model structures.
    - Now you can simulate predefined models by just passing the model structure codification. Check the notebook examples.

- Added 4 new notebooks in the example section.

- Added iterative notebooks. Now you can run the notebooks in Jupyter notebook section of the documentation in Colab.

- Fix code format issues.

- Added new tests for SimulatePolynomialNarmax and generate_data.

- Started changes related to numpy 1.19.4 update. There are still some Deprecation warnings that will be fixed in next update.

- Fix issues related to html on Jupyter notebooks examples on documentation.

- Updated Readme with examples of how to use.



v0.1.4
------

CONTRIBUTORS
~~~~~~~~~~~~

- wilsonrljr

CHANGES
~~~~~~~

- MAJOR: Introducing NARX Neural Network in SysIdentPy.
    - Now you can build NARX Neural Network on SysIdentPy.
    - This feature is built on top of Pytorch. See the docs for more details and examples of how to use.

- MAJOR: Introducing general estimators in SysIdentPy.
    - Now you are able to use any estimator that have Fit/Predict methods (estimators from Sklearn and Catboost, for example) and build NARX models based on those estimators.
    - We use the core functions of SysIdentPy and keep the Fit/Predict approach from those estimators to keep the process easy to use.
    - More estimators are coming soon like XGboost.

- Added notebooks to show how to build NARX neural Network.

- Added notebooks to show how to build NARX models using general estimators.

- Changed the default parameters of the plot_results function.

- NOTE: We will keeping improving the Polynomial NARX models (new model structure selection algorithms and multiobjective identification
is on our roadmap). These recent modifications will allow us to introduce new NARX models like PWARX models very soon.

- New template for the documentation site.

- Fix issues related to html on Jupyter notebooks examples on documentation.

- Updated Readme with examples of how to use.


v0.1.3
------

CONTRIBUTORS
~~~~~~~~~~~~

- wilsonrljr
- renard162

CHANGES
~~~~~~~

- Fixed a bug concerning the xlag and ylag in multiple input scenarios.
- Refactored predict function. Improved performance up to 87% depending on the number of regressors.
- You can set lags with different size for each input.
- Added a new function to get the max value of xlag and ylag. Work with int, list, nested lists.
- Fixed tests for information criteria.
- Added SysIdentPy logo.
- Refactored code of all classes following PEP 8 guidelines to improve readability.
- Added Citation information on Readme.
- Changes on information Criteria tests.
- Added workflow to run the tests when merge branch into master.
- Added new site domain.
- Updated docs.

Changes in SysIdentPy
=====================

v0.1.7
------

CONTRIBUTORS
~~~~~~~~~~~~

- wilsonrljr

CHANGES
~~~~~~~
- The update **v0.1.7**  has been released with major changes and additional features. There are several API modifications and you will need to change your code to have the new (and upcoming) features. All modifications are meant to make future expansion easier.

- On the user's side, the changes are not that disruptive, but in the background there are many changes that allowed the inclusion of new features and bug fixes that would be complex to solve without the changes. Check the `documentation page <http://sysidentpy.org/notebooks.html>`__

- Many classes were basically rebuild it from scratch, so I suggest to look at the new examples of how to use the new version.

- I will present the main updates below in order to highlight features and usability and then all API changes will be reported.

- MAJOR: NARX models with Fourier basis function `Issue63 <https://github.com/wilsonrljr/sysidentpy/issues/63>`__, `Issue64 <https://github.com/wilsonrljr/sysidentpy/issues/64>`__
    - The user can choose which basis they want by importing it from sysidentpy.basis_function. Check the notebooks with examples of how to use it.
    - Polynomial and Fourier are supported for now. New basis functions will be added in next releases.

- MAJOR: NAR models `Issue58 <https://github.com/wilsonrljr/sysidentpy/issues/58>`__
    - It was already possible to build Polynomial NAR models, but with some hacks. Now the user just need to pass model_type="NAR" to build NAR models.
    - The user doesn't need to pass a vector of zeros as input anymore.
    - Works for any model structure selection algorithm (FROLS, AOLS, MetaMSS)

- Major: NFIR models `Issue59 <https://github.com/wilsonrljr/sysidentpy/issues/59>`__
    - NFIR models are models where the output depends only on past inputs. It was already possible to build Polynomial NFIR models, but with a lot of code on the user's side (much more than NAR, btw). Now the user just need to pass model_type="NFIR" to build NFIR models.
    - Works for any model structure selection algorithm (FROLS, AOLS, MetaMSS)

- Major: Select the order for the residues lags to use in Extended Least Squares - elag
    - The user can select the maximum lag of the residues to be used in the Extended Least Squares algorithm. In previous versions sysidentpy used a predefined subset of residual lags.
    - The degree of the lags follows the degree of the basis function

- Major: Residual analysis methods `Issue60 <https://github.com/wilsonrljr/sysidentpy/issues/60>`__
    - There are now specific functions to calculate the autocorrelation of the residuals and cross-correlation for the analysis of the residuals. In previous versions the calculation was limited to just two inputs, for example, limiting user usability.

- Major: Plotting methods `Issue61 <https://github.com/wilsonrljr/sysidentpy/issues/61>`__
    - The plotting functions are now separated from the models objects, so there are more flexibility regarding what to plot.
    - Residual plots were separated from the forecast plot

- API Change: sysidentpy.polynomial_basis.PolynomialNarmax is deprecated. Use sysidentpy.model_structure_selection.FROLS instead. `Issue64 <https://github.com/wilsonrljr/sysidentpy/issues/62>`__
    - Now the user doesn't need to pass the number of inputs as a parameter.
    - Added the elag parameter for unbiased_estimator. Now the user can define the number of lags of the residues for parameter estimation using the Extended Least Squares algorithm.
    - model_type parameter: now the user can select the model type to be built. The options are "NARMAX", "NAR" and "NFIR". "NARMAX" is the default. If you want to build a NAR model without any "hack", just set model_type="NAR". The same for "NFIR" models.

- API Change: sysidentpy.polynomial_basis.MetaMSS is deprecated. Use sysidentpy.model_structure_selection.MetaMSS instead. `Issue64 <https://github.com/wilsonrljr/sysidentpy/issues/64>`__
    - Now the user doesn't need to pass the number of inputs as a parameter.
    - Added the elag parameter for unbiased_estimator. Now the user can define the number of lags of the residues for parameter estimation using the Extended Least Squares algorithm.

- API Change: sysidentpy.polynomial_basis.AOLS is deprecated. Use sysidentpy.model_structure_selection.AOLS instead. `Issue64 <https://github.com/wilsonrljr/sysidentpy/issues/64>`__

- API Change: sysidentpy.polynomial_basis.SimulatePolynomialNarmax is deprecated. Use sysidentpy.simulation.SimulateNARMAX instead.

- API Change: Introducing sysidentpy.basis_function. Because NARMAX models can be built on different basis function, a new module is added to make easier to implement new basis functions in future updates `Issue64 <https://github.com/wilsonrljr/sysidentpy/issues/64>`__.
    - Each basis function class must have a fit and predict method to be used in training and prediction respectively. 

- API Change: unbiased_estimator method moved to Estimators class.
    - added elag option
    - change the build_information_matrix method to build_output_matrix

- API Change (new): sysidentpy.narmax_base
    - This is the new base for building NARMAX models. The classes have been rewritten to make it easier to expand functionality.

- API Change (new): sysidentpy.narmax_base.GenerateRegressors
    - create_narmax_code: Creates the base coding that allows representation for the NARMAX, NAR, and NFIR models.
    - regressor_space: Creates the encoding representation for the NARMAX, NAR, and NFIR models.

- API Change (new): sysidentpy.narmax_base.ModelInformation
    - _get_index_from_regressor_code: Get the index of the model code representation in regressor space.
    - _list_output_regressor_code: Create a flattened array of output regressors.
    - _list_input_regressor_code: Create a flattened array of input regressors.
    - _get_lag_from_regressor_code: Get the maximum lag from array of regressors.
    - _get_max_lag_from_model_code: the name says it all.
    - _get_max_lag: Get the maximum lag from ylag and xlag.

- API Change (new): sysidentpy.narmax_base.InformationMatrix
    - _create_lagged_X: Create a lagged matrix of inputs without combinations.
    - _create_lagged_y: Create a lagged matrix of the output without combinations.
    - build_output_matrix: Build the information matrix of output values.
    - build_input_matrix: Build the information matrix of input values.
    - build_input_output_matrix: Build the information matrix of input and output values.

- API Change (new): sysidentpy.narmax_base.ModelPrediction
    - predict: base method for prediction. Support infinity_steps ahead, one-step ahead and n-steps ahead prediction and any basis function.
    - _one_step_ahead_prediction: Perform the 1-step-ahead prediction for any basis function.
    - _n_step_ahead_prediction: Perform the n-step-ahead prediction for polynomial basis.
    - _model_prediction: Perform the infinity-step-ahead prediction for polynomial basis.
    - _narmax_predict: wrapper for NARMAX and NAR models.
    - _nfir_predict: wrapper for NFIR models.
    - _basis_function_predict: Perform the infinity-step-ahead prediction for basis functions other than polynomial.
    - basis_function_n_step_prediction: Perform the n-step-ahead prediction for basis functions other than polynomial.

- API Change (new): sysidentpy.model_structure_selection.FROLS `Issue62 <https://github.com/wilsonrljr/sysidentpy/issues/62>`__, `Issue64 <https://github.com/wilsonrljr/sysidentpy/issues/64>`__
    - Based on the old sysidentpy.polynomial_basis.PolynomialNARMAX. The class has been rebuilt with new functions and optimized code.
    - Enforcing keyword-only arguments. This is an effort to promote clear and non-ambiguous use of the library.
    - Add support for new basis functions.
    - The user can choose the residual lags.
    - No need to pass the number of inputs anymore.
    - Improved docstring.
    - Fixed minor grammatical and spelling mistakes.
    - New prediction method.
    - many under the hood changes.
 
- API Change (new): sysidentpy.model_structure_selection.MetaMSS `Issue64 <https://github.com/wilsonrljr/sysidentpy/issues/64>`__
    - Based on the old sysidentpy.polynomial_basis.MetaMSS. The class has been rebuilt with new functions and optimized code.
    - Enforcing keyword-only arguments. This is an effort to promote clear and non-ambiguous use of the library.
    - The user can choose the residual lags.
    - Extended Least Squares support.
    - Add support for new basis functions.
    - No need to pass the number of inputs anymore.
    - Improved docstring.
    - Fixed minor grammatical and spelling mistakes.
    - New prediction method.
    - many under the hood changes.

- API Change (new): sysidentpy.model_structure_selection.AOLS `Issue64 <https://github.com/wilsonrljr/sysidentpy/issues/64>`__
    - Based on the old sysidentpy.polynomial_basis.AOLS. The class has been rebuilt with new functions and optimized code.
    - Enforcing keyword-only arguments. This is an effort to promote clear and non-ambiguous use of the library.
    - Add support for new basis functions.
    - No need to pass the number of inputs anymore.
    - Improved docstring.
    - Change "l" parameter to "L".
    - Fixed minor grammatical and spelling mistakes.
    - New prediction method.
    - many under the hood changes.

- API Change (new): sysidentpy.simulation.SimulateNARMAX
    - Based on the old sysidentpy.polynomial_basis.SimulatePolynomialNarmax. The class has been rebuilt with new functions and optimized code.
    - Fix the Extended Least Squares support.
    - Fix n-steps ahead prediction and 1-step ahead prediction.
    - Enforcing keyword-only arguments. This is an effort to promote clear and non-ambiguous use of the library.
    - The user can choose the residual lags.
    - Improved docstring.
    - Fixed minor grammatical and spelling mistakes.
    - New prediction method.
    - Do not inherit from the structure selection algorithm anymore, only from narmax_base. Avoid circular import and other issues.
    - many under the hood changes.

- API Change (new): sysidentpy.residues
    - compute_residues_autocorrelation: the name says it all.
    - calculate_residues: get the residues from y and yhat.
    - get_unnormalized_e_acf: compute the unnormalized autocorrelation of the residues.
    - compute_cross_correlation: compute cross correlation between two signals.
    - _input_ccf
    - _normalized_correlation: compute the normalized correlation between two signals.

- API Change (new): sysidentpy.utils.plotting
    - plot_results: plot the forecast
    - plot_residues_correlation: the name says it all.

- API Change (new): sysidentpy.utils.display_results
    - results: return the model regressors, estimated parameter and ERR index of the fitted model in a table.

- DOC: Air passenger benchmark
    - Added notebook with Air passenger forecasting benchmark.
    - We compare SysIdentPy against prophet, neuralprophet, autoarima, tbats and many more.

- DOC: Load forecasting benchmark
    - Added notebook with load forecasting benchmark.

- DOC: PV forecasting benchmark
    - Added notebook with PV forecasting benchmark.

- DOC: Presenting main functionality
    - Example rewritten following the new api.
    - Fixed minor grammatical and spelling mistakes.

- DOC: Multiple Inputs usage
    - Example rewritten following the new api
    - Fixed minor grammatical and spelling mistakes.

- DOC: Information Criteria - Examples
    - Example rewritten following the new api.
    - Fixed minor grammatical and spelling mistakes.

- DOC: Important notes and examples of how to use Extended Least Squares
    - Example rewritten following the new api.
    - Fixed minor grammatical and spelling mistakes.

- DOC: Setting specific lags
    - Example rewritten following the new api.
    - Fixed minor grammatical and spelling mistakes.

- DOC: Parameter Estimation
    - Example rewritten following the new api.
    - Fixed minor grammatical and spelling mistakes.

- DOC: Using the Meta-Model Structure Selection (MetaMSS) algorithm for building Polynomial NARX models
    - Example rewritten following the new api.
    - Fixed minor grammatical and spelling mistakes.

- DOC: Using the Accelerated Orthogonal Least-Squares algorithm for building Polynomial NARX models
    - Example rewritten following the new api.
    - Fixed minor grammatical and spelling mistakes.

- DOC: Example: F-16 Ground Vibration Test benchmark
    - Example rewritten following the new api.
    - Fixed minor grammatical and spelling mistakes.

- DOC: Building NARX Neural Network using Sysidentpy
    - Example rewritten following the new api.
    - Fixed minor grammatical and spelling mistakes.

- DOC: Building NARX models using general estimators
    - Example rewritten following the new api.
    - Fixed minor grammatical and spelling mistakes.

- DOC: Simulate a Predefined Model
    - Example rewritten following the new api.
    - Fixed minor grammatical and spelling mistakes.

- DOC: System Identification Using Adaptive Filters
    - Example rewritten following the new api.
    - Fixed minor grammatical and spelling mistakes.

- DOC: Identification of an electromechanical system
    - Example rewritten following the new api.
    - Fixed minor grammatical and spelling mistakes.

- DOC: Example: N-steps-ahead prediction - F-16 Ground Vibration Test benchmark
    - Example rewritten following the new api.
    - Fixed minor grammatical and spelling mistakes.

- DOC: Introduction to NARMAX models
    - Fixed grammatical and spelling mistakes.



v0.1.6
------

CONTRIBUTORS
~~~~~~~~~~~~

- wilsonrljr

CHANGES
~~~~~~~

- MAJOR: Meta-Model Structure Selection Algorithm (Meta-MSS).
    - A new method for build NARMAX models based on metaheuristics. The algorithm uses a Binary hybrid Particle Swarm Optimization and Gravitational Search Algorithm with a new cost function to build parsimonious models.
    
    - New class for the BPSOGSA algorithm. New algorithms can be adapted in the Meta-MSS framework.
	
    - Future updates will add NARX models for classification and multiobjective model structure selection.

- MAJOR: Accelerated Orthogonal Least-Squares algorithm.
    - Added the new class AOLS to build NARX models using the Accelerated Orthogonal Least-Squares algorithm.
    
    - At the best of my knowledge, this is the first time this algorithm is used in the NARMAX framework. The tests I've made are promising, but use it with caution until the results are formalized into a research paper.

- Added notebook with a simple example of how to use MetaMSS and a simple model comparison of the Electromechanical system.

- Added notebook with a simple example of how to use AOLS

- Added ModelInformation class. This class have methods to return model information such as max_lag of a model code.
    - added _list_output_regressor_code
    - added _list_input_regressor_code
    - added _get_lag_from_regressor_code
    - added _get_max_lag_from_model_code

- Minor performance improvement: added the argument "predefined_regressors" in build_information_matrix function on base.py
    to improve the performance of the Simulation method.

- Pytorch is now an optional dependency. Use pip install sysidentpy['full'] 

- Fix code format issues.

- Fixed minor grammatical and spelling mistakes.

- Fix issues related to html on Jupyter notebooks examples on documentation.

- Updated Readme with examples of how to use.

- Improved descriptions and comments in methods.

- metaheuristics.bpsogsa (detailed description on code docstring)
    - added evaluate_objective_function
    - added optimize
    - added generate_random_population
    - added mass_calculation
    - added calculate_gravitational_constant
    - added calculate_acceleration
    - added update_velocity_position

- FIX issue #52


v0.1.5
------

CONTRIBUTORS
~~~~~~~~~~~~

- wilsonrljr

CHANGES
~~~~~~~

- MAJOR: n-steps-ahead prediction.
    - Now you can define the numbers of steps ahead in the predict function.
	- Only for Polynomial models for now. Next update will bring this functionality to Neural NARX and General Estimators.

- MAJOR: Simulating predefined models.
    - Added the new class SimulatePolynomialNarmax to handle the simulation of known model structures.
    - Now you can simulate predefined models by just passing the model structure codification. Check the notebook examples.

- Added 4 new notebooks in the example section.

- Added iterative notebooks. Now you can run the notebooks in Jupyter notebook section of the documentation in Colab.

- Fix code format issues.

- Added new tests for SimulatePolynomialNarmax and generate_data.

- Started changes related to numpy 1.19.4 update. There are still some Deprecation warnings that will be fixed in next update.

- Fix issues related to html on Jupyter notebooks examples on documentation.

- Updated Readme with examples of how to use.



v0.1.4
------

CONTRIBUTORS
~~~~~~~~~~~~

- wilsonrljr

CHANGES
~~~~~~~

- MAJOR: Introducing NARX Neural Network in SysIdentPy.
    - Now you can build NARX Neural Network on SysIdentPy.
    - This feature is built on top of Pytorch. See the docs for more details and examples of how to use.

- MAJOR: Introducing general estimators in SysIdentPy.
    - Now you are able to use any estimator that have Fit/Predict methods (estimators from Sklearn and Catboost, for example) and build NARX models based on those estimators.
    - We use the core functions of SysIdentPy and keep the Fit/Predict approach from those estimators to keep the process easy to use.
    - More estimators are coming soon like XGboost.

- Added notebooks to show how to build NARX neural Network.

- Added notebooks to show how to build NARX models using general estimators.

- Changed the default parameters of the plot_results function.

- NOTE: We will keeping improving the Polynomial NARX models (new model structure selection algorithms and multiobjective identification
is on our roadmap). These recent modifications will allow us to introduce new NARX models like PWARX models very soon.

- New template for the documentation site.

- Fix issues related to html on Jupyter notebooks examples on documentation.

- Updated Readme with examples of how to use.


v0.1.3
------

CONTRIBUTORS
~~~~~~~~~~~~

- wilsonrljr
- renard162

CHANGES
~~~~~~~

- Fixed a bug concerning the xlag and ylag in multiple input scenarios.
- Refactored predict function. Improved performance up to 87% depending on the number of regressors.
- You can set lags with different size for each input.
- Added a new function to get the max value of xlag and ylag. Work with int, list, nested lists.
- Fixed tests for information criteria.
- Added SysIdentPy logo.
- Refactored code of all classes following PEP 8 guidelines to improve readability.
- Added Citation information on Readme.
- Changes on information Criteria tests.
- Added workflow to run the tests when merge branch into master.
- Added new site domain.
- Updated docs.

Changes in SysIdentPy
=====================

v0.1.6
------

CONTRIBUTORS
~~~~~~~~~~~~

- wilsonrljr

CHANGES
~~~~~~~

- MAJOR: Meta-Model Structure Selection Algorithm (Meta-MSS).
    - A new method for build NARMAX models based on metaheuristics. The algorithm uses a Binary hybrid Particle Swarm Optimization and Gravitational Search Algorithm with a new cost function to build parsimonious models.
    
    - New class for the BPSOGSA algorithm. New algorithms can be adapted in the Meta-MSS framework.
	
    - Future updates will add NARX models for classification and multiobjective model structure selection.

- MAJOR: Accelerated Orthogonal Least-Squares algorithm.
    - Added the new class AOLS to build NARX models using the Accelerated Orthogonal Least-Squares algorithm.
    
    - At the best of my knowledge, this is the first time this algorithm is used in the NARMAX framework. The tests I've made are promising, but use it with caution until the results are formalized into a research paper.

- Added notebook with a simple example of how to use MetaMSS and a simple model comparison of the Electromechanical system.

- Added notebook with a simple example of how to use AOLS

- Added ModelInformation class. This class have methods to return model information such as max_lag of a model code.
    - added _list_output_regressor_code
    - added _list_input_regressor_code
    - added _get_lag_from_regressor_code
    - added _get_max_lag_from_model_code

- Minor performance improvement: added the argument "predefined_regressors" in build_information_matrix function on base.py
    to improve the performance of the Simulation method.

- Pytorch is now an optional dependency. Use pip install sysidentpy['full'] 

- Fix code format issues.

- Fixed minor grammatical and spelling mistakes.

- Fix issues related to html on Jupyter notebooks examples on documentation.

- Updated Readme with examples of how to use.

- Improved descriptions and comments in methods.

- metaheuristics.bpsogsa (detailed description on code docstring)
    - added evaluate_objective_function
    - added optimize
    - added generate_random_population
    - added mass_calculation
    - added calculate_gravitational_constant
    - added calculate_acceleration
    - added update_velocity_position

- FIX issue #52


v0.1.5
------

CONTRIBUTORS
~~~~~~~~~~~~

- wilsonrljr

CHANGES
~~~~~~~

- MAJOR: n-steps-ahead prediction.
    - Now you can define the numbers of steps ahead in the predict function.
	- Only for Polynomial models for now. Next update will bring this functionality to Neural NARX and General Estimators.

- MAJOR: Simulating predefined models.
    - Added the new class SimulatePolynomialNarmax to handle the simulation of known model structures.
    - Now you can simulate predefined models by just passing the model structure codification. Check the notebook examples.

- Added 4 new notebooks in the example section.

- Added iterative notebooks. Now you can run the notebooks in Jupyter notebook section of the documentation in Colab.

- Fix code format issues.

- Added new tests for SimulatePolynomialNarmax and generate_data.

- Started changes related to numpy 1.19.4 update. There are still some Deprecation warnings that will be fixed in next update.

- Fix issues related to html on Jupyter notebooks examples on documentation.

- Updated Readme with examples of how to use.



v0.1.4
------

CONTRIBUTORS
~~~~~~~~~~~~

- wilsonrljr

CHANGES
~~~~~~~

- MAJOR: Introducing NARX Neural Network in SysIdentPy.
    - Now you can build NARX Neural Network on SysIdentPy.
    - This feature is built on top of Pytorch. See the docs for more details and examples of how to use.

- MAJOR: Introducing general estimators in SysIdentPy.
    - Now you are able to use any estimator that have Fit/Predict methods (estimators from Sklearn and Catboost, for example) and build NARX models based on those estimators.
    - We use the core functions of SysIdentPy and keep the Fit/Predict approach from those estimators to keep the process easy to use.
    - More estimators are coming soon like XGboost.

- Added notebooks to show how to build NARX neural Network.

- Added notebooks to show how to build NARX models using general estimators.

- Changed the default parameters of the plot_results function.

- NOTE: We will keeping improving the Polynomial NARX models (new model structure selection algorithms and multiobjective identification
is on our roadmap). These recent modifications will allow us to introduce new NARX models like PWARX models very soon.

- New template for the documentation site.

- Fix issues related to html on Jupyter notebooks examples on documentation.

- Updated Readme with examples of how to use.


v0.1.3
------

CONTRIBUTORS
~~~~~~~~~~~~

- wilsonrljr
- renard162

CHANGES
~~~~~~~

- Fixed a bug concerning the xlag and ylag in multiple input scenarios.
- Refactored predict function. Improved performance up to 87% depending on the number of regressors.
- You can set lags with different size for each input.
- Added a new function to get the max value of xlag and ylag. Work with int, list, nested lists.
- Fixed tests for information criteria.
- Added SysIdentPy logo.
- Refactored code of all classes following PEP 8 guidelines to improve readability.
- Added Citation information on Readme.
- Changes on information Criteria tests.
- Added workflow to run the tests when merge branch into master.
- Added new site domain.
- Updated docs.

Changes in SysIdentPy
=====================

v0.1.5
------

CONTRIBUTORS
~~~~~~~~~~~~

- wilsonrljr

CHANGES
~~~~~~~

- MAJOR: n-steps-ahead prediction.
    - Now you can define the numbers of steps ahead in the predict function.
	- Only for Polynomial models for now. Next update will bring this functionality to Neural NARX and General Estimators.

- MAJOR: Simulating predefined models.
    - Added the new class SimulatePolynomialNarmax to handle the simulation of known model structures.
    - Now you can simulate predefined models by just passing the model structure codification. Check the notebook examples.

- Added 4 new notebooks in the example section.

- Added iterative notebooks. Now you can run the notebooks in Jupyter notebook section of the documentation in Colab.

- Fix code format issues.

- Added new tests for SimulatePolynomialNarmax and generate_data.

- Started changes related to numpy 1.19.4 update. There are still some Deprecation warnings that will be fixed in next update.

- Fix issues related to html on Jupyter notebooks examples on documentation.

- Updated Readme with examples of how to use.


v0.1.4
------

CONTRIBUTORS
~~~~~~~~~~~~

- wilsonrljr

CHANGES
~~~~~~~

- MAJOR: Introducing NARX Neural Network in SysIdentPy.
    - Now you can build NARX Neural Network on SysIdentPy.
    - This feature is built on top of Pytorch. See the docs for more details and examples of how to use.

- MAJOR: Introducing general estimators in SysIdentPy.
    - Now you are able to use any estimator that have Fit/Predict methods (estimators from Sklearn and Catboost, for example) and build NARX models based on those estimators.
    - We use the core functions of SysIdentPy and keep the Fit/Predict approach from those estimators to keep the process easy to use.
    - More estimators are coming soon like XGboost.

- Added notebooks to show how to build NARX neural Network.

- Added notebooks to show how to build NARX models using general estimators.

- Changed the default parameters of the plot_results function.

- NOTE: We will keeping improving the Polynomial NARX models (new model structure selection algorithms and multiobjective identification
is on our roadmap). These recent modifications will allow us to introduce new NARX models like PWARX models very soon.

- New template for the documentation site.

- Fix issues related to html on Jupyter notebooks examples on documentation.

- Updated Readme with examples of how to use.


v0.1.3
------

CONTRIBUTORS
~~~~~~~~~~~~

- wilsonrljr
- renard162

CHANGES
~~~~~~~

- Fixed a bug concerning the xlag and ylag in multiple input scenarios.
- Refactored predict function. Improved performance up to 87% depending on the number of regressors.
- You can set lags with different size for each input.
- Added a new function to get the max value of xlag and ylag. Work with int, list, nested lists.
- Fixed tests for information criteria.
- Added SysIdentPy logo.
- Refactored code of all classes following PEP 8 guidelines to improve readability.
- Added Citation information on Readme.
- Changes on information Criteria tests.
- Added workflow to run the tests when merge branch into master.
- Added new site domain.
- Updated docs.


v0.1.5
------
CONTRIBUTORS
CHANGES
Changes in SysIdentPy
=====================

v0.1.4
------

CONTRIBUTORS
~~~~~~~~~~~~

- wilsonrljr

CHANGES
~~~~~~~

- MAJOR: Introducing NARX Neural Network in SysIdentPy.
    - Now you can build NARX Neural Network on SysIdentPy.
    - This feature is built on top of Pytorch. See the docs for more details and examples of how to use.

- MAJOR: Introducing general estimators in SysIdentPy.
    - Now you are able to use any estimator that have Fit/Predict methods (estimators from Sklearn and Catboost, for example) and build NARX models based on those estimators.
    - We use the core functions of SysIdentPy and keep the Fit/Predict approach from those estimators to keep the process easy to use.
    - More estimators are coming soon like XGboost.

- Added notebooks to show how to build NARX neural Network.

- Added notebooks to show how to build NARX models using general estimators.

- Changed the default parameters of the plot_results function.

- NOTE: We will keeping improving the Polynomial NARX models (new model structure selection algorithms and multiobjective identification
is on our roadmap). These recent modifications will allow us to introduce new NARX models like PWARX models very soon.

- New template for the documentation site.

- Fix issues related to html on Jupyter notebooks examples on documentation.

- Updated Readme with examples of how to use.


v0.1.3
------

CONTRIBUTORS
~~~~~~~~~~~~

- wilsonrljr
- renard162

CHANGES
~~~~~~~

- Fixed a bug concerning the xlag and ylag in multiple input scenarios.
- Refactored predict function. Improved performance up to 87% depending on the number of regressors.
- You can set lags with different size for each input.
- Added a new function to get the max value of xlag and ylag. Work with int, list, nested lists.
- Fixed tests for information criteria.
- Added SysIdentPy logo.
- Refactored code of all classes following PEP 8 guidelines to improve readability.
- Added Citation information on Readme.
- Changes on information Criteria tests.
- Added workflow to run the tests when merge branch into master.
- Added new site domain.
- Updated docs.


v0.1.5
------
CONTRIBUTORS
CHANGES
