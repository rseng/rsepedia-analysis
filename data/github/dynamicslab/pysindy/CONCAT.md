---
name: Bug report
about: Create a report to help us improve
title: "[BUG]"
labels: ''
assignees: ''

---

<!-- Please describe the issue in detail here, and fill in the fields below -->

### Reproducing code example:

<!-- A short code example that reproduces the problem/missing feature. It should be
self-contained, i.e., possible to run as-is via 'python myproblem.py' -->

```python
import pysindy
<< your code here >>
```

### Error message:

<!-- Full error message, if any (starting from line Traceback: ...) -->

### PySINDy/Python version information:

<!-- Output from 'import sys, pysindy; print(pysindy.__version__, sys.version)' -->
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: ''
assignees: ''

---

## Is your feature request related to a problem? Please describe.
<!-- A clear and concise description of what the problem is. Ex. I'm always frustrated when [...] -->

## Describe the solution you'd like
<!-- A clear and concise description of what you want to happen. -->


## Describe alternatives you've considered
<!-- A clear and concise description of any alternative solutions or features you've considered. -->


## Additional context
<!-- Add any other context or screenshots about the feature request here. -->

# Contributor Covenant Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, caste, color, religion, or sexual
identity and orientation.

We pledge to act and interact in ways that contribute to an open, welcoming,
diverse, inclusive, and healthy community.

## Our Standards

Examples of behavior that contributes to a positive environment for our
community include:

* Demonstrating empathy and kindness toward other people
* Being respectful of differing opinions, viewpoints, and experiences
* Giving and gracefully accepting constructive feedback
* Accepting responsibility and apologizing to those affected by our mistakes,
  and learning from the experience
* Focusing on what is best not just for us as individuals, but for the overall
  community

Examples of unacceptable behavior include:

* The use of sexualized language or imagery, and sexual attention or advances of
  any kind
* Trolling, insulting or derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or email address,
  without their explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

## Enforcement Responsibilities

Community leaders are responsible for clarifying and enforcing our standards of
acceptable behavior and will take appropriate and fair corrective action in
response to any behavior that they deem inappropriate, threatening, offensive,
or harmful.

Community leaders have the right and responsibility to remove, edit, or reject
comments, commits, code, wiki edits, issues, and other contributions that are
not aligned to this Code of Conduct, and will communicate reasons for moderation
decisions when appropriate.

## Scope

This Code of Conduct applies within all community spaces, and also applies when
an individual is officially representing the community in public spaces.
Examples of representing our community include using an official e-mail address,
posting via an official social media account, or acting as an appointed
representative at an online or offline event.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported to the community leaders responsible for enforcement at
alanakaptanoglu@gmail.com or briandesilva1@gmail.com.
All complaints will be reviewed and investigated promptly and fairly.

All community leaders are obligated to respect the privacy and security of the
reporter of any incident.

## Enforcement Guidelines

Community leaders will follow these Community Impact Guidelines in determining
the consequences for any action they deem in violation of this Code of Conduct:

### 1. Correction

**Community Impact**: Use of inappropriate language or other behavior deemed
unprofessional or unwelcome in the community.

**Consequence**: A private, written warning from community leaders, providing
clarity around the nature of the violation and an explanation of why the
behavior was inappropriate. A public apology may be requested.

### 2. Warning

**Community Impact**: A violation through a single incident or series of
actions.

**Consequence**: A warning with consequences for continued behavior. No
interaction with the people involved, including unsolicited interaction with
those enforcing the Code of Conduct, for a specified period of time. This
includes avoiding interactions in community spaces as well as external channels
like social media. Violating these terms may lead to a temporary or permanent
ban.

### 3. Temporary Ban

**Community Impact**: A serious violation of community standards, including
sustained inappropriate behavior.

**Consequence**: A temporary ban from any sort of interaction or public
communication with the community for a specified period of time. No public or
private interaction with the people involved, including unsolicited interaction
with those enforcing the Code of Conduct, is allowed during this period.
Violating these terms may lead to a permanent ban.

### 4. Permanent Ban

**Community Impact**: Demonstrating a pattern of violation of community
standards, including sustained inappropriate behavior, harassment of an
individual, or aggression toward or disparagement of classes of individuals.

**Consequence**: A permanent ban from any sort of public interaction within the
community.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 2.1, available at
[https://www.contributor-covenant.org/version/2/1/code_of_conduct.html][v2.1].

Community Impact Guidelines were inspired by
[Mozilla's code of conduct enforcement ladder][Mozilla CoC].

For answers to common questions about this code of conduct, see the FAQ at
[https://www.contributor-covenant.org/faq][FAQ]. Translations are available at
[https://www.contributor-covenant.org/translations][translations].

[homepage]: https://www.contributor-covenant.org
[v2.1]: https://www.contributor-covenant.org/version/2/1/code_of_conduct.html
[Mozilla CoC]: https://github.com/mozilla/diversity
[FAQ]: https://www.contributor-covenant.org/faq
[translations]: https://www.contributor-covenant.org/translations
---
title: 'PySINDy: A comprehensive Python package for robust sparse system identification'
tags:
  - Python
  - dynamical systems
  - sparse regression
  - model discovery
  - system identification
  - machine learning
authors:
  - name: Alan A. Kaptanoglu
    affiliation: 1
  - name: Brian M. de Silva
    affiliation: 2
  - name: Urban Fasel
    affiliation: 3
  - name: Kadierdan Kaheman
    affiliation: 3
  - name: Andy J. Goldschmidt
    affiliation: 1
  - name: Jared Callaham
    affiliation: 3
  - name: Charles B. Delahunt
    affiliation: 2
  - name: Zachary G. Nicolaou
    affiliation: 2
  - name: Kathleen Champion
    affiliation: 2
  - name: Jean-Christophe Loiseau
    affiliation: 4
  - name: J. Nathan Kutz
    affiliation: 2
  - name: Steven L. Brunton
    affiliation: 3
affiliations:
 - name: Department of Physics, University of Washington
   index: 1
 - name: Department of Applied Mathematics, University of Washington
   index: 2
 - name: Department of Mechanical Engineering, University of Washington
   index: 3
 - name: Arts et Métiers Institute of Technology, CNAM, DynFluid, HESAM Université
   index: 4
date: 21 October 2021
output: bookdown::html_document2
bibliography: paper.bib
---

# Summary
Automated data-driven modeling, the process of directly discovering the governing equations of a system from data, is increasingly being used across the scientific community. `PySINDy` is a Python package that provides tools for applying the sparse identification of nonlinear dynamics (SINDy) approach to data-driven model discovery. In this major update to `PySINDy`, we implement several advanced features that enable the discovery of more general differential equations from noisy and limited data. The library of candidate terms is extended for the identification of actuated systems, partial differential equations (PDEs), and implicit differential equations. Robust formulations, including the integral form of SINDy and ensembling techniques, are also implemented to improve performance for real-world data. Finally, we provide a range of new optimization algorithms, including several sparse regression techniques and algorithms to enforce and promote inequality constraints and stability. Together, these updates enable entirely new SINDy model discovery capabilities that have not been reported in the literature, such as constrained PDE identification and ensembling with different sparse regression optimizers.

# Statement of need
Traditionally, the governing laws and equations of nature have been derived from first principles and based on rigorous experimentation and expert intuition. 
In the modern era, cheap and efficient sensors have resulted in an unprecedented growth in the availability of measurement data, opening up the opportunity to perform automated model discovery using data-driven modeling. These data-driven approaches are also increasingly useful for processing and interpreting the information in these large datasets.
A number of such approaches have been developed in recent years, including the dynamic mode decomposition [@schmid2010dynamic;@Kutz2016book], Koopman theory [@Brunton2021koopman], nonlinear autoregressive algorithms [@Billings2013book], neural networks [@pathak2018model;@vlachas2018data;@Raissi2019jcp], Gaussian process regression [@raissi2017machine], operator inference and reduced-order modeling [@Benner2015siamreview;@peherstorfer2016data;@qian2020lift], genetic programming [@Bongard2007pnas;@schmidt_distilling_2009], and sparse regression [@brunton2016pnas].
These approaches have seen many variants and improvements over the years, so data-driven modeling software must be regularly updated to remain useful to the scientific community. The SINDy approach has experienced particularly rapid development, motivating this major update to aggregate these innovations into a single open-source tool that is transparent and easy to use for non-experts or scientists from other fields.

The original `PySINDy` code [@de2020pysindy] provided an implementation of the traditional SINDy method [@brunton2016pnas], which 
assumes that the dynamical evolution of a state variable $\mathbf{q}(t)\in\mathbb{R}^n$ follows an ODE described by a function $\mathbf{f}$,
\begin{equation}\label{eq:sindy_eq}
   \frac{d}{dt} \mathbf{q} = \mathbf{f}(\mathbf{q}).
\end{equation}
SINDy approximates the dynamical system $\mathbf{f}$ in Eq. \eqref{eq:sindy_eq} as a sparse combination of terms from a library of candidate basis functions $\boldsymbol{\theta}(\mathbf{q}) = [\theta_1(\mathbf{q}),\theta_2(\mathbf{q}),\dots,\theta_p(\mathbf{q})]$ 
\begin{equation}\label{eq:sindy_expansion}
\mathbf{f}(\mathbf{q})\approx \sum_{k=1}^{p}\theta_k(\mathbf{q})\boldsymbol\xi_k, \quad \text{or equivalently} \quad \frac{d}{dt}\mathbf{q} \approx \mathbf{\Theta}(\mathbf{q})\mathbf{\Xi},
\end{equation}
where $\boldsymbol{\Xi} = [\boldsymbol\xi_1,\boldsymbol\xi_2,\dots,\boldsymbol\xi_p]$ contain the sparse coefficients. In order for this strategy to be successful, a reasonably accurate approximation of $\mathbf{f}(\mathbf{q})$ should exist as a sparse expansion in the span of $\boldsymbol{\theta}$. Therefore, background scientific knowledge about expected terms in $\mathbf{f}(\mathbf{q})$ can be used to choose the library $\boldsymbol{\theta}$. 
To pose SINDy as a regression problem, we assume we have a set of state measurements sampled at time steps $t_1, ..., t_m$ and rearrange the data into the data matrix $\mathbf{Q} \in \mathbb{R}^{m\times n}$, \begin{eqnarray}\label{eq:Q_matrix}
\mathbf{Q} = \begin{bmatrix}
q_1(t_1) & q_2(t_1) & \cdots & q_n(t_1)\\
q_1(t_2) & q_2(t_2) & \cdots & q_n(t_2)\\
\vdots & \vdots & \ddots & \vdots \\
q_1(t_m) & q_2(t_m) & \cdots & q_n(t_m)
\end{bmatrix}
\label{Eq:DataMatrix}.
\end{eqnarray}
A matrix of derivatives in time, $\mathbf Q_t$, is defined similarly and can be numerically computed from $\mathbf{Q}$. PySINDy defaults to second order finite differences for computing derivatives, although a host of more sophisticated methods are now available, including arbitrary order finite differences, Savitzky-Golay derivatives (i.e. polynomial-filtered derivatives), spectral derivatives with optional filters, arbitrary order spline derivatives, and total variational derivatives [@ahnert2007numerical;@chartrand2011numerical;@tibshirani2011solution].

After $\mathbf Q_t$ is obtained, Eq. \eqref{eq:sindy_expansion} becomes $\mathbf Q_t \approx \mathbf{\Theta}(\mathbf{Q})\mathbf{\Xi}$ and the goal of the SINDy sparse regression problem is to choose a sparse set of coefficients $\mathbf{\Xi}$ that accurately fits the measured data in $\mathbf Q_t$. We can promote sparsity in the identified coefficients via a sparse regularizer $R(\mathbf{\Xi})$, such as the $l_0$ or $l_1$ norm, and use a sparse regression algorithm such as SR3 [@champion2020unified] to solve the resulting optimization problem,
\begin{equation}\label{eq:sindy_regression}
  \text{argmin}_{\boldsymbol\Xi}\|\mathbf Q_t - \boldsymbol\Theta(\mathbf{Q}) \boldsymbol\Xi\|^2 + R(\boldsymbol\Xi).
\end{equation}

The original `PySINDy` package was developed to identify a particular class of systems described by Eq. \eqref{eq:sindy_eq}.
Recent variants of the SINDy method are available that address systems with control inputs and model predictive control (MPC) [@Kaiser2018prsa;@fasel2021sindy], systems with physical constraints [@Loiseau2017jfm;@kaptanoglu2020physics], implicit ODEs [@mangan2016inferring;@kaheman2020sindy], PDEs [@Rudy2017sciadv;@Schaeffer2017prsa], and weak form ODEs and PDEs [@Schaeffer2017pre;@Reinbold2020pre;@messenger2021weakpde]. Other methods, such as ensembling and sub-sampling [@maddu2019stability;@reinbold2021robust;@delahunt2021toolkit], are often vital for making the identification of Eq. \eqref{eq:sindy_eq} more robust. 
In order to incorporate these new developments and accommodate the wide variety of possible dynamical systems, we have extended `PySINDy` to a more general setting and added significant new functionality. Our code\footnote{\url{https://github.com/dynamicslab/pysindy}} is thoroughly documented, contains extensive examples, and integrates a wide range of functionality, some of which may be found in a number of other local SINDy implementations\footnote{\url{https://github.com/snagcliffs/PDE-FIND}, \url{https://github.com/eurika-kaiser/SINDY-MPC},\\ \url{https://github.com/dynamicslab/SINDy-PI}, \url{https://github.com/SchatzLabGT/SymbolicRegression},\\ \url{https://github.com/dynamicslab/databook_python}, \url{https://github.com/sheadan/SINDy-BVP},\\ \url{https://github.com/sethhirsh/BayesianSindy}, \url{https://github.com/racdale/sindyr},\\ \url{https://github.com/SciML/DataDrivenDiffEq.jl}, \url{https://github.com/MathBioCU/WSINDy_PDE},\\ \url{https://github.com/pakreinbold/PDE_Discovery_Weak_Formulation}, \url{https://github.com/ZIB-IOL/CINDy}}. In contrast to some of these existing codes, `PySINDy` is completely open-source, professionally-maintained (for instance, providing unit tests and adhering to PEP8 stylistic standards), and minimally dependent on non-standard Python packages.

# New features
Given spatiotemporal data $\mathbf{Q}(\mathbf{x}, t) \in \mathbb{R}^{m\times n}$, and optional control inputs $\mathbf{u} \in \mathbb{R}^{m \times r}$ (note $m$ has been redefined here to be the product of the number of spatial measurements and the number of time samples), `PySINDy` can now approximate algebraic systems of PDEs (and corresponding weak forms) in an arbitrary number of spatial dimensions. Assuming the system is described by a function $\mathbf{g}$, we have
\begin{equation}\label{eq:pysindy_eq}
    \mathbf{g}(\mathbf{q},\mathbf q_t, \mathbf q_x, \mathbf q_y, \mathbf q_{xx}, ..., \mathbf{u}) = 0.
\end{equation}
ODEs, implicit ODEs, PDEs, and other dynamical systems are subsets of Eq. \eqref{eq:pysindy_eq}. We can accommodate control terms and partial derivatives in the SINDy library by adding them as columns in $\mathbf{\Theta}(\mathbf{Q})$, which becomes $\mathbf{\Theta}(\mathbf{Q}, \mathbf Q_t, \mathbf Q_x, ..., \mathbf{u})$. 

In addition, we have extended `PySINDy` to handle more complex modeling scenarios, including trapping SINDy for provably stable ODE models for fluids [@kaptanoglu2021promoting], models trained using multiple dynamic trajectories, and the generation of many models with sub-sampling and ensembling methods [@fasel2021ensemble] for cross-validation and probabilistic system identification. In order to solve Eq. \eqref{eq:pysindy_eq}, `PySINDy` implements several different sparse regression algorithms. Greedy sparse regression algorithms, including step-wise sparse regression (SSR) [@boninsegna2018sparse] and forward regression orthogonal least squares (FROLS) [@Billings2013book], are now available. For maximally versatile candidate libraries, the new `GeneralizedLibrary` class allows for tensoring, concatenating, and otherwise combining many different candidate libraries, along with optionally specifying a subset of the inputs to use for generating each of the libraries. \autoref{fig:package-structure} illustrates the `PySINDy` code structure, changes, and high-level goals for future work, and [`YouTube` tutorials](https://www.youtube.com/playlist?list=PLN90bHJU-JLoOfEk0KyBs2qLTV7OkMZ25) for this new functionality are available online.

`PySINDy` includes extensive Jupyter notebook tutorials that demonstrate the usage of various features of the package and reproduce nearly the entirety of the examples from the original SINDy paper [@brunton2016pnas], trapping SINDy paper [@kaptanoglu2021promoting], and the PDE-FIND paper [@Rudy2017sciadv]. 
We include an extended example for the quasiperiodic shear-driven cavity flow [@callaham2021role].
As a simple illustration of the new functionality, we demonstrate how SINDy can be used to identify the Kuramoto-Sivashinsky (KS) PDE from data. We train the model on the first 60\% of the data from Rudy et al. [@Rudy2017sciadv], which in total contains 1024 spatial grid points and 251 time steps. The KS model is identified correctly and the prediction for $\dot{\mathbf{q}}$ on the remaining testing data indicates strong performance in \autoref{fig:pde_id}. Lastly, we provide a useful flow chart in \autoref{fig:flow_chart} so that users can make informed choices about which advanced methods are suitable for their datasets. 

# Conclusion
The goal of the `PySINDy` package is to enable anyone with access to measurement data to engage in scientific model discovery. The package is designed to be accessible to inexperienced users, adhere to `scikit-learn` standards, include most of the existing SINDy variations in the literature, and provide a large variety of functionality for more advanced users. We hope that researchers will use and contribute to the code in the future, pushing the boundaries of what is possible in system identification.

# Acknowledgments
`PySINDy` is a fork of [`sparsereg`](https://github.com/Ohjeah/sparsereg) [@markus_quade_sparsereg].
SLB, AAK, KK, and UF acknowledge support from the Army Research Office (ARO  W911NF-19-1-0045). JLC acknowledges support from funding support from the Department of Defense (DoD) through the National Defense Science \& Engineering Graduate (NDSEG) Fellowship Program. ZGN is a Washington Research Foundation Postdoctoral Fellow.

![Summary of SINDy features organized by (a) `PySINDy` structure and (b) functionality. (a) Hierarchy from the sparse regression problem solved by SINDy, to the submodules of `PySINDy`, to the individual optimizers, libraries, and differentiation methods implemented in the code.
(b) Flow chart for organizing the SINDy variants and functionality in the literature. Bright color boxes indicate the features that have been implemented through this work, roughly organized by functionality. Semi-transparent boxes indicate features that have not yet been implemented.\label{fig:package-structure}](Fig1.png)

![`PySINDy` can now be used for PDE identification; we illustrate this new capability by accurately capturing a set of testing data from the Kuramoto-Sivashinsky system, described by $q_t = -qq_x - q_{xx} - q_{xxxx}$. The identified model is $q_t = -0.98qq_x -0.99q_{xx} - 1.0q_{xxxx}$.\label{fig:pde_id}](Fig2.png)

![This flow chart summarizes how `PySINDy` users can start with a dataset and systematically choose the proper candidate library and sparse regression optimizer that are tailored for a specific scientific task. \label{fig:flow_chart}](Fig3.png)

# References
---
title: 'PySINDy: A Python package for the sparse identification of nonlinear dynamical systems from data'
tags:
  - Python
  - dynamical systems
  - sparse regression
  - model discovery
  - system identification
  - machine learning
authors:
  - name: Brian M. de Silva
    affiliation: 1
  - name: Kathleen Champion
    affiliation: 1
  - name: Markus Quade
    affiliation: 2
  - name: Jean-Christophe Loiseau
    affiliation: 3
  - name: J. Nathan Kutz
    affiliation: 1
  - name: Steven L. Brunton
    affiliation: "4, 1"
affiliations:
 - name: Department of Applied Mathematics, University of Washington
   index: 1
 - name: Ambrosys GmbH
   index: 2
 - name: École Nationale Supérieure des Arts et Métiers
   index: 3
 - name: Department of Mechanical Engineering, University of Washington
   index: 4
date: 11 February 2020
bibliography: paper.bib
---

# Summary

Scientists have long quantified empirical observations by developing mathematical models that characterize the observations, have some measure of interpretability, and are capable of making predictions.
Dynamical systems models in particular have been widely used to study, explain, and predict system behavior in a wide range of application areas, with examples ranging from Newton's laws of classical mechanics to the Michaelis-Menten kinetics for modeling enzyme kinetics.
While governing laws and equations were traditionally derived by hand, the current growth of available measurement data and resulting emphasis on data-driven modeling motivates algorithmic approaches for model discovery.
A number of such approaches have been developed in recent years and have generated widespread interest, including Eureqa [@Schmidt81], sure independence screening and sparsifying operator [@PhysRevMaterials.2.083802], and the sparse identification of nonlinear dynamics (SINDy) [@brunton2016pnas].
Maximizing the impact of these model discovery methods requires tools to make them widely accessible to scientists across domains and at various levels of mathematical expertise.

`PySINDy` is a Python package for the discovery of governing dynamical systems models from data.
In particular, `PySINDy` provides tools for applying the SINDy approach to model discovery [@brunton2016pnas].
Given data in the form of state measurements $\mathbf{x}(t) \in \mathbb{R}^n$, the SINDy method seeks a function $\mathbf{f}$ such that
$$\frac{d}{dt}\mathbf{x}(t) = \mathbf{f}(\mathbf{x}(t)).$$
SINDy poses this model discovery as a sparse regression problem, wherein relevant terms in $\mathbf{f}$ are selected from a library of candidate functions.
Thus, SINDy models balance accuracy and efficiency, resulting in parsimonious models that avoid overfitting while remaining interpretable and generalizable.
This approach is straightforward to understand and can be readily customized using different sparse regression algorithms or library functions.

The `PySINDy` package is aimed at researchers and practitioners alike, enabling anyone with access to measurement data to engage in scientific model discovery.
The package is designed to be accessible to inexperienced practitioners, while also including options that allow more advanced users to customize it to their needs.
A number of popular SINDy variants are implemented, but `PySINDy` is also designed to enable further extensions for research and experimentation.
The package follows object-oriented design and is `scikit-learn` compatible.

The SINDy method has been widely applied for model identification in applications such as chemical reaction dynamics [@Hoffmann2018], nonlinear optics [@Sorokina2016oe], thermal fluids [@Loiseau2019data], plasma convection [@Dam2017pf], numerical algorithms [@Thaler2019jcp], and structural modeling [@lai2019sparse].
It has also been extended to handle  more complex modeling scenarios such as partial differential equations [@Schaeffer2017prsa;@Rudy2017sciadv], systems with inputs or control [@Kaiser2018prsa], corrupt or limited data [@tran2017exact;@schaeffer2018extracting], integral formulations [@Schaeffer2017pre;@Reinbold2020pre], physical constraints [@Loiseau2017jfm], tensor representations [@Gelss2019mindy], and stochastic systems [@boninsegna2018sparse].
However, there is not a definitive standard implementation or package for applying SINDy.
Versions of SINDy have been implemented within larger projects such as `sparsereg` [@markus_quade_sparsereg], but no specific implementation has emerged as the most widely adopted and most versions implement only a limited set of features.
Researchers have thus typically written their own implementations, resulting in duplicated effort and a lack of standardization.
This not only makes it more difficult to apply SINDy to scientific data sets, but also makes it more challenging to benchmark extensions to the method against the original and makes such extensions less accessible to end users.
The `PySINDy` package provides a dedicated central codebase where many of the basic SINDy features are implemented, allowing for easy use and standardization.
This also makes it straightforward for users to extend the package in a way such that new developments are available to a wider user base.


# Features
The core object in the `PySINDy` package is the `SINDy` model class, which is implemented as a `scikit-learn` estimator.
This design was chosen to make the package simple to use for a wide user base, as many potential users will be familiar with `scikit-learn`.
It also expresses the `SINDy` model object at the appropriate level of abstraction so that users can embed it into more complicated pipelines in `scikit-learn`, such as tools for parameter tuning and model selection.

Applying `SINDy` involves making several modeling decisions, namely: which numerical differentiation method is used, which functions make up the feature library, and which sparse regression algorithm is applied to learn the model.
The core `SINDy` object uses a set of default options but can be easily customized using a number of common approaches implemented in `PySINDy`.
The package provides a few standard options for numerical differentiation (finite difference and smoothed finite difference), feature libraries (polynomial and Fourier libraries, as well as a class for creating custom libraries), and sparse regression techniques (sequentially thresholded least squares [@brunton2016pnas], LASSO [@10.2307/2346178], and sparse relaxed regularized regression [@zheng2018ieee]).
Users can also create their own differentiation, sparse regression, or feature library objects for further customization.

The software package includes tutorials in the form of Jupyter notebooks.
These tutorials demonstrate the usage of various features in the package and reproduce the examples from the original SINDy paper [@brunton2016pnas].


# Acknowledgments

This project is a fork of [`sparsereg`](https://github.com/Ohjeah/sparsereg) [@markus_quade_sparsereg].
SLB acknowledges funding support from the Air Force Office of Scientific Research (AFOSR FA9550-18-1-0200) and the Army Research Office (ARO W911NF-19-1-0045).
JNK acknowledges support from the Air Force Office of Scientific Research (AFOSR FA9550-17-1-0329).
This material is based upon work supported by the National Science Foundation Graduate Research Fellowship under Grant Number DGE-1256082.

# References
PySINDy
=========

|BuildCI| |RTD| |PyPI| |Codecov| |JOSS1| |JOSS2| |DOI|

**PySINDy** is a sparse regression package with several implementations for the Sparse Identification of Nonlinear Dynamical systems (SINDy) method introduced in Brunton et al. (2016a), including the unified optimization approach of Champion et al. (2019), SINDy with control from Brunton et al. (2016b), Trapping SINDy from Kaptanoglu et al. (2021), SINDy-PI from Kaheman et al. (2020), PDE-FIND from Rudy et al. (2017), and so on. A comprehensive literature review is given in de Silva et al. (2020) and Kaptanoglu, de Silva et al. (2021).

.. contents:: Table of contents

System identification
---------------------
System identification refers to the process of leveraging measurement data to infer governing equations, in the form of dynamical systems, describing the data. Once discovered, these equations can make predictions about future states, can inform control inputs, or can enable the theoretical study using analytical techniques.
Dynamical systems are a flexible, well-studied class of mathematical objects for modeling systems evolving in time.
SINDy is a model discovery method which uses *sparse regression* to infer nonlinear dynamical systems from measurement data.
The resulting models are inherently *interpretable* and *generalizable*.

How it works
^^^^^^^^^^^^
Suppose, for some physical system of interest, we have measurements of state variables ``x(t)`` (a vector of length n) at different points in time. Examples of state variables include the position, velocity, or acceleration of objects; lift, drag, or angle of attack of aerodynamic objects; and concentrations of different chemical species. If we suspect that the system could be well-modeled by a dynamical system of the form

.. code-block:: text

    x'(t) = f(x(t)),

then we can use SINDy to learn ``f(x)`` from the data (``x'(t)`` denotes the time derivative of ``x(t)``). Note that both ``f(x)`` and ``x(t)`` are typically vectors. The fundamental assumption SINDy employs is that each component of ``f(x)``, ``f_i(x)`` can be represented as a *sparse* linear combination of basis functions ``theta_j(x)``

.. code-block:: text

    f_i(x) = theta_1(x) * xi_{1,i} + theta_2(x) * xi_{2,i} + ... + theta_k * xi{k,i}

Concatenating all the objects into matrices (denoted with capitalized names) helps to simplify things.
To this end we place all measurements of the state variables into a data matrix ``X`` (with a row per time measurement and a column per variable), the derivatives of the state variables into a matrix ``X'``, all basis functions evaluated at all points in time into a matrix ``Theta(X)`` (each basis function gets a column), and all coefficients into a third matrix ``Xi`` (one column per state variable).
The approximation problem to be solved can then be compactly written as

.. code-block:: text

    X' = Theta(X) * Xi.

Each row of this matrix equation corresponds to one coordinate function of ``f(x)``.
SINDy employs sparse regression techniques to find a solution ``Xi`` with sparse column vectors.
For a more in-depth look at the mathematical foundations of SINDy, please see our `introduction to SINDy <https://pysindy.readthedocs.io/en/latest/examples/2_introduction_to_sindy.html>`__.

Relation to PySINDy
^^^^^^^^^^^^^^^^^^^
The PySINDy package revolves around the ``SINDy`` class which consists of three primary components; one for each term in the above matrix approximation problem.

* ``differentiation_method``: computes ``X'``, though if derivatives are known or measured directly, they can be used instead
* ``feature_library``: specifies the candidate basis functions to be used to construct ``Theta(X)``
* ``optimizer``: implements a sparse regression method for solving for ``Xi``

Once a ``SINDy`` object has been created it must be fit to measurement data, similar to a ``scikit-learn`` model. It can then be used to predict derivatives given new measurements, evolve novel initial conditions forward in time, and more. PySINDy has been written to be as compatible with ``scikit-learn`` objects and methods as possible.

Example
^^^^^^^
Suppose we have measurements of the position of a particle obeying the following dynamical system at different points in time

.. code-block:: text

  x' = -2x
  y' = y

Note that this system of differential equations decouples into two differential equations whose solutions are simply ``x(t) = x_0 * exp(-2 * t)`` and ``y(t) = y_0 * exp(t)``, where ``x_0 = x(0)`` and ``y_0 = y(0)`` are the initial conditions.

Using the initial conditions ``x_0 = 3`` and ``y_0 = 0.5``, we construct the data matrix ``X``.

.. code-block:: python

  import numpy as np
  import pysindy as ps

  t = np.linspace(0, 1, 100)
  x = 3 * np.exp(-2 * t)
  y = 0.5 * np.exp(t)
  X = np.stack((x, y), axis=-1)  # First column is x, second is y

To instantiate a ``SINDy`` object with the default differentiation method, feature library, and optimizer and then fit it to the data, we invoke

.. code-block:: python

  model = ps.SINDy(feature_names=["x", "y"])
  model.fit(X, t=t)

We use the ``feature_names`` argument so that the model prints out the correct labels for ``x`` and ``y``. We can inspect the governing equations discovered by the model and check whether they seem reasonable with the ``print`` function.

.. code-block:: python

  model.print()

which prints the following

.. code-block:: text

  x' = -2.000 x
  y' = 1.000 y

PySINDy provides numerous other features not shown here. We recommend the `feature overview <https://pysindy.readthedocs.io/en/latest/examples/1_feature_overview.html>`__ section of the documentation for a more exhaustive summary of additional features.

Installation
------------

Installing with pip
^^^^^^^^^^^^^^^^^^^

If you are using Linux or macOS you can install PySINDy with pip:

.. code-block:: bash

  pip install pysindy

Installing from source
^^^^^^^^^^^^^^^^^^^^^^
First clone this repository:

.. code-block:: bash

  git clone https://github.com/dynamicslab/pysindy.git

Then, to install the package, run

.. code-block:: bash

  pip install .

If you do not have pip you can instead use

.. code-block:: bash

  python setup.py install

If you do not have root access, you should add the ``--user`` option to the above lines.

Caveats
^^^^^^^
If you would like to use the ``SINDy-PI`` optimizer, the ``Trapping SINDy`` optimizer (TrappingSR3), or the other SR3 optimizations with inequality constraints, you will also need to install the cvxpy package, e.g. with ``pip install cvxpy``.

To run the unit tests, example notebooks, or build a local copy of the documentation, you should install the additional dependencies in ``requirements-dev.txt``

.. code-block:: bash

  pip install -r requirements-dev.txt


Documentation
-------------
The documentation site for PySINDy can be found `here <https://pysindy.readthedocs.io/en/latest/>`__. There are numerous `examples <https://pysindy.readthedocs.io/en/latest/examples/index.html>`_ of PySINDy in action to help you get started. Examples are also available as `Jupyter notebooks <https://github.com/dynamicslab/pysindy/tree/master/examples>`__. A video overview of PySINDy can be found on `Youtube <https://www.youtube.com/watch?v=DvbbXX8Bd90>`__. We have also created a `video playlist <https://www.youtube.com/playlist?list=PLN90bHJU-JLoOfEk0KyBs2qLTV7OkMZ25>`__ with practical PySINDy tips.

PySINDy implements a lot of advanced functionality that may be overwhelming for new users or folks who are unfamiliar with these methods. Below (see here if image does not render https://github.com/dynamicslab/pysindy/blob/master/docs/JOSS2/Fig3.png), we provide a helpful flowchart for figuring out which methods to use, given the characteristics of your dataset:

.. image:: https://github.com/dynamicslab/pysindy/blob/master/docs/JOSS2/Fig3.png

This flow chart summarizes how `PySINDy` users can start with a dataset and systematically choose the proper candidate library and sparse regression optimizer that are tailored for a specific scientific task. The `GeneralizedLibrary` class allows for tensoring, concatenating, and otherwise combining many different candidate libraries.

Community guidelines
--------------------

Contributing examples
^^^^^^^^^^^^^^^^^^^^^
We love seeing examples of PySINDy being used to solve interesting problems! If you would like to contribute an example, reach out to us by creating an issue.

Contributing code
^^^^^^^^^^^^^^^^^
We welcome contributions to PySINDy. To contribute a new feature please submit a pull request. To get started we recommend installing the packages in ``requirements-dev.txt`` via

.. code-block:: bash

    pip install -r requirements-dev.txt

This will allow you to run unit tests and automatically format your code. To be accepted your code should conform to PEP8 and pass all unit tests. Code can be tested by invoking

.. code-block:: bash

    pytest

We recommend using ``pre-commit`` to format your code. Once you have staged changes to commit

.. code-block:: bash

    git add path/to/changed/file.py

you can run the following to automatically reformat your staged code

.. code-block:: bash

    pre-commit

Note that you will then need to re-stage any changes ``pre-commit`` made to your code.

There are a number of SINDy variants and advanced functionality that would be great to implement in future releases:

1. Bayesian SINDy, for instance that from Hirsh, Seth M., David A. Barajas-Solano, and J. Nathan Kutz. "Sparsifying Priors for Bayesian Uncertainty Quantification in Model Discovery." arXiv preprint arXiv:2107.02107 (2021).

2. Tensor SINDy, using the methods in Gelß, Patrick, et al. "Multidimensional approximation of nonlinear dynamical systems." Journal of Computational and Nonlinear Dynamics 14.6 (2019).

3. Stochastic SINDy, using the methods in Brückner, David B., Pierre Ronceray, and Chase P. Broedersz. "Inferring the dynamics of underdamped stochastic systems." Physical review letters 125.5 (2020): 058103.

4. Integration of PySINDy with a Python model-predictive control (MPC) code.

5. The PySINDy weak formulation is based on the work in Reinbold, Patrick AK, Daniel R. Gurevich, and Roman O. Grigoriev. "Using noisy or incomplete data to discover models of spatiotemporal dynamics." Physical Review E 101.1 (2020): 010203. It might be useful to additionally implement the weak formulation from Messenger, Daniel A., and David M. Bortz. "Weak SINDy for partial differential equations." Journal of Computational Physics (2021): 110525. The weak formulation in PySINDy is also fairly slow and computationally intensive, so finding ways to speed up the code would be great. 

6. The blended conditional gradients (BCG) algorithm for solving the constrained LASSO problem, Carderera, Alejandro, et al. "CINDy: Conditional gradient-based Identification of Non-linear Dynamics--Noise-robust recovery." arXiv preprint arXiv:2101.02630 (2021).

Reporting issues or bugs
^^^^^^^^^^^^^^^^^^^^^^^^
If you find a bug in the code or want to request a new feature, please open an issue.

Getting help
^^^^^^^^^^^^
For help using PySINDy please consult the `documentation <https://pysindy.readthedocs.io/en/latest/>`__ and/or our `examples <https://github.com/dynamicslab/pysindy/tree/master/examples>`__, or create an issue.

Citing PySINDy
--------------
PySINDy has been published in the Journal of Open Source Software (JOSS). The paper can be found `here <https://joss.theoj.org/papers/10.21105/joss.02104>`__.

If you use PySINDy in your work, please cite it using the following two references:

Brian M. de Silva, Kathleen Champion, Markus Quade, Jean-Christophe Loiseau, J. Nathan Kutz, and Steven L. Brunton., (2020). *PySINDy: A Python package for the sparse identification of nonlinear dynamical systems from data.* Journal of Open Source Software, 5(49), 2104, https://doi.org/10.21105/joss.02104

Alan A. Kaptanoglu, Brian M. de Silva, Urban Fasel, Kadierdan Kaheman, Andy J. Goldschmidt, Jared L. Callaham,   Charles  B.  Delahunt,   Zachary G. Nicolaou,   Kathleen  Champion,   Jean-Christophe  Loiseau,J. Nathan Kutz, and Steven L. Brunton. *PySINDy:  A comprehensive Python package for robust sparse system identification.* arXiv preprint arXiv:2111.08481, 2021.

Bibtex:

.. code-block:: text

    @article{desilva2020,
    doi = {10.21105/joss.02104},
    url = {https://doi.org/10.21105/joss.02104},
    year = {2020},
    publisher = {The Open Journal},
    volume = {5},
    number = {49},
    pages = {2104},
    author = {Brian de Silva and Kathleen Champion and Markus Quade and Jean-Christophe Loiseau and J. Kutz and Steven Brunton},
    title = {PySINDy: A Python package for the sparse identification of nonlinear dynamical systems from data},
    journal = {Journal of Open Source Software}
    }

Bibtex:

.. code-block:: text

      @article{kaptanoglu2021pysindy,
      title={PySINDy: A comprehensive Python package for robust sparse system identification},
      author={Alan A. Kaptanoglu and Brian M. de Silva and Urban Fasel and Kadierdan Kaheman and Andy J. Goldschmidt and Jared L. Callaham and Charles B. Delahunt and Zachary G. Nicolaou and Kathleen Champion and Jean-Christophe Loiseau and J. Nathan Kutz and Steven L. Brunton},
      year={2021},
	  Journal = {arXiv preprint arXiv:2111.08481},
      }

References
----------------------
-  de Silva, Brian M., Kathleen Champion, Markus Quade,
   Jean-Christophe Loiseau, J. Nathan Kutz, and Steven L. Brunton.
   *PySINDy: a Python package for the sparse identification of
   nonlinear dynamics from data.* arXiv preprint arXiv:2004.08424 (2020)
   `[arXiv] <https://arxiv.org/abs/2004.08424>`__

-  Kaptanoglu, Alan A., Brian M. de Silva, Urban Fasel, Kadierdan Kaheman, Andy J. Goldschmidt
   Jared L. Callaham, Charles B. Delahunt, Zachary G. Nicolaou, Kathleen Champion, 
   Jean-Christophe Loiseau, J. Nathan Kutz, and Steven L. Brunton.
   *PySINDy: A comprehensive Python package for robust sparse system identification.*
   arXiv preprint arXiv:2111.08481 (2021).
   `[arXiv] <https://arxiv.org/abs/2111.08481>`__

-  Brunton, Steven L., Joshua L. Proctor, and J. Nathan Kutz.
   *Discovering governing equations from data by sparse identification
   of nonlinear dynamical systems.* Proceedings of the National
   Academy of Sciences 113.15 (2016): 3932-3937.
   `[DOI] <http://dx.doi.org/10.1073/pnas.1517384113>`__

-  Champion, K., Zheng, P., Aravkin, A. Y., Brunton, S. L., & Kutz, J. N. (2020).
   *A unified sparse optimization framework to learn parsimonious physics-informed
   models from data.* IEEE Access, 8, 169259-169271.
   `[DOI] <https://doi.org/10.1109/ACCESS.2020.3023625>`__

-  Brunton, Steven L., Joshua L. Proctor, and J. Nathan Kutz.
   *Sparse identification of nonlinear dynamics with control (SINDYc).*
   IFAC-PapersOnLine 49.18 (2016): 710-715.
   `[DOI] <https://doi.org/10.1016/j.ifacol.2016.10.249>`__

-  Kaheman, K., Kutz, J. N., & Brunton, S. L. (2020).
   *SINDy-PI: a robust algorithm for parallel implicit sparse identification
   of nonlinear dynamics.* Proceedings of the Royal Society A, 476(2242), 20200279.
   `[DOI] <https://doi.org/10.1098/rspa.2020.0279>`__

-  Kaptanoglu, A. A., Callaham, J. L., Aravkin, A., Hansen, C. J., & Brunton, S. L. (2021).
   *Promoting global stability in data-driven models of quadratic nonlinear dynamics.*
   Physical Review Fluids, 6(9), 094401.
   `[DOI] <https://doi.org/10.1103/PhysRevFluids.6.094401>`__


Related packages
----------------
* `Deeptime <https://github.com/deeptime-ml/deeptime>`_ - A Python library for the analysis of time series data with methods for dimension reduction, clustering, and Markov model estimation.
* `PyDMD <https://github.com/mathLab/PyDMD/>`_ - A Python package using the Dynamic Mode Decomposition (DMD) for a data-driven model simplification based on spatiotemporal coherent structures. DMD is a great alternative to SINDy.
* `PySINDyGUI <https://github.com/hyumo/pysindy-gui>`_ - A slick-looking GUI for PySINDy.
* `SEED <https://github.com/M-Vause/SEED2.0>`_ - Software for the Extraction of Equations from Data: a GUI for many of the methods provided by PySINDy.

Contributors
------------
Thanks to the members of the community who have contributed to PySINDy!

+-------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------+
| `billtubbs <https://github.com/kopytjuk>`_            | Bug fix `#68 <https://github.com/dynamicslab/pysindy/issues/68>`_                                                                                          |
+-------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------+
| `kopytjuk <https://github.com/kopytjuk>`_             | Concatenation feature for libraries `#72 <https://github.com/dynamicslab/pysindy/pull/72>`_                                                                |
+-------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------+
| `andgoldschmidt <https://github.com/andgoldschmidt>`_ | `derivative <https://derivative.readthedocs.io/en/latest/>`_ package for numerical differentiation `#85 <https://github.com/dynamicslab/pysindy/pull/85>`_ |
+-------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------+

.. |BuildCI| image:: https://github.com/dynamicslab/pysindy/workflows/Build%20CI/badge.svg
    :target: https://github.com/dynamicslab/pysindy/actions?query=workflow%3A%22Build+CI%22

.. |RTD| image:: https://readthedocs.org/projects/pysindy/badge/?version=latest
    :target: https://pysindy.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. |PyPI| image:: https://badge.fury.io/py/pysindy.svg
    :target: https://badge.fury.io/py/pysindy

.. |Codecov| image:: https://codecov.io/gh/dynamicslab/pysindy/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/dynamicslab/pysindy

.. |JOSS1| image:: https://joss.theoj.org/papers/82d080bbe10ac3ab4bc03fa75f07d644/status.svg
    :target: https://joss.theoj.org/papers/82d080bbe10ac3ab4bc03fa75f07d644
    
.. |JOSS2| image:: https://joss.theoj.org/papers/10.21105/joss.03994/status.svg
    :target: https://doi.org/10.21105/joss.03994

.. |DOI| image:: https://zenodo.org/badge/186055899.svg
   :target: https://zenodo.org/badge/latestdoi/186055899
PySINDy Examples
================

This directory showcases the following examples of PySINDy in action.

`Feature overview <https://pysindy.readthedocs.io/en/latest/examples/1_feature_overview.html>`_
-----------------------------------------------------------------------------------------------------------
This notebook gives an almost exhaustive overview of the different features available in PySINDy. It's a good reference for how to set various options and work with different types of datasets.

`Introduction to SINDy <https://pysindy.readthedocs.io/en/latest/examples/2_introduction_to_sindy.html>`_
---------------------------------------------------------------------------------------------------------------------
We recommend that people new to SINDy start here. We give a gentle introduction to the SINDy method and how different steps in the algorithm are represented in PySINDy. We also show how to use PySINDy to learn a model for a simple linear differential equation.

`Original paper <https://pysindy.readthedocs.io/en/latest/examples/3_original_paper.html>`_
-------------------------------------------------------------------------------------------------------
This notebook uses PySINDy to reproduce the examples in the `original SINDy paper <https://www.pnas.org/content/pnas/113/15/3932.full.pdf>`_. Namely, it applies PySINDy to the following problems:

* Linear 2D ODE
* Cubic 2D ODE
* Linear 3D ODE
* Lorenz system
* Fluid wake behind a cylinder
* Logistic map
* Hopf system

`Scikit-learn compatibility <https://pysindy.readthedocs.io/en/latest/examples/4_scikit_learn_compatibility.html>`_
-------------------------------------------------------------------------------------------------------------------------------
Shows how PySINDy interfaces with various Scikit-learn objects.

* Cross-validation
* Sparse regressors

`Differentiation <https://pysindy.readthedocs.io/en/latest/examples/5_differentation.html>`_
---------------------------------------------------------------------------------------------------------
Explore the differentiation methods available in PySINDy on pure differentiation problems and as components in the SINDy algorithm.

`Deeptime compatibility <https://pysindy.readthedocs.io/en/latest/examples/6_deeptime_compatibility.html>`_
------------------------------------------------------------------------------------------------------------------------
See a demonstration of PySINDy objects designed to conform to the `Deeptime <https://deeptime-ml.github.io/latest/index.html>`_ API.

`Plasma physics <https://pysindy.readthedocs.io/en/latest/examples/7_plasma_example.html>`_
----------------------------------------------------------------------------------------------
Use the ``ConstrainedSR3`` optimizer to build a constrained model for the temporal POD modes of a plasma simulation.


`Trapping SINDy <https://pysindy.readthedocs.io/en/latest/examples/8_trapping_sindy_paper_examples.html>`_
----------------------------------------------------------------------------------------------
This notebook applies the ``TrappingSR3`` optimizer to various canonical fluid systems., proposed in this paper: Kaptanoglu, Alan A., et al. "Promoting global stability in data-driven models of quadratic nonlinear dynamics." Physical Review Fluids 6.9 (2021): 094401. A preprint is found here `<https://arxiv.org/abs/2105.01843>`_.

`SINDyPI <https://pysindy.readthedocs.io/en/latest/examples/9_sindypi_with_sympy.html>`_
----------------------------------------------------------------------------------------------
This notebook applies the ``SINDyPI`` optimizer to a simple implicit ODE and was originally proposed in this paper: Kaheman, Kadierdan, J. Nathan Kutz, and Steven L. Brunton. "SINDy-PI: a robust algorithm for parallel implicit sparse identification of nonlinear dynamics." Proceedings of the Royal Society A 476.2242 (2020): 20200279. 

`PDEFIND <https://pysindy.readthedocs.io/en/latest/examples/10_PDEFIND_examples.html>`_
----------------------------------------------------------------------------------------------
This notebook applies the PDEFIND algorithm (SINDy for PDE identification) to a number of PDEs, and was originally proposed in this paper: Rudy, Samuel H., et al. "Data-driven discovery of partial differential equations." Science Advances 3.4 (2017): e1602614.

`Greedy Algorithms <https://pysindy.readthedocs.io/en/latest/examples/11_SSR_FROLS_examples.html>`_
----------------------------------------------------------------------------------------------
This notebook uses the step-wise sparse regression (SSR) and forward-regression orthogonal least-squares (FROLS) algorithms, which are greedy algorithms that iteratively truncate (or add) one nonzero coefficient at each algorithm iteration. 

`Weak formulation SINDy <https://pysindy.readthedocs.io/en/latest/examples/12_weakform_SINDy_examples.html>`_
----------------------------------------------------------------------------------------------
This notebook uses SINDy to identify the weak-formulation of a system of ODEs or PDEs, adding significant robustness against noise in the data.

`Model ensembles <https://pysindy.readthedocs.io/en/latest/examples/13_ensembling.html>`_
----------------------------------------------------------------------------------------------
This notebook uses sub-sampling of the data and sub-sampling of the SINDy library to generate many models, and the user can choose how to average or otherwise combine these models together. This tends to make SINDy more robust against noisy data.

`Cavity flow <https://pysindy.readthedocs.io/en/latest/examples/14_cavity_flow.html>`_
----------------------------------------------------------------------------------------------
Demonstrates the use of SINDy to learn a model for the quasiperiodic dynamics in a shear-driven cavity at Re=7500, following Callaham, Brunton, and Loiseau (2021), preprint available here `<https://arxiv.org/pdf/2106.02409>`_.


Full table of contents
----------------------
Practical tips
==============

Here we provide pragmatic advice for using PySINDy effectively. We discuss potential pitfalls and strategies for overcoming them. We also specify how to incorporate custom methods not implemented natively in PySINDy, where applicable. The information presented here is derived from a combination of experience and theoretical considerations.

Numerical differentiation
-------------------------

Numerical differentiation is one of the core components of the SINDy method. Derivatives of measurement variables provide the targets (left-hand side :math:`\dot{X}`) for the sparse regression problem solved by SINDy:

.. math::

	\dot{X} \approx \Theta(X)\Xi.

If care is not taken in computing these derivatives, the quality of the learned model is likely to suffer.

By default, a second order finite difference method is used to differentiate input data. Finite difference methods tend to amplify noise in data. If the data are smooth (at least twice differentiable), then finite difference methods give accurate derivative approximations. When the data are noisy, they give derivative estimates with *more* noise than the original data. The following figure visualizes the impact of noise on numerical derivatives. Note that even a small amount of noise in the data can produce noticeable degradation in the quality of the numerical derivative.

.. figure:: figures/noisy_differentiation.png
	:align: center
	:alt: A toy example illustrating the effect of noise on derivatives computed with a second order finite difference method
	:figclass: align-center

	A toy example illustrating the effect of noise on derivatives computed with a second order finite difference method. Left: The data to be differentiated; :math:`y=\sin(x)` with and without a small amount of additive noise (normally distributed with mean 0 and standard deviation 0.01). Right: Derivatives of the data; the exact derivative :math:`\cos(x)` (blue), the finite difference derivative of the exact data (black, dashed), and the finite difference derivative of the noisy data.

One way to mitigate the effects of noise is to smooth the measurements before computing derivatives. The :code:`SmoothedFiniteDifference` method can be used for this purpose.
A numerical differentiation scheme with total variation regularization has also been proposed [Chartrand_2011]_ and recommended for use in SINDy [Brunton_2016]_.

Users wishing to employ their own numerical differentiation schemes have two ways of doing so. Derivatives of input measurements can be computed externally with the method of choice and then passed directly into the :code:`SINDy.fit` method via the :code:`x_dot` keyword argument. Alternatively, users can implement their own differentiation methods and pass them into the :code:`SINDy` constructor using the :code:`differentiation_method` argument. In this case, the supplied class need only have implemented a :code:`__call__` method taking two arguments, :code:`x` and :code:`t`.

Library selection
-----------------

The SINDy method assumes dynamics can be represented as a *sparse* linear combination of library functions. If this assumption is violated, the method is likely to exhibit poor performance. This issue tends to manifest itself as numerous library terms being active, often with weights of vastly different magnitudes, still resulting in poor model error.

Typically, prior knowledge of the system of interest and its dynamics should be used to make a judicious choice of basis functions. When such information is unavailable, the default class of library functions, polynomials, are a good place to start, as smooth functions have rapidly converging Taylor series. Brunton et al. [Brunton_2016]_ showed that, equipped with a  polynomial library, SINDy can recover the first few terms of the (zero-centered) Taylor series of the true right-hand side function :math:`\mathbf{f}(x)`. If one has reason to believe the dynamics can be sparsely represented in terms of Chebyshev polynomials rather than monomials, then the library should include Chebyshev polynomials.

PySINDy includes the :code:`CustomLibrary` and :code:`IdentityLibrary` objects to allow for flexibility in the library functions. When the desired library consists of a set of functions that should be applied to each measurement variable (or pair, triplet, etc. of measurement variables) in turn, the :code:`CustomLibrary` class should be used. The :code:`IdentityLibrary` class is the most customizable, but transfers the work of computing library functions over to the user. It expects that all the features one wishes to include in the library have already been computed and are present in :code:`X` before :code:`SINDy.fit` is called, as it simply applies the identity map to each variable that is passed to it. 
It is best suited for situations in which one has very specific instructions for how to apply library functions (e.g. if some of the functions should be applied to only some of the input variables).

As terms are added to the library, the underlying sparse regression problem becomes increasingly ill-conditioned. Therefore it is recommended to start with a small library whose size is gradually expanded until the desired level of performance is achieved. 
For example, a user may wish to start with a library of linear terms and then add quadratic and cubic terms as necessary to improve model performance.  
For the best results, the strength of regularization applied should be increased in proportion to the size of the library to account for the worsening condition number of the resulting linear system.

Users may also choose to implement library classes tailored to their applications. To do so one should have the new class inherit from our :code:`BaseFeatureLibrary` class. See the documentation for guidance on which functions the new class is expected to implement.

Optimization
------------
PySINDy uses various optimizers to solve the sparse regression problem. For a fixed differentiation method, set of inputs, and candidate library, there is still some variance in the dynamical system identified by SINDY, depending on which optimizer is employed.

The default optimizer in PySINDy is the sequentially-thresholded least-squares algorithm (:code:`STLSQ`). In addition to being the method originally proposed for use with SINDy, it involves a single, easily interpretable hyperparameter, and it exhibits good performance across a variety of problems.

The sparse relaxed regularized regression (:code:`SR3`) [Zheng_2018]_ [Champion_2019]_ algorithm can be used when the results of :code:`STLSQ` are unsatisfactory. It involves a few more hyperparameters that can  be tuned for improved accuracy. In particular, the :code:`thresholder` parameter controls the type of regularization that is applied. For optimal results, one may find it useful to experiment with :math:`L^0`, :math:`L^1`, and clipped absolute deviation (CAD) regularization. The other hyperparameters can be tuned with cross-validation.

Custom or third party sparse regression methods are also supported. Simply instantiate an instance of the custom object and pass it to the :code:`SINDy` constructor using the :code:`optimizer` keyword. Our implementation is compatible with any of the linear models from Scikit-learn (e.g. :code:`RidgeRegression`, :code:`Lasso`, and :code:`ElasticNet`).
See the documentation for a list of methods and attributes a custom optimizer is expected to implement. There you will also find an example where the Scikit-learn :code:`Lasso` object is used to perform sparse regression.

Regularization
--------------
Regularization, in this context, is a technique for improving the conditioning of ill-posed problems. Without regularization, one often obtains highly unstable results, with learned parameter values differing substantially for slightly different inputs. SINDy seeks weights that express dynamics as a *sparse* linear combination of library functions. When the columns of the measurement data or the library are statistically correlated, which is likely for  large libraries, the SINDy inverse problem can quickly become ill-posed. Though the sparsity constraint is a type of regularization itself, for many problems another form of regularization is needed for SINDy to learn a robust dynamical model.

In some cases regularization can be interpreted as enforcing a prior distribution on the model parameters [Bishop_2016]_.
Applying strong regularization biases the learned weights *away* from the values that would allow them to best fit the data and *toward* the values preferred by the prior distribution (e.g. :math:`L^2` regularization corresponds to a Gaussian prior).
Therefore once a sparse set of nonzero coefficients is discovered, our methods apply an extra  "unbiasing" step where *unregularized* least-squares is used to find the values of the identified nonzero coefficients.
All of our built-in methods use regularization by default.

Some general best practices regarding regularization follow. Most problems will benefit from some amount of regularization. Regularization strength should be increased as the size of the candidate right-hand side library grows. If warnings about ill-conditioned matrices are generated when :code:`SINDy.fit` is called, more regularization may help. We also recommend setting :code:`unbias` to :code:`True` when invoking the :code:`SINDy.fit` method, especially when large amounts of regularization are being applied. Cross-validation can be used to select appropriate regularization parameters for a given problem.


.. [Chartrand_2011] R. Chartrand, “Numerical differentiation of noisy, nonsmooth data,” *ISRN Applied Mathematics*, vol. 2011, 2011.

.. [Brunton_2016] S. L. Brunton, J. L. Proctor, and J. N. Kutz, “Discovering governing equations from data by sparse identification of nonlinear dynamical systems,” *Proceedings of the National Academy of Sciences*, vol. 113, no. 15, pp. 3932–3937, 2016.

.. [Zheng_2018] P. Zheng, T. Askham, S. L. Brunton, J. N. Kutz, and A. Y. Aravkin, “A unified framework for sparse relaxed regularized regression: Sr3,” *IEEE Access*, vol. 7, pp. 1404–1423, 2018.

.. [Champion_2019] K. Champion, P. Zheng, A. Y. Aravkin, S. L. Brunton, and J. N. Kutz, “A unified sparse optimization framework to learn parsimonious physics-informed models from data,” *arXiv preprint arXiv:1906.10612*, 2019.

.. [Bishop_2016] C. M. Bishop, Pattern recognition and machine learning. Springer, 2006... include:: ../README.rst


.. toctree::
   :maxdepth: 1
   :caption: User Guide

   API Documentation <api/pysindy>
   Examples <examples/index>
   Practical tips <tips>

.. toctree::
   :maxdepth: 1
   :caption: Useful links


   PySINDy @ PyPI <https://pypi.org/project/PySINDy/>
   Issue Tracker <https://github.com/dynamicslab/pysindy/issues>
