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
