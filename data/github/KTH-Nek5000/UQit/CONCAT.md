![UQit](https://raw.githubusercontent.com/KTH-Nek5000/UQit/master/docsrc/source/_static/uqit_logo.png)

## A Python Package for Uncertainty Quantification (UQ) in Computational Fluid Dynamics (CFD)
SimEx/FLOW, Engineering Mechanics, KTH Royal Institute of Technology, Stockholm, Sweden <br/>
#

### Features:
* **Sampling**:
  - Various stochastic and spectral types of samples

* **Uncertainty propagation or UQ forward problem**: 
  - generalized Polynomial Chaos Expansion (gPCE)
  - Probabilistic PCE (PPCE)

* **Global sensitivity analysis (GSA)**:
  - Sobol sensitivity indices

* **Surrogates**:
  - Lagrange interpolation
  - gPCE
  - Gaussian process regression (GPR) 

## Release Notes
### Release v 1.0.2, 27.10.2020
Source code, documentation, tests and notebooks are provided for the above-listed features. 

![`UQit`](./docsrc/source/_static/uqit_logo.png?style=centerme)

## A Python Package for Uncertainty Quantification (UQ) in Computational Fluid Dynamics (CFD)
SimEx/FLOW, Engineering Mechanics, KTH Royal Institute of Technology, Stockholm, Sweden <br/>
#

### Features:
* **Sampling**:
  - Various stochastic and spectral types of samples

* **Uncertainty propagation or UQ forward problem**: 
  - generalized Polynomial Chaos Expansion (gPCE)
  - Probabilistic PCE (PPCE)

* **Global sensitivity analysis (GSA)**:
  - Sobol sensitivity indices

* **Surrogates**:
  - Lagrange interpolation
  - gPCE
  - Gaussian process regression (GPR) 

### Installation:
`pip install UQit`

### Documentation:
The html documentation is on [this GitHub page](https://kth-nek5000.github.io/UQit/)

### How to run a test:
The tests associated to different UQ techniques are provided in `./tests/`. 
Most of the tests have also been used in the notebooks and documentation. 
For instance, to run a test from `pce_tests.py` in `./tests/`, do the following:

`python3 -c 'import pce_tests as X;X.pce_1d_test()'`

### Required libraries:
 * General  
   - [`numpy`](https://numpy.org/)
   - [`scipy`](https://www.scipy.org/)
   - [`matplotlib`](https://matplotlib.org/)
 * Optional
   - [`cvxpy`](https://www.cvxpy.org/) (for compressed sensing in PCE)
   - [`PyTorch`](https://pytorch.org/) (for GPR)
   - [`GPyTorch`](https://gpytorch.ai/) (for GPR)

### Bugs/Questions
* In case there is a bug, please feel free to open an issue on Github. 

* Qestions/comments:
  - Saleh Rezaeiravesh, salehr@kth.se <br/>
  - Ricardo Vinuesa, rvinuesa@mech.kth.se 
  - Philipp Schlatter, pschlatt@mech.kth.se <br/>

### To cite `UQit`:
* [Rezaeiravesh S., Vinuesa R., Schlatter P., UQit: A Python package for uncertainty quantification (UQ) in computational fluid dynamics (CFD). Journal of Open Source Software, 6(60), 2871, 2021.](https://joss.theoj.org/papers/10.21105/joss.02871)

### Publications using `UQit`:
* [Rezaeiravesh S., Vinuesa R., Schlatter P., On numerical uncertainties in scale-resolving simulations of canonical wall turbulence, Computers & Fluids, 227:105024, 2021.](https://www.sciencedirect.com/science/article/pii/S0045793021001900)

* [Rezaeiravesh S., Vinuesa R., Schlatter P., An Uncertainty-Quantification Framework for Assessing Accuracy, Sensitivity, and Robustness in Computational Fluid Dynamics, arXiv:2007.07071, 2020.](https://arxiv.org/abs/2007.07071)

## Release Notes
### Release v 1.0.2, 27.10.2020
Source code, documentation, tests and notebooks are provided for the above-listed features. 

`UQit` is planned to be continuously developed.
Therefore, any contribution in the code development and documentation is warmly welcomed!

* In case you are interested, please contact the authors listed in `README.md`.

* A general guide for how to contribute in an open-source code is found in this [GitHub link](https://github.com/github/docs/blob/main/CONTRIBUTING.md).
---
title: 'UQit: A Python package for uncertainty quantification (UQ) in computational fluid dynamics (CFD)'
tags:
  - Python
  - uncertainty quantification (UQ)
  - computational fluid dynamics (CFD)
authors:
  - name: Saleh Rezaeiravesh^[Corresponding author]
    orcid: 0000-0002-9610-9910
    affiliation: "1,2"
  - name: Ricardo Vinuesa
    affiliation: "1,2"
  - name: Philipp Schlatter
    affiliation: "1,2"
affiliations:
 - name: SimEx/FLOW, Engineering Mechanics, KTH Royal Institute of Technology,
   index: 1
 - name: Swedish e-Science Research Centre (SeRC), Stockholm, Sweden
   index: 2
date: 9 November 2020
bibliography: paper.bib
---


# Introduction
In computational physics, mathematical models are numerically solved and as a result, realizations for the quantities of interest (QoIs) are obtained. 
Even when adopting the most accurate numerical methods for deterministic mathematical models, the QoIs can still be up to some extent uncertain. 
Uncertainty is defined as the lack of certainty and it originates from the lack, impropriety or insufficiency of knowledge and information [@Smith:2013;@Ghanem:2017].
It is important to note that for a QoI, uncertainty is different from error which is defined as the deviation of a realization from a reference (true) value. 
In computational models, various sources of uncertainties may exist.
These include, but not limited to, the fidelity of the mathematical model (i.e., the extent by which the model can reflect the truth), the parameters in the models, initial data and boundary conditions, finite sampling time when computing the time-averaged QoIs, the way numerical errors interact and evolve, computer arithmetic, coding bugs, geometrical uncertainties, etc. 
Various mathematical and statistical techniques gathered under the umbrella of uncertainty quantification (UQ) can be exploited to assess the uncertainty in different models and their QoIs [@Smith:2013;@Ghanem:2017]. 
The UQ techniques not only facilitate systematic evaluation of validation and verification metrics, but also play a vital role in evaluation of the confidence and reliability of the data acquired in computations and experiments. 
Note that accurate accounting for such confidence intervals is crucial in data-driven engineering designs. 


In general, uncertainties can be divided into two main categories: aleatoric and epistemic [@Smith:2013]. 
The aleatoric uncertainties are random, inherent in the models, and hence, cannot be removed.
In contrast, epistemic uncertainties originate from using simplified models, insufficient data, etc.
Therefore,  they can be reduced, for instance, through improving the models. 
As a general strategy in UQ, we try to reformulate the epistemic uncertainties in terms of aleatoric uncertainties so that probabilistic approaches can be applied. 
To implement the resulting framework, we can adopt a non-intrusive point of view which does not require the computational codes, hereafter simulators, to be modified.
As a result, the UQ techniques can be combined with different features of computer experiments [@Santner:2003].


These strategies constitute the foundations of developing `UQit`, a Python package for uncertainty quantification in computational physics, in general, and computational fluid dynamics (CFD), in particular. 
In CFD, the Navier-Stokes equations are numerically integrated as a model of fluid flows. 
The flows are, in general, three-dimensional and time-dependent (unsteady) and at most of the Reynolds numbers relevant to practical applications, turbulent. 
A wide range of approaches has been used for numerical modeling of turbulence [@Sagaut:2013].
Moving from low- toward high-fidelity approaches, the nature of the uncertainties inherent in the simulations change from model-based to numerical-driven. 
Regardless of the approach, we may need to study the influence of different factors on the simulations QoIs, where UQ techniques are beneficial. 


# Statement of need \& Design
Performing different types of UQ analyses in CFD is so important that it has been considered as one of the required technologies in the NASA CFD vision 2030 [@Slotnick:2014].
In this regard, `UQit` can be seen as a good match noting that it can be (one of) the first Python-based open-source packages for UQ in CFD.
In fact, there are many similarities as well as connections between UQ and the techniques in the fields of machine learning and data sciences in which Python libraries are rich. 
These besides the flexible design of `UQit` provide a good potential for further development of `UQit` in response to different needs coming up in particular applications. 
Due to the non-intrusive nature of the implemented UQ techniques, `UQit` treats the CFD simulator as a blackbox, therefore it can be linked to any CFD simulator conditioned on having an appropriate interface.
As a possible future development, a Python VTK interface can be considered for the purpose of in-situ UQ analyses which will be suitable for large-scale simulations of fluid flows on supercomputers without the need of storing large data sets.

The documentation for each UQ technique in `UQit` starts from providing an overview of the theoretical background and introducing the main relevant references. 
These are followed by the details of implementation, instructions on how to use the method, and notebooks.
The examples in each notebook are exploited not only as a user guide, but also as a way to verify and validate the implementations in `UQit` through comparison of the results with reference data. 
Considering these aspects, `UQit` provides an appropriate environment for pedagogical purposes when it comes to practical guides to UQ approaches.  


# Features
Here, a short summary of the main UQ techniques implemented in `UQit` is given. 
In general, the methods are implemented at the highest required flexibility and they can be applied to any number of uncertain parameters. 
For the theoretical background, further details, and different applications in CFD, see our recent paper [@Rezaeiravesh:2020].

**Surrogates** play a key role in conducting non-intrusive UQ analyses and computer experiments.
They establish a functional relation between the simulator outputs (or QoIs) and the model inputs and parameters. 
The surrogates are constructed based on a limited number of training data and once constructed, they are much less expensive to run than the simulators. 
`UQit` uses different approaches to construct surrogates, including Lagrange interpolation, polynomial chaos expansion [@Xiu:2002;@Xiu:2010], and Gaussian process regression [@Rasmussen:2005;@Gramacy:2020]. 
In developing `UQit`, a high level of flexibility in constructing GPR surrogates has been considered especially when it comes to incorporating the observational uncertainties.


The goal of **uncertainty propagation** or **UQ forward problem** is to estimate how the known uncertainties in the inputs and parameters propagate into the QoIs. 
In `UQit`, these problems are efficiently handled using non-intrusive generalized polynomial chaos expansion (PCE) [@Xiu:2002;@Xiu:2010]. 
For constructing a PCE, `UQit` offers a diverse set of options for the schemes of truncating the expansion, types of parameter samples, and methods to compute the coefficients in the expansion.
For the latter, regression and projection methods can be adopted. 
As a useful feature for computationally expensive CFD simulations, the compressed sensing method can be utilized when the number of training samples is less than the number of terms in the expansion. 
By combining standard PCE and GPR, `UQit` provides the novel probabilistic PCE which is applicable to many CFD applications. 
    
    
**Global sensitivity analysis** is performed to quantify the sensitivity of the QoIs with respect to the uncertain inputs and parameters. 
Contrary to the local sensitivity analysis, in GSA all parameters are allowed to vary simultaneously and no linearization is involved in computing sensitivities. 
In `UQit`, Sobol Sensitivity Indices (main, interaction, and total) [@Sobol:2001] are computed as indicators of GSA. 

Driven by the needs, different features will be developed and added to `UQit` in future.


# Acknowledgments
This work has been supported by the EXCELLERAT project which has received funding from the European Union's Horizon 2020 research and innovation programme under grant agreement No 823691.
The financial support by the Linn&#233; FLOW Centre at KTH for SR is gratefully acknowledged.
PS and SR acknowledge financial support by the Knut and Alice Wallenberg Foundation as part of the Wallenberg Academy Fellow programme.

 
# References
============
Bibliography
============

.. [Smith13] `R. C. Smith. Uncertainty Quantification: Theory, Implementation, and Applications. Society for Industrial and Applied Mathematics, Philadelphia, PA, USA, 2013. <https://rsmith.math.ncsu.edu/UQ_TIA/>`_

.. [Ghanem17] `R. Ghanem, D. Higdon, and H. Owhadi, editors. Handbook of Uncertainty Quantification. Springer International Publishing, 2017. <https://www.springer.com/gp/book/9783319123844>`_

.. [Xiu02] `D. Xiu and G. E. Karniadakis. The Wiener–Askey polynomial chaos for stochastic differential equations. SIAM Journal on Scientific Computing, 24(2):619–644, 2002. <https://epubs.siam.org/doi/10.1137/S1064827501387826>`_

.. [Xiu05] `D. Xiu and J. S. Hesthaven. High-order collocation methods for differential equations with random inputs. SIAM Journal on Scientific Computing, 27(3):1118–1139, 2005. <https://epubs.siam.org/doi/10.1137/040615201>`_

.. [Xiu07] `D. Xiu. Efficient Collocational Approach for ParametricUncertainty Analysis. Communications in Computational Physics, 2(2): 293-309, 2007 <http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.324.2923&rep=rep1&type=pdf>`_

.. [Rezaeiravesh18] `S. Rezaeiravesh, R. Vinuesa, M. Liefvendahl, and P. Schlatter. Assessment of uncertainties in hot-wire anemometry and oil-film interferometry measurements for wall-bounded turbulent flows. European Journal of Mechanics - B/Fluids, 72:57-73, 2018. <https://www.sciencedirect.com/science/article/abs/pii/S099775461730496X>`_

.. [Rezaeiravesh20] `S. Rezaeiravesh, R. Vinuesa and P. Schlatter, An Uncertainty-Quantification Framework for Assessing Accuracy, Sensitivity, and Robustness in Computational Fluid Dynamics, arXiv:2007.07071, 2020. <https://arxiv.org/abs/2007.07071>`_

.. [Diamond16] `S. Diamond and S. Boyd. CVXPY: A Python-embedded modeling language for convex optimization. Journal of Machine Learning Research, 17(83):1{5, 2016. <https://www.cvxpy.org/index.html>`_

.. [Rasmussen05] `C. E. Rasmussen and C. K. I. Williams. Gaussian Processes for Machine Learning (Adaptive Computation and Machine Learning). The MIT Press, 2005. ISBN 026218253X. <http://www.gaussianprocess.org/gpml/>`_

.. [Goldberg98] `. W. Goldberg, C. K. I. Williams, and C. M. Bishop. Regression with input-dependent noise: A gaussian process treatment. In Proceedings of the 1997 Conference on Advances in Neural Information Processing Systems 10, NIPS 97, page 493-499, Cambridge, MA, USA, 1998. MIT Press. ISBN 0262100762. <https://www.microsoft.com/en-us/research/publication/regression-with-input-dependent-noise-a-gaussian-process-treatment/>`_

.. [Sobol01] `Sobol, I. Global sensitivity indices for nonlinear mathematical models and their monte carlo estimates. Mathematics and Computers in Simulation, 55(1):271 – 280, 2001. <https://www.sciencedirect.com/science/article/abs/pii/S0378475400002706>`_

.. [Gramacy20] `R. B. Gramacy. Surrogates: Gaussian Process Modeling, Design and Optimization for the Applied Sciences. Chapman Hall/CRC, Boca Raton, Florida, 2020. <https://bookdown.org/rbg/surrogates/>`_

.. [Santner03] `T. J. Santner, B. J. Williams, and W. I. Notz. The Design and Analysis of Computer Experiments. Springer New York, 2003. <https://www.springer.com/gp/book/9781441929921>`_

.. [Owen17] `N. E. Owen. A comparison of polynomial chaos and Gaussian process emulation for uncertainty quantification in computer experiments. PhD thesis, University of Exeter, UK, 2017. <https://ore.exeter.ac.uk/repository/handle/10871/29296?show=full>`_

.. [Schobi15] `R. Schobi, B. Sudret, and J. Wiart. Polynomial-chaos-based Kriging. International Journal for Uncertainty Quantification, 5(2):171{193, 2015. <http://www.dl.begellhouse.com/journals/52034eb04b657aea,65319583582efa6d,26fcd479064bfbc7.html>`_

.. [Eldred09] `M. Eldred and J. Burkardt. Comparison of non-intrusive polynomial chaos and stochastic collocation methods for uncertainty quantification. In 47th AIAA Aerospace Sciences Meeting including The New Horizons Forum and Aerospace Exposition. American Institute of Aeronautics and Astronautics, Jan. 2009. <https://arc.aiaa.org/doi/10.2514/6.2009-976>`_

.. [Gardner18] `J. R. Gardner, G. Pleiss, D. Bindel, K. Q. Weinberger, and A. G. Wilson. Gpytorch: Blackbox matrix-matrix gaussian process inference with GPU acceleration. CoRR, abs/1809.11165, 2018. <https://arxiv.org/abs/1809.11165>`_ 

.. [Canuto87] `Canuto C., Hussaini M. Y., Quarteroni A., Tang T. A., Spectral Methods in Fluid Dynamics, Springer-Verlag 1987. <https://link.springer.com/book/10.1007/978-3-642-84108-8>`_

.. UQit documentation master file, created by
   sphinx-quickstart on Wed Aug  5 18:11:47 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. figure:: ./_static/uqit_logo.png
   :scale: 70%
   :align: center
   :alt: A Python Package for Uncertainty Quantification in CFD
   
   **A Python Package for Uncertainty Quantification (UQ) in Computational Fluid Dynamics (CFD)**


:code:`UQit` is a Python package for uncertainty quantification (UQ) in computational physics, in general, and in computational fluid dynamics (CFD), in particular.
Different techniques are included to address various types of UQ analyses [Smith13]_, [Ghanem17]_ particularly arising in CFD.
The target CFD problems are in general, three-dimensional, unsteady, and computationally expensive.
These put constraints on the nature of UQ techniques which could be appropriate. 
:code:`UQit` is designed to be non-intrusively linked to any CFD solver through appropriate interfaces. 
Another important design concept in :code:`UQit` is to facilitate adding new techniques upon need and also provide the possibility of combining different UQ tools with each other and also with machine learning and data science techniques which can be easily added to :code:`UQit`.
Some of the main features in the current version of :code:`UQit` are listed below with a general overview and terminology in :ref:`overview-sect`.
Note that for each technique listed below, there is a short theory section followed by implementation details, example and a notebook. 



Licence
-------
:code:`UQit` is distributed under the terms of this `LICENSE <../../../LICENSE>`_. 

Contact
-------
We would appreciate to hear your comments, ideas, and feedbacks about :code:`UQit`. 
In case there is a bug, please feel free to open an issue on `Github`. 


.. toctree::
   :maxdepth: 3
   :caption: Contents:

   ./UQit_core_/instl_dep
   ./UQit_core_/codes_list
   ./UQit_core_/terminology
   ./UQit_core_/sampling
   ./UQit_core_/surrogate
   ./UQit_core_/uqFWD
   ./UQit_core_/gsa
   ./UQit_core_/others
   ./bib

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

=======================
Overview \& Terminology 
=======================

.. figure:: ../_static/uqFrame.png
   :scale: 40%
   :align: center
   :alt: Schematic of different UQ problems, taken from [Rezaeiravesh18]_.

   Schematic of different UQ problems, taken from [Rezaeiravesh18]_.

.. _overview-sect:  

Overview
--------
A computational model or code depends on different types of inputs and parameters which according to [Santner03]_ can be categorized as controlled, environmental and uncertain. 
The focus of the uncertainty quantification (UQ) techniques is mainly on the last two. 

:code:`UQit` is designed mainly based on the needs for UQ in the CFD (computational fluid dynamics) community. 
Its connection with the CFD solvers is non-intrusive where the CFD code is treated as a blackbox. 
As a result, we always deal with discrete data which comprises of parameter samples and corresponding responses acquired by running the simulator. 
A good overview over different UQ approaches can be found for instance in [Smith13]_ and [Ghanem17]_.
Moreover, the terminology and some of materials in this documentation are taken from [Rezaeiravesh20]_.
Below, we list some of the main features of :code:`UQit`.


* **Uncertainty propagation or UQ forward problem:**

Estimates how the known uncertainties in the inputs and parameters propagate into the quantities of interest (QoIs).
These problems can be efficiently handled using non-intrusive generalized polynomial chaos expansion (PCE), see [Xiu02]_, [Xiu07]_.
In :code:`UQit`, for constructing PCE both regression and projection methods are implemented.
Using compressed sensing method, PCE can be constructed using a small number of training samples.
Samples from the parameter space can be taken using different methods implemented in :ref:`sampling_sect` module.
See the details in :ref:`uqFwd-sect`.

* **Global sensitivity analysis (GSA):**


GSA is performed to quantify the sensitivity of the QoIs to the simultaneous variation of the inputs/parameters.
Contrary to local sensitivity analysis (LSA), in GSA all parameters are allowed to vary simultaneously and no linearization is involved in computing sensitivities.
In :code:`UQit`, GSA is performed by :ref:`sobol-sect` [Sobol01]_.

* **Surrogates:**

:code:`UQit` uses different approaches including Lagrange interpolation, polynomial chaos expansion and more importantly Gaussian process regression [Rasmussen05]_, [Gramacy20]_ to construct :ref:`surrogates-sect` which connect the QoIs to the inputs/parameters.
Surrogates are the pillars for conducting computer experiments [Santner03]_.
In particular, highest possible flexibility in constructing GPR surrogates have been considered when it comes to incorporating the observational uncertainties.



Nomenclature
------------
Throughout this documentation, we adopt the terminologies and nomenclature from [Rezaeiravesh20]_, as summarized in the following table. 

======================== =============================================
      **Symbol**                       **Definition**
------------------------ ---------------------------------------------
QoI                      Quantity of Interest
:math:`f(\cdot)`         Model function or simulator
:math:`\tilde{f}(\cdot)` Surrogate
:math:`\chi`             Controlled parameter
:math:`q_i`              i-th uncertain parameter (single-variate)
:math:`\mathbf{q}`       Multivariate uncertain parameter
:math:`\mathbf{q}^{(j)}` j-th sample of :math:`\mathbf{q}`
:math:`p`                Dimension of :math:`\mathbf{q}`
:math:`\mathbb{Q}`       Admissible space of :math:`\mathbf{q}`
:math:`\mathbb{Q}_i`     Admissible space of :math:`q_i`
:math:`r`                Model response, output or QoI
:math:`\bigotimes`       Tensor product
:math:`\mathcal{U}`      Uniform distribution
:math:`\mathcal{N}`      Normal (Gaussian) distribution
======================== =============================================

=================================
Global Sensitivity Analysis (GSA)
=================================
Global sensitivity analysis (GSA) aims at quantifying the sensitivity of the model response
or quantity of interest (QoI) with respect to the variation of the uncertain parameters and inputs. 
In other words, the influence of each of the parameters in the propagated uncertainty in the QoI is measured.
In contrast to the local sensitivity analysis, in GSA all the parameters are allowed to vary simultaneously over their admissible space, [Smith13]_. 
In :code:`UQit`, the Sobol sensitivity indices [Sobol01]_ are computed to measure GSA. 

.. _sobol-sect:

Sobol Sensitivity Indices
-------------------------

Theory
~~~~~~
The following short description has been taken from [Rezaeiravesh20]_. 
For more in depth review, the reader is referred to [Sobol01]_, Chapter 15 in [Smith13]_, and [Ghanem17]_.

By analysis of variance (`ANOVA <https://en.wikipedia.org/wiki/Analysis_of_variance>`_) or Sobol decomposition [Sobol01]_, a model function is decomposed as,

.. math::
  f(\mathbf{q}) =
  f_0 +\sum_{i=1}^p f_i(q_i) + \sum_{1\leq i<j\leq p} f_{ij}(q_i,q_j)+\cdots\,,
  \label{eq:anova_f}\tag{1}

where, :math:`f_0` is the mean of :math:`f(\mathbf{q})`, :math:`f_i(q_i)` specify the contribution of each parameter, :math:`f_{ij}(q_i,q_j)` denote effects of interaction between each pair of parameters, and so on for other interactions.
These contributors are defined as,

.. math::
   \begin{eqnarray*}
   f_0 &=& \mathbb{E}_\mathbf{q}[f(\mathbf{q})] \,, \\
   f_i(q_i) &=&\mathbb{E}_\mathbf{q}[f(\mathbf{q})|q_i] - f_0 \,, \\
   f_{ij}(q_{i},q_j) &=& \mathbb{E}_\mathbf{q}[f(\mathbf{q})|q_i,q_j] -f_i(q_i) -f_j(q_j) - f_0 \,.
   \end{eqnarray*}

Here, :math:`\mathbb{E}_\mathbf{q}[f(\mathbf{q})|q_i]`, for instance, denotes the expected value of :math:`f(\mathbf{q})` conditioned on fixed values of :math:`q_i`.
Similar to Eq. \eqref{eq:anova_f}, the total variance of :math:`f(\mathbf{q})`, denoted by :math:`D`, is decomposed as,

.. math::
   \begin{equation}
   \mathbb{V}_\mathbf{q}[f(\mathbf{q})] = D=\sum_{i=1}^p D_i + \sum_{1\leq i<j\leq p} D_{ij} + \cdots \,,
   \end{equation}

where, :math:`D_i=\mathbb{V}_\mathbf{q}[f_i(q_i)]`, :math:`D_{ij}=\mathbb{V}_\mathbf{q}[f_{ij}(q_i,q_j)]`, and so on.
The main Sobol indices are eventually defined as the contribution of each of :math:`D_i`, :math:`D_{ij}`, ... in the total variance :math:`D`:

.. math::
   \begin{equation}
   S_i=D_i/D\,,\quad
   S_{ij}=D_{ij}/D \,,\, \ldots \,, \quad i,j=1,2,\cdots,p
   \label{eq:sobol} \tag{2}
   \end{equation}

Example
~~~~~~~
Given samples :code:`q` with associated :code:`pdf` and response values :code:`f`, the Sobol indices in :code:`UQit` are computed by, 

.. code-block:: python

   sobol_=sobol(q,f,pdf)
   Si=sobol_.Si     #1st-order main indices
   STi=sobol_.STi   #1st-order total indices
   Sij=sobol_.Sij   #2nd-order interactions
   SijName=sobol_.SijName  #Names of Sij

Implementation
~~~~~~~~~~~~~~

.. automodule:: sobol
   :members:

Notebook
~~~~~~~~
Try this `GSA Notebook <../examples/sobol.ipynb>`_ to see how to use :code:`UQit` to compute Sobol indices. The provided examples can also be seen as a way to validate the implementation of the methods in :code:`UQit`.  

.. _surrogates-sect:

==========
Surrogates
==========

A surrogate or metamodel is an approximation of the actual model function or simulator over the parameter space.
Running a surrogate is computationally much less expensive than the actual simulator, a characteristic that makes the use of surrogates inevitable in different UQ problems.
However, the predictions by the surrogate should be accurate enough compared to the actual predictions by the simulator. 
Using an additive error model, the following relation can be considered between the model function (simulator) and its observations:

.. math::
   r=f(\chi,\mathbf{q})+\mathbf{\varepsilon}\,,

where :math:`\mathbf{\varepsilon}` expresses the bias and random noise.
Our task is to construct a surrogate :math:`\tilde{f}(\chi,\mathbf{q})` for the actual model function (simulator).

Treating the simulator as a blackbox, a set of training data :math:`\mathcal{D}=\{(\mathbf{q}^{(i)},r^{(i)})\}_{i=1}^n` is obtained. 
There are different techniques to construct a surrogate, for instance see [Smith13]_, [Ghanem17]_, and [Gramacy20]_. 
Some of the techniques relevant to CFD applications have been implemented in :code:`UQit`.
Here, we provide a short overview to the theory behind these techniques, explain their implementation in :code:`UQit`, and provide examples to show how to use them. 


Non-intrusive Polynomial Chaos Expansion
-----------------------------------------
As a type of stochastic collocation (SC) methods, see e.g. Chapter 20 in [Ghanem17]_, non-intrusive PCE [Xiu05]_, [Xiu07]_ can be used to construct a surrogate. 

.. math::
   \tilde{f}(\chi,\mathbf{q}) = \sum_{k=0}^K \hat{f}_k(\chi) \Psi_{k}(\xi) \,.

There is a one-to-one correspondence between any sample of :math:`\mathbf{q}\in \mathbb{Q}` and :math:`\xi\in\Gamma`, where :math:`\mathbb{Q}=\bigotimes_{i=1}^p \mathbb{Q}_i` and :math:`\Gamma=\bigotimes_{i=1}^p \Gamma_i`. 
Note that :math:`\mathbb{Q}_i` is the admissible space of the i-th parameter which can be mappd onto :math:`\Gamma_i` based on the gPCE rules, see [Xiu02]_, [Eldred09]_.
For the details of the non-intrusive PCE method refer to :ref:`gPCE-sect`.



Lagrange Interpolation
----------------------

Theory
~~~~~~
As another form of SC-based surrogates, Lagrange interpolation can be considered:

.. math::
   \tilde{f}(\chi,\mathbf{q}) = \sum_{k=1}^n \hat{f}_k(\chi,\mathbf{q}) L_k(\mathbf{q}) \,,

where :math:`\hat{f}_k(\chi,\mathbf{q})=f(\chi,\mathbf{q}^{(k)})=r^{(k)}` are the training model outputs.
If the :math:`n_i` samples taken from the :math:`i`-th parameter space are represented by :math:`Q_{n_i}`, then the use of tensor-product leads to the nodal set :math:`Q_n` of size :math:`n=\prod_{i=1}^p n_i`, where,

.. math::
   Q_n= Q_{n_1} \bigotimes Q_{n_2}\bigotimes \ldots \bigotimes Q_{n_p} \,.

Correspondingly, the Lagrange bases :math:`L_k(\mathbf{q})` are constructed using the tensor-product of the bases in each of the parameter spaces: 

.. math::
   L_k(\mathbf{q})=L_{k_1}(q_1) \bigotimes L_{k_2}(q_2) \bigotimes \ldots \bigotimes L_{k_p}(q_p) \,,

where,

.. math::
   L_{k_i}(q_i) = \prod_{\substack{{k_i=1}\\{k_i\neq j}}}^{n_i} 
   \frac{q_i - q_i^{(k_i)}}{q_i^{(k_i)}-q_i^{(j)}} \,,\quad i=1,2,\ldots,p \,.

Note that the Lagrange basis satisfies :math:`L_{k_i}(q_i^{(j)})=\delta_{k_{i}j}`, where :math:`\delta` represents the Kronecker delta. 


Example
~~~~~~~
* For :math:`p=1` (one-dimensional parameter :math:`q`):

.. code-block:: python

    fInterp=lagInt(fNodes=fNodes,qNodes=[qNodes],qTest=[qTest]).val

* For :math:`p>1` (multi-dimensional parameter :math:`\mathbf{q}`):

.. code-block:: python

    fInterp=lagInt(fNodes=fNodes,qNodes=qNodes,qTest=qTestList,liDict={'testRule':'tensorProd'}).val

Implementation
~~~~~~~~~~~~~~
.. automodule:: lagInt
   :members:

Notebook
~~~~~~~~
Try this `LagInt notebook <../examples/lagInt.ipynb>`_ to see how to use :code:`UQit` for Lagrange interpolation over a parameter space. 



Gaussian Process Regression
---------------------------
Theory
~~~~~~
Consider the simulator :math:`f(\mathbf{q})` where :math:`\mathbf{q}\in \mathbb{Q}`. 
The observations are assumed to be generated from the following model,

.. math::
   y = f(\mathbf{q}) + \varepsilon  \,.

Since the exact simulator :math:`f(\mathbf{q})` is not known, we can put a prior on it, which is in the form of a Gaussian process, see [Rasmussen05]_, [Gramacy20]_. 
Based on the training data :math:`\mathcal{D}`, the posterior of the :math:`f(q)`, denoted by :math:`\tilde{f}(\mathbf{q})`, is inferred. 
Without loss of generality we assume :math:`\varepsilon\sim\mathcal{N}(0,\sigma^2)`. 
Contrary to the common use of Gaussian process regression (GPR) where :math:`\sigma` is assumed to be fixed for all observations (homoscedastic noise), we are interested in cases where :math:`\sigma` is observation-dependent (heteroscedastic noise).
In the latter, we need to have a Gaussian process to infer the noise parameters, see [Goldberg98]_.
Eventually, the posteriors of :math:`\tilde{f}(\mathbf{q})` and response :math:`y` can be sampled over the parameter space, see [Rezaeiravesh20]_ and the references therein for the details. 


Example
~~~~~~~
Given the training data including the observational noise, A GPR is constructed in :code:`UQit` as,

.. code-block:: python

   gpr_=gpr(xTrain,yTrain[:,None],noiseSdev,xTest,gprOpts)
   post_f=gpr_.post_f
   post_obs=gpr_.post_y


Implementation
~~~~~~~~~~~~~~
In :code:`UQit`, the GPR is implemented using the existing Python library :code:`GPyTorch` [Gardner18]_. 
The user can similarly use any other available library for GPR as long as the code structure is kept consistent with :code:`UQit`. 

.. automodule:: gpr_torch
   :members:



Notebook
~~~~~~~~
Try `GPR notebook <../examples/gpr.ipynb>`_ to see how to use :code:`UQit` for Gaussian process regression over a parameter space. 

=============================
Installation and Dependencies
=============================


Installation
------------
:code:`UQit` can be found on `PyPI`, see `here <https://pypi.org/project/UQit/>`_. 

To install :code:`UQit`:

.. code-block::

   pip install UQit


The source code is found at this `GitHub repository <https://github.com/KTH-Nek5000/UQit>`_.


Documentation
-------------
The `html` documentation is found at this `GitHub page <https://kth-nek5000.github.io/UQit/>`_.   


Dependencies
------------



* Required:

  - `numpy <https://numpy.org/>`_
  - `scipy <https://www.scipy.org/>`_
  - `matplotlib <https://matplotlib.org/>`_


* Optional:

  - `cvxpy <https://www.cvxpy.org/>`_ (for compressed sensing in PCE)
  - `GPyTorch <https://gpytorch.ai/>`_ (for GPR)
  - `PyTorch <https://pytorch.org/>`_ (for GPR)




.. _sampling_sect:

=========
Sampling 
=========

Sampling
--------
In :code:`UQit`, different types of samples can be taken from the parameter space. 
From one point of view, the parameter samples are divided into training and test. 
To construct a surrogate or perform a UQ forward problem, we need to take training samples from a mapped or standardized space :math:`\Gamma` and then map them onto the parameter admissible space :math:`\mathbb{Q}`.
In contrast, the test samples which are, for instance, used to evaluate the constructed surrogates, are taken from :math:`\mathbb{Q}` and then are mapped onto :math:`\Gamma`.

Available types of training samples:

 * :code:`GQ`: Gauss-Quadrature nodes

   can be used with distributions :code:`Unif`, :code:`Norm`
 * :code:`GLL`: Gauss-Lobatto-Legendre nodes   
 * :code:`unifSpaced`: Uniformly-spaced samples   
 * :code:`unifRand`: Uniformly distributed random samples   
 * :code:`normRand`: Gaussian distributed random samples
 * :code:`Clenshaw`: Clenshaw points
 * :code:`Clenshaw-Curtis`: Clenshaw-Curtis points

Available types of test samples:

 * :code:`GLL`: Gauss-Lobatto-Legendre nodes
 * :code:`unifSpaced`: Uniformly-spaced points
 * :code:`unifRand`: Uniformly distributed random
 * :code:`normRand`: Gaussian distributed random

Note that the argument :code:`qInfo` appearing in sampling methods:
 * :code:`qInfo=[a,b]`, if the parameter is :code:`Unif` over range :math:`[a,b]`, i.e. :math:`q\sim\mathcal{U}[a,b]`
 * :code:`qInfo=[m,s]` contains the mean :math:`m` and standard-deviation :math:`s`, if the parameter is :code:`Norm`, i.e. :math:`q\sim \mathcal{N}(m,s^2)`


Example
~~~~~~~

.. code-block:: python

   tr_=trainSample(sampleType='GQ',GQdistType='Unif',qInfo=[2,3],nSamp=10) 
   tr_=trainSample(sampleType='NormRand',qInfo=[2,3],nSamp=10) 
   tr_=trainSample(sampleType='GLL',qInfo=[2,3],nSamp=10)

.. code-block:: python

   ts_=testSample(sampleType='unifRand',GQdistType='Unif',qBound=[-1,3],nSamp=10) 
   ts_=testSample(sampleType='unifRand',qBound=[-1,3],nSamp=10) 
   ts_=testSample(sampleType='normRand',GQdistType='Norm',qBound=[-1,3],qInfo=[0.5,2],nSamp=10) 
   ts_=testSample(sampleType='unifSpaced',GQdistType='Norm',qBound=[-1,3],qInfo=[0.5,2],nSamp=10) 
   ts_=testSample(sampleType='unifSpaced',GQdistType='Unif',qBound=[-1,3],nSamp=10) 
   ts_=testSample(sampleType='GLL',qBound=[-1,3],nSamp=10)


Implementation
~~~~~~~~~~~~~~

.. automodule:: sampling
   :members: 


Nodes
-----
Some of the sampling methods rely on generating nodes from mathematical polynomials, for instance see [Canuto87]_.
The associated methods are implemented in :code:`nodes.py`.

.. automodule:: nodes
   :members:
.. _uqFwd-sect:

==================
UQ Forward Problem
==================

In UQ forward (uncertainty propagation) problems, the aim is to estimate the propagation of uncertainties in the model response or QoIs, from known uncertain parameters/inputs. 
In many situations, we are mostly interested in approximately estimating the statistical moments of the model outputs and QoIs. 
In particular, among the moments, our main focus is on the expected value :math:`\mathbb{E}_\mathbf{q}[r]` and variance :math:`\mathbb{V}_\mathbf{q}[r]` of a QoI :math:`r`, with :math:`\mathbf{q}` denoting the uncertain parameters.
Different approaches can be used for this purpose, see [Smith13]_, [Ghanem17]_.
The Monte Carlo approaches rely on taking independent samples from the parameters and running the model simulator as a blackbox. 
However, the convergence rate of these methods can be low, see e.g. [Ghanem17]_.
As a more efficient technique, the spectral-based method non-intrusive generalized polynomial chaos expansion (gPCE) [Xiu02]_, [Xiu05]_, [Xiu07]_ is employed in :code:`UQit` for UQ forward problems. 
For an overview of the technique with the same notations as in this document, refer to [Rezaeiravesh20]_.

.. _gPCE-sect:

Standard Polynomial Chaos Expansion
-----------------------------------

Theory
~~~~~~
The generalized polynomial chaos expansion for :math:`\tilde{f}(\chi,\mathbf{q})` is written as,

.. math::
   \tilde{f}(\chi,\mathbf{q}) = \sum_{k=0}^K \hat{f}_k(\chi) \Psi_{k}(\xi) \,,

where the basis function for the multi-variate parameter :math:`\xi\in\Gamma` is defined as :math:`\Psi_{k}(\xi)=\prod_{i=1}^p \psi_{k_i}(\xi_i)`.
In the framework of generalized PCE, see [Xiu02]_, given the distribution of the single-variate parameter :math:`\xi_i\in \Gamma_i`, a set of orthogonal basis functions are chosen for :math:`\psi_{k_i}(\xi_i)`, for :math:`i=1,2,\ldots,p`. 
Note that, there is a one-to-one correspondence between any sample of :math:`\mathbf{q}\in \mathbb{Q}` and :math:`\xi\in\Gamma`, where :math:`\mathbb{Q}=\bigotimes_{i=1}^p \mathbb{Q}_i` and :math:`\Gamma=\bigotimes_{i=1}^p \Gamma_i`.
The mapped space :math:`\Gamma_i` is known based on the gPCE rule, see [Xiu02]_, [Eldred09]_.


Given a set of training data :math:`\mathcal{D}=\{(\mathbf{q}^{(i)},r^{(i)})\}_{i=1}^n`, there are two main steps to construct the above expansion.
First, a truncation scheme is needed to handle :math:`p`-dimensional parameter space and determine :math:`K`.
Currently, tensor-product and total-order schemes are available in :code:`UQit`. 
Second, the coefficients :math:`\{\hat{f}_k(\chi)\}_{k=0}^K` have to be determined. 
In :code:`UQit`, two different approaches can be used for this purpose: projection and regression method, see [Rezaeiravesh20]_ and the references therein. 
In case the number of training data is less than :math:`K`, compressed sensing method can be adopted which is implemented in :code:`UQit` through the external Python library :code:`cvxpy` [Diamond16]_.

Once the coefficients :math:`\{\hat{f}_k(\chi)\}_{k=0}^K` are obtained, the PCE can be used as a surrogate for the actual unobserved :math:`f(\chi,\mathbf{q})`.
A main advantage of the PCE method is that the approximate estimation of the statistical moments of the :math:`f(\chi,\mathbf{q})` or response :math:`r` is a natural outcome of the surrogate construction. 
Using gPCE, the mean and variance of the simulator are estimated by,

.. math::
   \mathbb{E}_{\mathbf{q}}[f(\chi,\mathbf{q})] = \hat{f}_0(\chi),

.. math::
   \mathbb{V}_{\mathbf{q}}[f(\chi,\mathbf{q})] = \sum_{k=1}^K \hat{f}^2_k(\chi) \gamma_k, 

where :math:`\gamma_k` is the inner-product of the polynomial basis.


Example
~~~~~~~
In :code:`UQit`, to construct and estimate expected value and variance of :math:`f(\mathbf{q})` for :math:`\mathbf{q}\in\mathbb{Q}\subset \mathbb{R}^p`, we have:

.. code-block:: python

   pce_=pce(fVal=fVal,xi=xiGrid,pceDict=pceDict,nQList=nQ)
   fMean=pce_.fMean       
   fVar=pce_.fVar         
   pceCoefs=pce_.coefs     
   kSet=pce_.kSet

To evaluate a constructed PCE at a set of test parameter samples taken from :math:`\Gamma`, we write:

.. code-block:: python

   pcePred_=pce.pceEval(coefs=pceCoefs,xi=xiTest,distType=distType,kSet=kSet)
   fPCE=pcePred_.pceVal

As for instance described in [Rezaeiravesh20]_, as an a-posteriori measure of the convergence of the PCE terms, we can evaluate the following indicator,

.. math::
   \vartheta_\mathbf{k} = |\hat{f}_\mathbf{k}| \, \|\Psi_{\mathbf{k}}(\mathbf{\xi})\|_2/|\hat{f}_0|
   
at different multi-indices :math:`\mathbf{k}`. 
In :code:`UQit` this is done through running,

.. code-block:: python

   pce.convPlot(coefs=pceCoefs,distType=distType)



Implementation
~~~~~~~~~~~~~~
In :code:`UQit`, the methods required for standard PCE are implemented in :code:`pce.py`. 

.. automodule:: pce
   :members:

Notebook
~~~~~~~~
Try the `PCE notebook <../examples/pce.ipynb>`_ to see how to use :code:`UQit` to perform standard polynomial chaos expansion (PCE). The provided examples can also be seen as a way to validate the implementation of the methods in :code:`UQit`.

Probabilistic Polynomial Chaos Expansion
----------------------------------------

Theory
~~~~~~
The standard PCE (polynomial chaos expansion) and GPR (Gaussian process regression) are two powerful approaches for surrogate construction and metamodeling in UQ. 
Combining these two approaches, probabilistic PCE is derived. 
There are at least two different views to this derivation which can be found in Schobi et al. [Schobi15]_ and Owen [Owen17]_. 
In :code:`UQit`, a generalization of the latter is implemented which is detailed in [Rezaeiravesh20]_. 


Example
~~~~~~~
Given training parameter samples :code:`qTrain` and associated responses :code:`yTrain` with observation noise :code:`noiseSdev`, the :code:`ppce` is constructed as follows. 

.. code-block:: python

   ppce_=ppce.ppce(qTrain,yTrain,noiseSdev,ppceDict)
   optOut=ppce_.optOut
   fMean_samples=ppce_.fMean_samps
   fVar_samples=ppce_.fVar_samps


Implementation
~~~~~~~~~~~~~~
.. automodule:: ppce
   :members:


Notebook
~~~~~~~~
Try this `PPCE notebook <../examples/ppce.ipynb>`_ to see how to use :code:`UQit` to perform probabilistic polynomial chaos expansion (PPCE). 


==================
List of core codes
==================

* :code:`analyticTestFuncs.py`: 
  Analytical model functions to test implementation of different techniques
* :code:`gpr_torch.py`:
  Gaussian Process Regression (GPR) using GPyTorch library
* :code:`lagInt.py`:
  Lagrange interpolation
* :code:`linAlg.py`:
  Tools for linear algebra 
* :code:`nodes.py`:
  Spectral nodes from different types of polynomials
* :code:`pce.py`: 
  generalized Polynomial Chaos Expansion
* :code:`stats.py`:
  Statistical tools
* :code:`ppce.py`:
  Probabilistic generalized Polynomial Chaos Expansion (PPCE)
* :code:`reshaper.py`:
  Tools for converting and reshaping arrays and lists  
* :code:`sampling.py`: 
  Sampling from parameter spaces
* :code:`sobol.py`:
  Sobol sensitivity indices
* :code:`surr2surr.py`:
  Interpolation from a surrogate to another surrogate
* :code:`writeUQ.py`:  
  Tools for printing or writing data in file
================
Other Core Codes
================

Analytical Test Functions
-------------------------
Analytical model functions to test and validate implementation of different UQ techniques.

.. automodule:: analyticTestFuncs
   :members:


Surrogate to Surrogate
----------------------
Interpolate values from one surrogate to another surrogate.

.. automodule:: surr2surr
   :members:


Statistical Tools
-----------------

.. automodule:: stats
   :members:


Linear Algebra
--------------
Tools for linear algebra.

To solve a linear system which is under-determined, the compressed sensing method is used. 
The required optimization is handled by :code:`cxvpy` [Diamond16]_. 
Different solvers can be used for this purpose, a list of which can be obtained by 
`cvxpy.installed_solvers()`. 
The required options for each solver can be found in `this cvxpy page <https://www.cvxpy.org/tutorial/advanced/index.html?highlight=installed_solvers>`_.
Note that the default solver is directly specified in :code:`linAlg.myLinearRegress()`.


.. automodule:: linAlg
   :members:

Reshaping Tools
---------------
Tools for converting and reshaping arrays and lists.

.. automodule:: reshaper
   :members:

Tools for Printing and Writing
------------------------------
Tools for printing or writing data in file

.. automodule:: write
   :members:

