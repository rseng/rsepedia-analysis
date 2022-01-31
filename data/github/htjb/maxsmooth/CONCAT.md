---
title: 'maxsmooth: Derivative Constrained Function Fitting'
tags:
  - Python
  - astrophysics
  - cosmology
authors:
  - name: Harry T. J. Bevins
    orcid: 0000-0002-4367-3550
    affiliation: "1"
affiliations:
 - name: Astrophysics Group, Cavendish Laboratory, J.J.Thomson Avenue, Cambridge, CB3 0HE, United Kingdom
   index: 1
date: 31 July 2020
bibliography: paper.bib

---

# Summary

``maxsmooth`` is an optimisation routine written in Python (supporting version $\geq 3.6$)
for fitting Derivative Constrained Functions (DCFs) to data.
DCFs are a family of functions which have derivatives that do not cross
zero in the band of interest. Two special cases of DCF are Maximally Smooth Functions
(MSFs) which have derivatives with order $m \geq 2$ constrained and Completely Smooth
Functions (CSFs) with $m \geq 1$ constrained. Alternatively, we can constrain an
arbitrary set of derivatives and we
refer to these models as Partially Smooth Functions.
Due to their constrained nature, DCFs can produce perfectly smooth fits to
data and reveal non-smooth signals of interest in the residuals.

# Statement of Need

The development of ``maxsmooth`` has been motivated by the problem
of foreground modelling in Global 21-cm experiments [@EDGES_LB; @LEDA; @SARAS; @REACH].
Global 21-cm cosmology is the study of the spin temperature of hydrogen gas and
its relative magnitude compared to the Cosmic Microwave Background
during the Cosmic Dawn (CD) and Epoch of Reionisation (EoR). During the CD and EoR
the first stars formed and the properties of the hydrogen gas
changed as it interacted with radiation from these stars
[@Furlanetto2006; @Pritchard2012; @Barkana2016].

The goal of Global 21-cm experiments is to detect this
structure in the sky averaged radio spectrum between $\nu = 50$ and $200$ MHz.
However, the signal of interest is expected to be approximately $\leq 250$ mK and masked by
foregrounds $10^4 - 10^5$ times brighter [@Cohen; @Cohen2019].

Modelling and removal of the foreground is essential for detection of
the Global 21-cm signal. DCFs provide a powerful alternative to unconstrained
polynomials for accurately modelling the smooth synchrotron/free-free emission
foreground from the Galaxy and extragalactic radio sources.

To illustrate the abilities of ``maxsmooth`` we produce mock
21-cm experiment data and model and remove the foreground using an MSF.
We add to a mock foreground, $\nu^{-2.5}$,
a Gaussian noise with standard deviation of $0.02$ K and a
Gaussian signal with amplitude $0.23$ K, central frequency of $100$ MHz
and standard deviation of $10$ MHz.

The figure below shows the residuals (bottom panel, green) when fitting
and removing an MSF from the data (top panel) compared to the injected signal
(bottom panel, red). While the removal of the foreground does not
directly recover the injected signal, rather a smooth baseline subtracted version,
we see the remnant of the signal in
the residuals (see @Bevins for more details and examples).

![**Top panel:** The mock 21-cm data used to illustrate the abilities of
``maxsmooth``. **Bottom Panel:** The residuals, green, when removing an MSF
model of the 21-cm foreground from the data showing a clear remnant of
the signal, red. When jointly fitting an MSF and signal model we find that
we can accurately recover the signal itself (see @Bevins).](example.png)

``maxsmooth`` is applicable to any experiment in which the signal of interest
has to be separated
from comparatively high magnitude smooth signals or foregrounds.

# ``maxsmooth``

DCFs can be fitted with routines such as Basin-hopping [@Basinhopping] and
Nelder-Mead [@Nelder-Mead] and this has
been the practice for 21-cm cosmology [@MSFCD; @MSF-EDGES].
However, ``maxsmooth`` employs quadratic programming via
[CVXOPT](https://pypi.org/project/cvxopt/) to
rapidly and efficiently fit DCFs which are constrained such that

$$  \frac{d^my}{dx^m}\geq0 ~~ \textnormal{or} ~~ \frac{d^my}{dx^m}~\leq0. $$

An example DCF from the ``maxsmooth`` library is given by

$$ y ~ = ~ \sum_{k=0}^{N} ~ a_{k} ~ x^k, $$

where $x$ and $y$ are the independent and dependent variables respectively and $N$
is the order of the fit with powers from $0 - (N-1)$. The library is intended
be extended by future contributions from users.
We find that the use of quadratic programming makes ``maxsmooth``
approximately two orders of magnitude quicker than a Basin-hopping/Nelder-Mead approach.

``maxsmooth`` rephrases the above condition such that

$$ \pm_m ~ \frac{d^my}{dx^m} ~ \leq0, $$

where the $\pm$ applies to a given $m$. This produces a set of sign spaces
with different combinations of constrained positive and negative derivatives. In each sign space
the associated minimum is found using quadratic programming and then ``maxsmooth``
identifies the optimum combination of signs, $\mathbf{s}$. To
summarise the minimisation problem we have

$$\min_{a,~s}~~\frac{1}{2}~\mathbf{a}^T~\mathbf{Q}~\mathbf{a}~+~\mathbf{q}^T~\mathbf{a},$$
$$\mathrm{s.t.}~~\mathbf{G(s)~a} \leq \mathbf{0},$$

where we are minimising $\chi^2$, $\mathbf{G(s)a}$ is a stacked matrix of derivative evaluations and $\mathbf{a}$
is the matrix of parameters we are optimising for. $\mathbf{Q}$ and $\mathbf{q}$
are given by

$$\mathbf{Q}~=~ \mathbf{\Phi}^T~\mathbf{\Phi}~~\textnormal{and}~~ \mathbf{q}^T~=~-\mathbf{y}^T~\mathbf{\Phi},$$

here $\mathbf{\Phi}$ is a matrix of basis function evaluations and $\mathbf{y}$
is a column matrix of the dependent data.

The discrete spaces can be searched in their entirety quickly and efficiently or
a sign navigating algorithm can be invoked using ``maxsmooth``
reducing the fitting time. Division of the parameter space into
sign spaces allows for a more complete exploration
when compared to Basin-hopping/Nelder-Mead based algorithms.

The sign navigating approach uses a cascading algorithm to identify a candidate
optimum $\mathbf{s}$ and $\mathbf{a}$. The algorithm starts with a randomly generated $\mathbf{s}$. Each
individual sign is then flipped, from the lowest order derivative first, until the
objective function decreases in value. The signs associated with the lower
$\chi^2$ value become the optimum set and the process is repeated until
$\chi^2$ stops decreasing. This is followed by a limited exploration
of the neighbouring sign spaces to identify the true global minimum.

![The time taken to fit polynomial data following an approximate $x^a$ power law
using both ``maxsmooth`` quadratic programming methods and for comparison a method
based in Basin-hopping and Nelder-Mead routines. We show the results using the later method
up to $N = 7$ after which the method begins to fail without adjustments to the routine parameters.
For $N = 3 - 7$ we find a maximum difference of $0.04\%$ between the optimum ``maxsmooth`` $\chi^2$
values and the Basin-hopping results. Figure taken from @Bevins.](times.png)

Documentation for ``maxsmooth`` is available at [ReadTheDocs](maxsmooth.readthedocs.io/)
and the code can be
found on [Github](https://github.com/htjb/maxsmooth). The code is also pip installable
([PyPI](https://pypi.org/project/maxsmooth/)). Continuous
integration is performed with [Travis](https://travis-ci.com/github/htjb/maxsmooth)
and [CircleCi](https://circleci.com/gh/htjb/maxsmooth). The
associated code coverage can be found at [CodeCov](https://codecov.io/gh/htjb/maxsmooth).

# Acknowledgements

Discussions on the applications of the software were provided by Eloy de Lera Acedo,
Will Handley and Anastasia Fialkov. The author is supported by the Science and
Technology Facilities Council (STFC) via grant number ST/T505997/1.

# References
==================================================
maxsmooth: Derivative Constrained Function Fitting
==================================================

Introduction
------------

:maxsmooth: Derivative Constrained Function Fitting
:Author: Harry Thomas Jones Bevins
:Version: 1.2.1
:Homepage: https://github.com/htjb/maxsmooth
:Documentation: https://maxsmooth.readthedocs.io/

.. image:: https://github.com/htjb/maxsmooth/workflows/CI/badge.svg?event=push
   :target: https://github.com/htjb/maxsmooth/actions
   :alt: github CI
.. image:: https://codecov.io/gh/htjb/maxsmooth/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/htjb/maxsmooth
   :alt: Test Coverage Status
.. image:: https://readthedocs.org/projects/maxsmooth/badge/?version=latest
   :target: https://maxsmooth.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status
.. image:: https://img.shields.io/badge/license-MIT-blue.svg
   :target: https://github.com/htjb/maxsmooth/blob/master/LICENSE
   :alt: License information
.. image:: https://pypip.in/v/maxsmooth/badge.svg
   :target: https://pypi.org/project/maxsmooth/#description
   :alt: Latest PyPI version
.. image:: https://img.shields.io/badge/ascl-2008.018-blue.svg?colorB=262255
   :target: http://ascl.net/2008.018
   :alt: Astrophysics Source Code Library
.. image:: https://joss.theoj.org/papers/7f53a67e2a3e8f021d4324de96fb59c8/status.svg
   :target: https://joss.theoj.org/papers/7f53a67e2a3e8f021d4324de96fb59c8
   :alt: JOSS paper
.. image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/htjb/maxsmooth/master?filepath=example_notebooks%2F
.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.4059339.svg
   :target: https://doi.org/10.5281/zenodo.4059339

Installation
~~~~~~~~~~~~
In the following two sections we highlight the purpose of ``maxsmooth`` and
show an example. To install the software follow these instructions:

The software can be pip installed from the PYPI repository like so,

.. code::

 pip install maxsmooth

or alternatively it can be installed from the git repository via,

.. code::

 git clone https://github.com/htjb/maxsmooth.git
 cd maxsmooth
 python setup.py install --user

Derivative Constrained Functions and ``maxsmooth``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``maxsmooth`` is an open source software, written in Python (supporting version 3 upwards),
for fitting derivative constrained
functions (DCFs) such as Maximally Smooth Functions
(MSFs) to data sets. MSFs are functions for which there are no zero
crossings in derivatives of order m >= 2 within the domain of interest.
More generally for DCFs the minimum
constrained derivative order, m can take on any value or a set of
specific high order derivatives can be constrained.
They are designed to prevent the loss of
signals when fitting out dominant smooth foregrounds or large magnitude signals that
mask signals of interest. Here "smooth" means that the foregrounds follow power
law structures in the band of interest.
In some cases DCFs can be used to
highlight systematics in the data.

``maxsmooth`` uses quadratic programming implemented with ``CVXOPT`` to fit
data subject to a fixed linear constraint, Ga <= 0, where the product
Ga is a matrix of derivatives.
The constraint on an MSF are not explicitly
linear and each constrained derivative can be positive or negative.
``maxsmooth`` is, however, designed to test the <= 0 constraint multiplied
by a positive or negative sign. Where a positive sign in front of the m\ :sup:`th`
order derivative forces the derivative
to be negative for all x. For an N\ :sup:`th` order polynomial ``maxsmooth`` can test
every available sign combination but by default it implements a sign navigating algorithm.
This is detailed in the ``maxsmooth`` paper (see citation), is summarized
below and in the software documentation.

The available sign combinations act as discrete parameter spaces all with
global minima and ``maxsmooth`` is capable of finding the minimum of these global
minima by implementing a cascading algorithm which is followed by a directional
exploration. The cascading routine typically finds an approximate to the global
minimum and then the directional exploration is a complete search
of the sign combinations in the neighbourhood
of that minimum. The searched region is limited by factors
that encapsulate enough of the neighbourhood to confidently return the global minimum.

The sign navigating method is reliant on the problem being "well defined" but this
is not always the case and it is in these instances it is possible to run the code testing
every available sign combination on the constrained derivatives. For a definition of
a "well defined" problem and it's counter part see the ``maxsmooth`` paper and the
documentation.

``maxsmooth`` features a built in library of DCFs or
allows the user to define their own. The addition of possible inflection points
and zero crossings in higher order derivatives is also available to the user.
The software has been designed with these two
applications in mind and is a simple interface.

Example Fit
~~~~~~~~~~~

Shown below is an example MSF fit performed with ``maxsmooth`` to data that
follows a y = x\ :sup:`-2.5` power law with a randomly generated Gaussian
noise with a standard deviation 0.02. The top panel shows the data and the
bottom panel shows the residual
after subtraction of the MSF fit alongside the actual noise in the data.
The software using the default built-in DCF model is shown to be
capable of recovering the random noise.

.. image:: https://github.com/htjb/maxsmooth/raw/master/docs/images/README.png
  :width: 400
  :align: center

Further examples can be found in the Documentation (https://maxsmooth.readthedocs.io/)
and in the github repository in the files 'example_codes/' and
'example_notebooks/' (notebooks can also be accessed online
`here <https://mybinder.org/v2/gh/htjb/maxsmooth/master?filepath=example_notebooks%2F>`__).

Licence and Citation
~~~~~~~~~~~~~~~~~~~~

The software is free to use on the MIT open source license. However if you use
the software for academic purposes we request that you cite the ``maxsmooth``
papers. They are detailed below.

MNRAS paper (referred to throughout the documentation as the ``maxsmooth``
paper),

  H. T. J. Bevins et al., `maxsmooth: Rapid maximally smooth function fitting with
  applications in Global 21-cm cosmology <https://academic.oup.com/mnras/advance-article/doi/10.1093/mnras/stab152/6105349>`__,
  Monthly Notices of the Royal Astronomical Society, 2021;, stab152, https://doi.org/10.1093/mnras/stab152

Below is the BibTex citation,

.. code:: bibtex

  @article{10.1093/mnras/stab152,
    author = {Bevins, H T J and Handley, W J and Fialkov, A and Acedo, E de Lera and Greenhill, L J and Price, D C},
    title = "{maxsmooth: rapid maximally smooth function fitting with applications in Global 21-cm cosmology}",
    journal = {Monthly Notices of the Royal Astronomical Society},
    year = {2021},
    month = {01},
    issn = {0035-8711},
    doi = {10.1093/mnras/stab152},
    url = {https://doi.org/10.1093/mnras/stab152},
    note = {stab152},
    eprint = {https://academic.oup.com/mnras/advance-article-pdf/doi/10.1093/mnras/stab152/35931358/stab152.pdf},
  }

JOSS paper,

  Bevins, H. T., (2020). maxsmooth: Derivative Constrained Function Fitting. Journal of Open Source Software, 5(54), 2596, https://doi.org/10.21105/joss.02596

and the BibTex,

.. code:: bibtex

  @article{Bevins2020,
      doi = {10.21105/joss.02596},
      url = {https://doi.org/10.21105/joss.02596},
      year = {2020},
      publisher = {The Open Journal},
      volume = {5},
      number = {54},
      pages = {2596},
      author = {Harry T. j. Bevins},
      title = {maxsmooth: Derivative Constrained Function Fitting},
      journal = {Journal of Open Source Software}
  }


Contributing
~~~~~~~~~~~~

Contributions to ``maxsmooth`` are welcome and can be made via:

- Opening an issue to purpose new features/report bugs.
- Making a pull request. Please consider opening an issue to discuss
  any proposals beforehand and ensure that your PR will be accepted.

An example contribution may be the addition of a basis function into the
standard library.

Documentation
~~~~~~~~~~~~~
The documentation is available at: https://maxsmooth.readthedocs.io/

Alternatively, it can be compiled locally from the git repository and requires
`sphinx <https://pypi.org/project/Sphinx/>`__ to be installed.
You can do this via:

.. code::

  cd docs/
  make SOURCEDIR=source html

or

.. code::

  cd docs/
  make SOURCEDIR=source latexpdf

The resultant docs can be found in the docs/_build/html/ and docs/_build/latex/
respectively.

Requirements
~~~~~~~~~~~~

To run the code you will need the following additional packages:

- `matplotlib <https://pypi.org/project/matplotlib/>`__
- `numpy <https://pypi.org/project/numpy/>`__
- `CVXOPT <https://pypi.org/project/cvxopt/>`__
- `scipy <https://pypi.org/project/scipy/>`__
- `progressbar <https://pypi.org/project/progressbar/>`__

When installing via pip or from source using the setup.py file
the above packages will also be installed if absent.

To compile the documentation locally you will need:

- `sphinx <https://pypi.org/project/Sphinx/>`__
- `numpydoc <https://pypi.org/project/numpydoc/>`__

To run the test suit you will need:

- `pytest <https://pypi.org/project/pytest/>`__

Basin-hopping/Nelder-Mead Code
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the ``maxsmooth`` MNRAS paper and JOSS paper we provide a comparison of
``maxsmooth`` to a Basin-hopping/Nelder-Mead approach for fitting DCFs. For
completeness we provide in this repo the code used to make this comparison
in the file 'Basin-hopping_Nelder_Mead/'.

The code times_chis.py is used to call ``maxsmooth`` and the Basin-hopping
methods (in the file 'BHNM/'). It will plot the recorded times and objective
function evaluations.

The Basin-hopping/Nelder-Mead code is designed to fit MSFs and is not
generalised to all types of DCF. It is also not documented, however there are
minor comments in the script and it should be self explanatory. Questions
on this are welcome and can be posted as an issue or by contacting the author.
This example code illustrates how to define your own basis function for the
DCF model.
It implements a modified version of the built in normalized polynomial model
but the structure is the same for more elaborate models.

As always we need to import the data, define an order :math:`{N}`
and import the function fitting routine, smooth().

.. code::

  import numpy as np
  from maxsmooth.DCF import smooth

  x = np.load('Data/x.npy')
  y = np.load('Data/y.npy')

  N=10

There are several requirements needed to define a new basis function completely
for ``maxsmooth`` to be able to fit it. They are as summarized below and then
examples of each are given in more detail,

    * **args:** Additional non-standard  arguments needed in the definition of the
      basis. The standard arguments are the data (x and y), the order of the fit N,
      the pivot point about which a model can be fit,
      the derivative order :math:`{m}` and the params. While the
      pivot point is not strictly needed it is a required argument for the
      functions defining a new basis to help the user in their definition.

    * **basis_functions:** This function defines the basis of the DCF model,
      :math:`{\phi}` where the model can be generally defined as,

      .. math::

          y = \sum_{k = 0}^N a_k \phi_k(x)

      where :math:`{a_k}` are the fit parameters.

    * **model:** This is the function described by the equation above.

    * **derivative:** This function defines the :math:`{m^{th}}` order derivative.

    * **derivative_pre:** This function defines the prefactors,
      :math:`{\mathbf{G}}` on the derivatives where ``CVXOPT``, the quadratic
      programming routine used, evaluates the constraints as,

      .. math::

          \mathbf{Ga} \leq \mathbf{h}

      where :math:`{\mathbf{a}}` is the matrix of parameters and :math:`{\mathbf{h}}`
      is the matrix of constraint limits. For more details on this see the ``maxsmooth``
      paper.


We can begin defining our new basis function by defining the additional arguments
needed to fit the model as a list,

.. code::

  arguments = [x[-1]*10, y[-1]*10]

The next step is to define the basis functions :math:`{\phi}`. This needs to be
done in a function that has the arguments *(x, y, pivot_point, N, \*args)*. 'args'
is optional but since we need them for this basis we are passing it in.

The basis functions, :math:`{\phi}`, should be an array of dimensions len(x)
by N and consequently evaluated at each N and x data point as shown below.

.. code::

  def basis_functions(x, y, pivot_point, N, *args):

      phi = np.empty([len(x), N])
      for h in range(len(x)):
          for i in range(N):
              phi[h, i] = args[1]*(x[h]/args[0])**i

      return phi

We can define the model that we are fitting in a function like that shown below.
This is used for evaluating :math:`{\chi^2}` and returning the optimum fitted model
once the code has finished running. It requires the arguments
*(x, y, pivot_point, N, params, \*args)* in that order and again where 'args' is optional.
'params' is the parameters of the fit, :math:`{\mathbf{a}}` which should have length
:math:`{N}`.

The function should return the fitted estimate of y.

.. code::

  def model(x, y, pivot_point, N, params, *args):

      y_sum = args[1]*np.sum([
          params[i]*(x/args[0])**i
          for i in range(N)], axis=0)

      return y_sum

Next we have to define a function for the derivatives of the model which
takes arguments *(m, x, y, N, pivot_point, params, *args)* where :math:`{m}` is
the derivative order. The function should return the :math:`{m^{th}}` order
derivative evaluation and is used for checking that the constraints have been
met and returning the derivatives of the optimum fit to the user.

.. code::

  def derivative(m, x, y, N, pivot_point, params, *args):

      mth_order_derivative = []
      for i in range(N):
          if i <= m - 1:
              mth_order_derivative.append([0]*len(x))
      for i in range(N - m):
              mth_order_derivative_term = args[1]*np.math.factorial(m+i) / \
                  np.math.factorial(i) * \
                  params[int(m)+i]*(x)**i / \
                  (args[0])**(i + 1)
              mth_order_derivative.append(
                  mth_order_derivative_term)

      return mth_order_derivative

Finally we have to define :math:`{\mathbf{G}}` which is used by ``CVXOPT`` to
build the derivatives and constrain the functions. It takes arguments
*(m, x, y, N, pivot_point, \*args)* and should return the prefactor on the
:math:`{m^{th}}` order derivative. For a more thorough definition of the
prefactor on the derivative and an explanation of how the problem is
constrained in quadratic programming see the ``maxsmooth`` paper.

.. code::

  def derivative_pre(m, x, y, N, pivot_point, *args):

      mth_order_derivative = []
      for i in range(N):
          if i <= m - 1:
              mth_order_derivative.append([0]*len(x))
      for i in range(N - m):
              mth_order_derivative_term = args[1]*np.math.factorial(m+i) / \
                  np.math.factorial(i) * \
                  (x)**i / \
                  (args[0])**(i + 1)
              mth_order_derivative.append(
                  mth_order_derivative_term)

      return mth_order_derivative

With our functions and additional arguments defined we can pass these
to the ``maxsmooth`` smooth() function as is shown below. This overwrites the
built in DCF model but you are still able to modify the fit type i.e. testing all
available sign combinations or sampling them.

.. code::

  result = smooth(x, y, N,
      basis_functions=basis_functions, model=model,
      derivatives=derivative, der_pres=derivative_pre, args=arguments)

The output of the fit can be accessed as before,

.. code::

  print('Objective Funtion Evaluations:\n', result.optimum_chi)
  print('RMS:\n', result.rms)
  print('Parameters:\n', result.optimum_params)
  print('Fitted y:\n', result.y_fit)
  print('Sign Combinations:\n', result.optimum_signs)
  print('Derivatives:\n', result.derivatives)
This example will walk the user through implementing DCF fits to data sets with
turning points and inflection points. It builds on the details in the
'Simple Example Code' and uses the 'constraints' keyword argument introduced
there. The 'constraints' keyword argument is used to adjust the type of DCF that
is being fitted. Recall that by default ``maxsmooth`` implements a Maximally
Smooth Function or MSF with constraints=2 i.e. derivatives of order :math:`{m \geq 2}`
constrained so that they do not cross zero. This allows for turning points in the
DCF as illustrated below.

We start by generating some noisy data that we know will include a turning point
and defining the order of the DCF we would like to fit.

.. code:: bash

    import numpy as np

    x = np.linspace(-10, 10, 100)
    noise = np.random.normal(0, 0.02, 100)
    y = x**(2) + noise

    N = 10

We can go ahead and plot this data just to double check it features a turning
point.

.. code:: bash

    import matplotlib.pyplot as plt

    plt.plot(x, y)
    plt.xlabel('x', fontsize=12)
    plt.ylabel('y', fontsize=12)
    plt.show()

.. image:: https://github.com/htjb/maxsmooth/raw/master/docs/images/turning_point_example.png
  :width: 400
  :align: center

As already stated ``maxsmooth`` does not constrain the first derivative of the
DCF by default so we can go ahead and fit the data.

.. code:: bash

    from maxsmooth.DCF import smooth

    res = smooth(x, y, N)

If we than plot the resultant residuals we will see that despite the data
having a turning point present we have recovered the Gaussian noise.

.. code:: bash

    plt.plot(x, y- res.y_fit, label='Recovered Noise')
    plt.plot(x, noise, label='Actual Noise')
    plt.ylabel(r'$\delta y$', fontsize=12)
    plt.xlabel('x', fontsize=12)
    plt.legend()
    plt.show()

.. image:: https://github.com/htjb/maxsmooth/raw/master/docs/images/turning_point_example_res.png
  :width: 400
  :align: center

To illustrate what happens when there is an inflection point in the data we can
define some sinusoidal data as so.

.. code:: bash

    x = np.linspace(1, 5, 100)
    noise = np.random.normal(0, 0.02, 100)
    y = np.sin(x) + noise

    N = 10

    plt.plot(x, y)
    plt.xlabel('x', fontsize=12)
    plt.ylabel('y', fontsize=12)
    plt.show()

.. image:: https://github.com/htjb/maxsmooth/raw/master/docs/images/inflection_point_example.png
  :width: 400
  :align: center

If we proceed to fit this with smooth() in its default settings we will get a
poor fit as by default the second derivative is constrained. We need to lift this
constraint to allow for the prominent inflection point to be modelled. We do this
by setting the keyword argument constraints=3 creating a Partially Smooth Function
or PSF.

.. code:: bash

    res_msf = smooth(x, y, N)
    res_psf = smooth(x, y, N, constraints=3)

    plt.plot(x, y, label='Data')
    plt.plot(x, res_msf.y_fit, label=r'MSF fit, $m \geq 2$')
    plt.plot(x, res_psf.y_fit, label=r'PSF fit, $m \geq 3$')
    plt.xlabel('x', fontsize=12)
    plt.ylabel('y', fontsize=12)
    plt.legend()
    plt.show()

.. image:: https://github.com/htjb/maxsmooth/raw/master/docs/images/inflection_point_example_fits.png
  :width: 400
  :align: center

Finally, we can plot the residuals to further see that by lifting the constraint on the
second derivative we have allowed an inflection point in the data.

.. code::

    plt.plot(x, y- res_psf.y_fit, label='Recovered Noise')
    plt.plot(x, noise, label='Actual Noise')
    plt.ylabel(r'$\delta y$', fontsize=12)
    plt.xlabel('x', fontsize=12)
    plt.legend()
    plt.show()

.. image:: https://github.com/htjb/maxsmooth/raw/master/docs/images/inflection_point_example_res.png
  :width: 400
  :align: center
Unreleased changes are not yet included in the pip install but are pushed to the
github.

Unrealeased
~~~~~~~~~~~

Version 1.1.0
~~~~~~~~~~~~~

- Two bug fixes in param_plotter()
- Extension of param_plotter() function to plot data, fit and residuals
  alongside the parameter space if required.
- Extension of param_plotter() to allow for highlighting of central
  regions in each panel if required.
- Inclusion of some theory into the documentation

Version 1.2.0
~~~~~~~~~~~~~

- Minor bug fix in param_plotter()
- Extension of the basis_test() function to allow users to compare different
  types of DCF not just MSFs. 
This section has been adapted from section 4 of the ``maxsmooth`` paper
in order to explain how the algorithm works. What follows is a discussion of
the fitting problem and the
``maxsmooth`` algorithm. To state concisely the problem being fitted we have

.. math::

        &\min_{a,~s}~~\frac{1}{2}~\mathbf{a}^T~\mathbf{Q}~\mathbf{a}~+~\mathbf{q}^T~\mathbf{a}, \\
        &\mathrm{s.t.}~~\mathbf{G(s)~a} \leq \mathbf{0}.

where :math:`{\mathbf{s}}` are the ``maxsmooth`` signs corresponding to the
signs on the derivatives. :math:`{\mathbf{G}}` is a matrix of prefactors on the derivatives,
:math:`{\mathbf{a}}` are the parameters we are optimising for and their
product gives the derivatives we are constraining with each fit.
:math:`{\mathbf{Q}}` is the dot product of the matrix of basis functions and
its transpose and :math:`\mathbf{q}` is the negative of the transposed data,
:math:`\mathbf{y}` dotted with the basis functions. For more details on this
equation see the ``maxsmooth`` paper.
A `problem' in this context is the combination of the data, order, basis
function and constraints on the DCF.

With ``maxsmooth`` we can test all possible sign combinations on the constrained derivatives.
This is a
reliable method and, provided the problem can be solved with quadratic programming,
will always give the correct global minimum. When the problem we are interested
in is "well defined", we can develop a quicker algorithm that searches or navigates
through the discrete ``maxsmooth`` sign spaces to find the global minimum.
Each sign space is a discrete parameter space with its own global minimum.
Using quadratic programming on a fit with a specific sign combination will
find this global minimum, and we are interested in finding the minimum
of these global minima.

A "well defined" problem is one in which the discrete sign spaces have large
variance in their minimum :math:`{\chi^2}` values and the sign space for the
global minimum is easily identifiable. In contrast we can have an "ill defined"
problem in which the variance in minimum :math:`{\chi^2}` across all sign
combinations is small. This concept of "well defined" and "ill defined" problems
is explored further in the following two sections.

Well Defined Problems and Discrete Sign Space Searches
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The :math:`{\chi^2}` Distribution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We investigate the distribution of :math:`{\chi^2}` values, shown in the figure below,
for a 10 :math:`{^{th}}` order y-log(x) space MSF fit to a :math:`{y = x^{-2.5}}`
power law plus gaussian noise.

In the figure, a combination of all positive derivatives~(negative signs) and
all negative derivatives~(positive signs) corresponds to sign combination numbers
255 and 0 respectively. Specifically, the ``maxsmooth`` signs, :math:`{\mathbf{s}}`,
are related to the sign combination number by its :math:`{C}` bit binary representation,
here :math:`{C = (N -2)}`. In binary the sign combination numbers run from
00000000 to 11111111. Each bit represents the sign on the :math:`{m^{th}}`
order derivative with a 1 representing a negative ``maxsmooth`` sign.

.. image:: https://github.com/htjb/maxsmooth/raw/master/docs/images/chi_dist_theory.png
  :width: 400
  :align: center

The distribution appears to be composed of smooth steps or shelves; however,
when each shelf is studied closer, we find a series of peaks and troughs. This can
be seen in the subplot of the above figure which shows the distribution in the
neighbourhood of the global minimum found in the large or `global' well. This type
of distribution with a large variance in :math:`{\chi^2}` is characteristic of a "well defined"
problem. We use this example :math:`{\chi^2}` distribution to motivate the ``maxsmooth``
algorithm outlined in the following section.

The ``maxsmooth`` Sign Navigating Algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Exploration of the discrete sign spaces for high :math:`{N}` can be achieved by
exploring the spaces around an iteratively updated optimum sign combination.
The ``maxsmooth`` algorithm begins with a randomly generated set of signs for
which the objective function is evaluated and the optimum parameters are found.
We flip each individual sign one at a time beginning with the lowest order
constrained derivative first. When the objective function is evaluated to be lower
than that for the optimum sign combination, we replace it with the new set and repeat
the process in a `cascading' routine until the objective function stops decreasing in value.

The local minima shown in the :math:`{\chi^2}` distribution above mean that the
cascading algorithm is not sufficient to consistently find the global minimum.
We can demonstrate this by performing 100 separate runs of the cascading
algorithm on :math:`{y = x^{-2.5} + \mathrm{noise}}`, and we use a y-log(x) space
:math:`{10^{th}}` order MSF again. We find the true global minimum 79
times and a second local minimum 21 times.

To prevent the routine terminating in a local minimum we perform a complete search
of the sign spaces surrounding the minimum found after the cascading routine.
We refer to this search as a directional exploration and impose limits on its
extent. In each direction we limit the number of sign combinations to explore and
we limit the maximum allowed increase in :math:`{\chi^2}` value. These limits can
be modified by the user. We prevent repeated calculations of the minimum for given
signs and treat the minimum of all tested signs as the global minimum.

We run the consistency test again, with the full ``maxsmooth`` algorithm, and find
that for all 100 trial fits we find the same :math:`{\chi^2}` found when testing
all sign combinations. In the figure below, the red arrows show the approximate path
taken through the discrete sign spaces against the complete distribution of :math:`{\chi^2}`.
Point (1a) shows the random starting point in the algorithm, and point (1b) shows a rejected sign
combination evaluated during the cascade from point (1a) to (2). Point (2), therefore,
corresponds to a step through the cascade. Point (3) marks the end of the cascade
and the start of the left directional exploration. Finally, point (4) shows the end
of the right directional exploration where the calculated :math:`{\chi^2}`
value exceeds the limit on the directional exploration.

.. image:: https://github.com/htjb/maxsmooth/raw/master/docs/images/routine.png
  :width: 400
  :align: center

The global well tends to be associated with signs that are all positive,
all negative or alternating. We see this in the figure above where the minimum falls
at sign combination number 169 and number 170, characteristic of the derivatives for
a :math:`{x^{-2.5}}` power law, corresponds to alternating positive and negative
derivatives from order :math:`{m = 2}`. Standard patterns of derivative signs can be seen
for all data following approximate power laws. All positive derivatives, all negative
and alternating signs correspond to data following the approximate power laws
:math:`{y\approx x^{k}}`, :math:`{y\approx -x^{k}}`, :math:`{y\approx x^{-k}}` and
:math:`{y\approx -x^{-k}}`.

The ``maxsmooth`` algorithm assumes that the global well is present in the :math:`{\chi^2}`
distribution and this is often the case. The use of DCFs is primarily driven by a
desire to constrain previously proposed polynomial models to foregrounds. As a result
we would expect that the data being fitted could be described by one of the four
approximate power laws highlighted above and that the global minimum will fall
around an associated sign combination. In rare cases the global well is not clearly
defined and this is described in the following subsection.

Ill Defined Problems and their Identification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We can illustrate an "ill defined" problem, with a small variation in
:math:`{\chi^2}` across the ``maxsmooth`` sign spaces, by adding a non-smooth signal
of interest into the foreground model, :math:`{x^{-2.5}}` and fitting this with
a 10 :math:`{^{th}}` order log(y)-log(x) space MSF. We add an additional noise of
:math:`{0.020}` to the mock data. The resultant :math:`{\chi^2}` distribution with its
global minimum is shown in the top panel of the figure below.

The global minimum, shown as a black data point, cannot be found using the
``maxsmooth`` algorithm. The cascading algorithm may terminate in any of the
approximately equal minima and the directional exploration will then quickly
terminate because of the limits imposed.

.. image:: https://github.com/htjb/maxsmooth/raw/master/docs/images/combined_chi.png
  :width: 400
  :align: center

If we repeat the above fit and perform it with a y-x space MSF we find that the
problem is well defined with a larger :math:`{\chi^2}` variation across sign
combinations. This is shown in the bottom panel of the above figure. The results,
when using the log(y)-log(x) space MSF, are significantly better than when using
y-x space MSF meaning it is important to be able to solve "ill defined" problems.
This can be done by testing all ``maxsmooth`` signs but knowing when this is
necessary is important if you are expecting to run multiple DCF fits to the
same data set. We can focus on diagnosing whether a DCF fit to the data is
"ill defined" because a joint fit to the same data set of a DCF and signal
of interest will also feature an "ill defined" :math:`{\chi^2}` distribution.

We can identify an "ill defined" problem by producing the equivalent of
the above figure using ``maxsmooth`` and visually assessing the :math:`{\chi^2}`
distribution for a DCF fit. Alternatively, we can use the parameter space plots,
detailed in the ``maxsmooth`` paper and later in this documentation,
to identify whether the constraints are weak or not, and if a local minima is
returned from the sign navigating routine then the minimum in these plots
will appear off centre.

Assessment of the first derivative of the data can also help to identify an
"ill defined" problem. For the example problem this is shown in the figure below
where the derivatives have been approximated using :math:`{\Delta y/ \Delta x}`.
Higher order derivatives of the data will have similarly complex or simplistic
structures in the respective spaces. There are many combinations of parameters
that will provide smooth fits with similar :math:`{\chi^2}` values in logarithmic
space leading to the presence of local minima. This issue will also be present
in any data set where the noise or signal of interest are of a similar magnitude
to the foreground in y - x space.

.. image:: https://github.com/htjb/maxsmooth/raw/master/docs/images/Gradients_fits.png
  :width: 400
  :align: center
This function can be used to identify which of the built in DCFs
fits the data best before running joint fits.

To use it we begin by loading in the data,

.. code::

  import numpy as np

  x = np.load('Data/x.npy')
  y = np.load('Data/y.npy')

and then importing the basis_test() function.

.. code::

  from maxsmooth.best_basis import basis_test

To call the function we use,

.. code::

  basis_test(x, y, base_dir='examples/', N=np.arange(3, 16, 1))

The function only requires the data but we can provide it with a base directory,
fit type and range of DCF orders to test. By defualt it uses the sign navigating
algorithm and tests :math:`{N = 3 - 13}`. Here we test the range
:math:``{N = 3 - 15}``.
The resultant graph is saved in the
base directory and the example generated here is shown below.

.. image:: https://github.com/htjb/maxsmooth/raw/master/docs/images/Basis_functions.png
  :width: 400
  :align: center

The graph shows us which basis is the optimum for solving this problem from the
built in library (that which can reach the minimum :math:``{\chi^2}``). If we
were to go to higher N we would also find that the :math:``{\chi^2}`` value
would stop decreasing in value. The value of N for which this occurs at is the
optimum DCF order. (See the ``maxsmooth`` paper for a real world application
of this concept.)

We can also provide this function with additional arguments such as the
fit type, minimum constrained derivative, directional exploration limits
ect. (see the ``maxsmooth`` Functions section).
.. highlight:: python

In order to run the ``maxsmooth`` software using the built
in DCF models for a simple fit the user can follow the simple structure detailed here.

An important point to make is that by default ``maxsmooth`` fits a
Maximally Smooth Function or MSF to the data. An MSF, as stated in
the introduction to the documentation, is a function which has
derivatives of order :math:`{m \geq 2}` constrained so that they do not cross
0. This means that they do not have inflection points or non smooth
structure produced by higher order derivatives. More generally a DCF
follows the constraint,

.. math:

    \frac{\delta^m y}{\delta x^m} \leq 0 ~~\mathrm{or}~~ \frac{\delta^m y}{\delta x^m} \geq 0 $

for every constrained order :math:`{m}`. The set of :math:`{m}` can be any set of
derivative orders as long as those derivatives exist for the function.

This means we can use ``maxsmooth`` to produce different DCF
models. MSFs are one of two special cases of DCF and we can also
have a Completely Smooth Function (CSF) with orders :math:`{m \geq 1}`
constrained. Alternatively we can have Partially Smooth Functions
(PSF) which are much more general and can have arbitrary sets of
derivatives constrained. We illustrate how this is implemented
towards the end of this example but we begin with the default case
fitting a MSF.

The user should begin by importing the `smooth` class from `maxsmooth.DCF`.

.. code::

    from maxsmooth.DCF import smooth

The user should then import the data they wish to fit.

.. code::

    import numpy as np

    x = np.load('Data/x.npy')
    y = np.load('Data/y.npy') + np.random.normal(0, 0.02, len(x))

and define the polynomial orders they wish to fit.

.. code::

    N = [3, 4, 5, 6, 7, 8, 9, 10, 11]
    for i in range(len(N)):
        `act on N[i]`

or for example,

.. code::

    N = 15

We can also plot the data to illustrate what is happening.
Here the data is a scaled :math:`{x^{-2.5}}` power law and I have added gaussian
noise in with a standard deviation of 0.02.

.. code:: bash

    import matplotlib.pyplot as plt

    plt.plot(x, y)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()

.. image:: https://github.com/htjb/maxsmooth/raw/master/docs/images/simple_program_data.png
  :width: 400
  :align: center

`smooth` can be called as is shown below. It takes the x and y data as standard
inputs as well as the order of the fit. There are a set of keyword arguments
also available that change the type of function being fitted and these are
detailed in the documentation.

.. code::

    result = smooth(x, y, N)

and it's resulting attributes can be accessed by writing
:code:`result.attribute_name`. For example printing the outputs is done like
so,

.. code::

    print('Objective Funtion Evaluations:\n', result.optimum_chi)
    print('RMS:\n', result.rms)
    print('Parameters:\n', result.optimum_params)
    print('Fitted y:\n', result.y_fit)
    print('Sign Combinations:\n', result.optimum_signs)
    print('Derivatives:\n', result.derivatives)

    plt.plot(x, y - result.y_fit)
    plt.xlabel('x', fontsize=12)
    plt.ylabel(r'$\delta y$', fontsize=12)
    plt.tight_layout()
    plt.show()

.. image:: https://github.com/htjb/maxsmooth/raw/master/docs/images/simple_program_msf_residuals.png
  :width: 400
  :align: center

To fit the data with a CSF we can use the 'constraints' keyword
argument in smooth(). 'constraints' sets the minimum constrained
derivative for the function which for a CSF we want to be one.

.. code:: bash

    res = smooth(x, y, N, constraints=1)

Note in the printed results the number of constrained derivatives has
increased by 1 and the only derivative that is allowed to cross through 0
(Zero Crossings Used?) is the the :math:`{0^{th}}` order i.e. the data.

.. code:: bash

    plt.plot(x, y - res.y_fit)
    plt.xlabel('x', fontsize=12)
    plt.ylabel(r'$\delta y$', fontsize=12)
    plt.tight_layout()
    plt.show()

.. image:: https://github.com/htjb/maxsmooth/raw/master/docs/images/simple_program_csf_residuals.png
  :width: 400
  :align: center

A Partially Smooth Function can have derivatives constrained via :math:`{m \geq a}`
where :math:`{a}` is
any order above 2 or it can have a set of derivatives that are allowed to cross
zero. For the first case we can once again use the 'constraints' keyword
argument. For example we can constrain derivatives with orders :math:`{\geq 3}` which will
allow the :math:`{1^{st}}` and :math:`{2^{nd}}` order derivatives to cross zero.
This is useful when our
data features an inflection point we want to model with our fit.

.. code:: bash

   res = smooth(x, y, N, constraints=3)

   plt.plot(x, y - res.y_fit)
   plt.xlabel('x', fontsize=12)
   plt.ylabel(r'$\delta y$', fontsize=12)
   plt.tight_layout()
   plt.show()

.. image:: https://github.com/htjb/maxsmooth/raw/master/docs/images/simple_program_psf1_residuals.png
  :width: 400
  :align: center

To allow a particular set of derivatives to cross zero we use the
'zero_crossings' keyword. In the example below we are lifting the constraints
on the :math:`{3^{rd}}`, :math:`{4^{th}}` and :math:`{5^{th}}` order derivatives
but our minimum constrained derivative is still set at the default 2. Therefore
this PSF has derivatives of order :math:`{m = [2, 6, 7, 8, 9]}`
constrained via the condition at the begining of this example code.

.. code::

   res = smooth(x, y, N, zero_crossings=[3, 4, 5])

   plt.plot(x, y - res.y_fit)
   plt.xlabel('x', fontsize=12)
   plt.ylabel(r'$\delta y$', fontsize=12)
   plt.tight_layout()
   plt.show()

.. image:: https://github.com/htjb/maxsmooth/raw/master/docs/images/simple_program_psf2_residuals.png
  :width: 400
  :align: center

While PSFs can seem like an attractive way to improve the quality of fit they
are less 'smooth' than a MSF or CSF and consequently they can introduce
additional turning points in to your residuals obscuring any signals of
intrest.
.. toctree::
   :maxdepth: 6

``maxsmooth`` Theory and Algorithm
----------------------------------

.. include:: theory.rst

``maxsmooth`` Example Codes
---------------------------

This section is designed to introduce the user to the software and the form
in which it is run. It provides basic examples of data fitting with a built in
MSF model and a user defined model.

There are also examples of functions that can be used pre-fitting and post-fitting
for various purposes including; determination of the best DCF model from the
built in library for the problem being fitted, analysis of the :math:`{\chi^2}`
distribution as a function of the discrete sign spaces and analysis of the
parameter space surrounding the optimum results.

The data used for all of this examples is available
`here <https://github.com/htjb/maxsmooth/tree/master/example_codes/Data>`__.

The example codes can be found
`here <https://github.com/htjb/maxsmooth/tree/master/example_codes>`__ and
corresponding Jupyter Notebooks are provided
`here <https://mybinder.org/v2/gh/htjb/maxsmooth/master?filepath=example_notebooks%2F>`__.

Simple Example code
~~~~~~~~~~~~~~~~~~~
.. include:: simple_program.rst

Turning Points and Inflection Points
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. include:: turning_points.rst

New Basis Example
~~~~~~~~~~~~~~~~~

.. include:: new_basis_example.rst

Best Basis Example
~~~~~~~~~~~~~~~~~~

.. include:: best_basis.rst

:math:`{\chi^2}` Distribution Example
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. include:: chi_dist_example.rst

Parameter Plotter Example
~~~~~~~~~~~~~~~~~~~~~~~~~

.. include:: param_plotter_example.rst

``maxsmooth`` Functions
-----------------------

This section details the specifics of the built in functions in ``maxsmooth`` including
the relevant keyword arguments and default parameters for all. Where keyword arguments
are essential for the functions to run this is stated.

smooth()
~~~~~~~~

.. automodule:: maxsmooth.DCF
   :members: smooth

best_basis()
~~~~~~~~~~~~

.. automodule:: maxsmooth.best_basis
  :members: basis_test

chidist_plotter()
~~~~~~~~~~~~~~~~~

.. automodule:: maxsmooth.chidist_plotter
   :members: chi_plotter

parameter_plotter()
~~~~~~~~~~~~~~~~~~~

.. automodule:: maxsmooth.parameter_plotter
  :members: param_plotter

Change Log
----------

.. include:: CHANGELOG.rst
This example will show you how to generate a plot of the :math:`{\chi^2}`
distribution as a function of the discrete sign combinations on the constrained
derivatives.

First you will need to import your data and fit this using ``maxsmooth`` as
was done in the simple example code.

.. code::

  import numpy as np

  x = np.load('Data/x.npy')
  y = np.load('Data/y.npy')

  from maxsmooth.DCF import smooth

  N = 10
  result = smooth(x, y, N, base_dir='examples/',
    data_save=True, fit_type='qp')

Here we have used some additional keyword arguments for the 'smooth' fitting
function. 'data_save' ensures that the files containing the tested sign combinations
and the corresponding objective function evaluations exist in the base directory
which we have changed to 'base_dir='examples/''. These files are essential for
the plotting the :math:`{\chi^2}` distribution and are not saved by ``maxsmooth``
without 'data_save=True'. We have also set the 'fit_type' to 'qp' rather than the
default 'qp-sign_flipping'. This ensures that all of the available sign
combinations are tested rather than a sampled set giving us a full picture of the
distribution when we plot it. We have used the default DCF model to fit this data.

We can import the 'chi_plotter' like so,

.. code::

  from maxsmooth.chidist_plotter import chi_plotter

and produce the fit which gets placed in the base directory with the following
code,

.. code::

  chi_plotter(N, base_dir='examples/', fit_type='qp')

We pass the same 'base_dir' as before so that the plotter can find the correct output
files. We also give the function the same 'fit_type' used for the fitting which
ensures that the files can be read.

The resultant plot is shown below and the yellow star shows the global minimum.
This can be used to determine how well
the sign sampling approach using a descent and directional exploration
can find the global minimum. If the distribution looks like noise then it is
unlikely the sign sampling algorithm will consistently find the global minimum.
Rather it will likely repeatedly return the local minima found after the descent
algorithm and you should use the 'qp' method testing all available sign combinations
in any future fits to the data with this DCF model.

.. image:: https://github.com/htjb/maxsmooth/raw/master/docs/images/chi_distribution.png
  :width: 400
  :align: center
Welcome to the ``maxsmooth`` documentation!
===========================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   Introduction <intro>
   maxsmooth <maxsmooth>
.. title:: Introduction

.. include:: ../../README.rst
We can assess the parameter space around the optimum solution
found using ``maxsmooth`` with the param_plotter() function.
This can help us identify how well a problem can be solved using the
sign navigating approach employed by ``maxsmooth`` or simply
be used to identify correlations between the foreground parameters.
For more details on this see the ``maxsmooth`` paper.

We begin by importing and fitting the data as with the chi_plotter()
function illustrated above.

.. code::

  import numpy as np

  x = np.load('Data/x.npy')
  y = np.load('Data/y.npy')

  from maxsmooth.DCF import smooth

  N = 5
  result = smooth(x, y, N, base_dir='examples/', fit_type='qp')

We have changed the order of the fit to 5 to illustrate that for
order :math:`{N \leq 5}` and fits with derivatives :math:`{m \geq 2}` constrained
the function will plot each region of the graph corresponding to
different sign combinations in a different colourmap. Recall that
by default the function smooth() fits a maximally smooth function (MSF) with
derivatives of order :math:`{m \geq 2}`. If the constraints are
different or the order is greater than 5 then the viable regions will have
a single colourmap. Invalid regions are plotted as black shaded colourmaps
and the contour lines are contours of :math:`{\chi^2}`.

Specifically, invalid regions violate the condition

.. math::

  \pm_m \frac{\delta^m y}{\delta x^m} \leq 0

where :math:`{m}` represents the derivative order, :math:`{y}` is the dependent
variable and :math:`{x}` is the independent variable. Violation of the
condition means that one or more of the constrained derivatives crosses 0 in the
band of interest. For an MSF, as mentioned, :math:`{m \geq 2}` and the sign :math:`{\pm_m}`
applies to specific derivative orders. For this specific example there are
3 constrained derivatives, :math:`{m = 2, 3, 4}` and consequently 3 signs to
optimise for alongside the parameters :math:`{a_k}`. The coloured valid regions
therefore correspond to a specific combination of :math:`{\pm_m}` for the problem.
:math:`{\pm_m}` is also referred to as :math:`{\mathbf{s}}` in the theory
section and the ``maxsmooth`` paper.

We can import the function like so,

.. code::

  from maxsmooth.parameter_plotter import param_plotter

and access it using,

.. code::

  param_plotter(result.optimum_params, result.optimum_signs,
      x, y, N, base_dir='examples/')

The function takes in the optimum parameters and signs found after the fit
as well as the data and order of the fit. There are a number of keyword arguments
detailed in the following section and the resultant fit is shown below. The
function by default samples the parameter ranges 50% either side of the optimum
and calculates 50 samples for each parameter. In each panel the two
labelled parameters are varied while the others are maintained at their optimum
values.

.. image:: https://github.com/htjb/maxsmooth/raw/master/docs/images/Parameter_plot.png

We are also able to plot the data, fit and residuals alongside the parameter
plot and this can be done by setting data_plot=True. We can also highlight the
central region in each panel of the parameter space by setting center_plot=True.

.. code::

  param_plotter(result.optimum_params, result.optimum_signs,
      x, y, N, base_dir='examples/', data_plot=True, center_plot=True)

which gives us the graph below.

.. image:: https://github.com/htjb/maxsmooth/raw/master/docs/images/Parameter_plot_extended.png
