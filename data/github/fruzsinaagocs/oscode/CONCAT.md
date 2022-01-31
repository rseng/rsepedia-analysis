---
title: '(py)oscode: fast solutions of oscillatory ODEs'
tags:
  - Python
  - C++
  - numerical methods
  - ordinary differential equations
authors:
  - name: Fruzsina Julia Agocs
    orcid: 0000-0002-1763-5884
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
affiliations:
 - name: Astrophysics Group, Cavendish Laboratory, J. J. Thomson Avenue, Cambridge, CB3 0HE, UK
   index: 1
 - name: Kavli Institute for Cosmology, Madingley Road, Cambridge, CB3 0HA, UK
   index: 2
date: 7 October 2020
bibliography: paper.bib

---

# Summary

Oscillatory differential equations are ubiquitous in physics, chemistry and beyond. They arise in
quantum mechanics, electrical circuitry, suspension systems, molecular dynamics,
and in models of gravitational and electromagnetic waves.
The numerical solution of such systems however can be a computational bottleneck when tackled with conventional methods
available from numerical libraries. 

We present `(py)oscode`, a general-purpose numerical routine for solving a class of highly
oscillatory ordinary differential equations (ODEs) efficiently. The package has
been designed to solve equations which describe a single harmonic oscillator
with a time-dependent frequency and damping term, i.e. are of the form
\begin{equation}\label{eq:eom}
y'' + 2\gamma(x) y' + \omega^2(x) y = 0.
\end{equation}
The frequency $\omega(x)$ and damping $\gamma(x)$ terms do not need
to be explicit functions of $x$ (they can instead be e.g. the result of another
numerical solution of an ODE), as they are supplied as sequences $\omega_j,
\gamma_j$ evaluated at $x_i \leq x_j \leq x_f$, where $(x_i, x_f)$ is the
integration range.

`(py)oscode` is written in C++, but comes with a Python wrapper.
Its Python interface was designed to be similar to those included in `SciPy`'s [@scipy] numerical ODE solution
modules. This is demonstrated in the example below whose output is shown in
\autoref{fig:airy}.

```python
import numpy as np
import scipy.special as sp
import pyoscode

# Set up the Airy equation as an example: y'' + xy = 0
xs = np.linspace(0,40.0,5000)
ws = np.sqrt(xs)
gs = np.zeros_like(xs)
# Initial conditions
xi = 1.0
xf = 40.0
yi = sp.airy(-xi)[0]
dyi = -sp.airy(-xi)[1]
# Get dense output at the following points
t_eval = np.linspace(15,35,600)
# Solve the equation
solution = pyoscode.solve(xs, ws, gs, xi, xf, yi, dyi, t_eval=t_eval)
```

![Numerical solution of the Airy equation, $y'' + xy = 0$, with `pyoscode`. The
increase in step-size of `pyoscode`'s internal steps (orange dots) is due to the
algorithm switching from using the RK method to the WKB approximation in the presence of high-frequency
oscillations. The orange segment shows dense output, the solution at these
points was computed at no additional evaluations of terms in the differential
equation. \label{fig:airy}](../examples/images/airy.png)

# Statement of need 

Even if the terms in \autoref{eq:eom} change slowly, if the frequency of
oscillations in the solution is high enough, standard numerical methods struggle
to solve such equations quickly. Traditional methods have to trace every
oscillation in the solution, taking many steps in $x$ at an enormous
computational cost. The algorithm underlying `(py)oscode`, published in
@oscode and based on @rkwkb-handley, can detect when the solution is oscillatory and switch to a method
based on an analytic approximation (Wentzel--Kramers--Brillouin, WKB) suited for
oscillatory functions, otherwise using a Runge--Kutta (RK) method. Using the WKB
approximation allows the algorithm to skip over several wavelengths of
oscillation in a single step, reducing the number of steps taken drastically. It
adaptively updates its step-size to keep the local numerical error within a
user-specified tolerance. `(py)oscode` is capable of producing a solution estimate
at an arbitrary value of $x$, not just at its internal steps, therefore it can
be used to generate a "continuous" solution, or dense output [@dense-output]. 

# Related research and software

`(py)oscode`'s development was motivated by the need for a significantly more
efficient solver for the evolution of early-universe quantum fluctuations. These
perturbations are thought to have been stretched to macroscopic scales by a
phase of accelerated expansion of the universe (cosmic inflation), to later become the
large-scale structure we see today. To understand the origins of structure it
is therefore essential to model the perturbations and understand the physics
involved in inflation. `(py)oscode` has been used to speed up the numerical evolution of quantum
fluctuations in the early universe, enabling the exploration of models beyond
the standard model of cosmology [@pps-curved]. It served as inspiration for
other numerical methods aiming to extend the range of oscillatory ODEs to solve
[@beyond-rkwkb]. 

The efficient solution of oscillatory ODEs is a long-standing
numerical analysis problem with many existing methods to handle certain
sub-classes of equations. Examples include successful methods by Petzold [@petzold], reviewed in @petzold-review with many references therein, 
Iserles et al. [@condon-deano-iserles; @deano-integrals; @condon-et-al-circuits], and Bremer [@bremer], with code available from @bremer-code.

# Acknowledgements

I thank Lukas Hergt for invaluable discussions during the early development of
`(py)oscode` and his ongoing support. Construction of the algorithm would not have been possible
without the help and guidance of Will Handley, Mike Hobson, and Anthony Lasenby. 
I was supported by the Science and Technology Facilities Council (STFC).

# References
========================================================================
oscode: Oscillatory ordinary differential equation solver
========================================================================

.. image:: https://codecov.io/gh/fruzsinaagocs/oscode/branch/joss-paper/graph/badge.svg
        :target: https://codecov.io/gh/fruzsinaagocs/oscode
        :alt: codecov status
.. image:: https://travis-ci.org/fruzsinaagocs/oscode.svg?branch=master
        :target: https://travis-ci.org/fruzsinaagocs/oscode
        :alt: Travis CI build status
.. image:: https://readthedocs.org/projects/oscode/badge/?version=latest
        :target: https://oscode.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status
.. image:: https://badges.gitter.im/oscode-help/community.svg
        :target: https://gitter.im/oscode-help/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge
        :alt: Chat on gitter
.. image:: https://img.shields.io/badge/License-BSD%203--Clause-blue.svg
        :target: https://opensource.org/licenses/BSD-3-Clause
        :alt: BSD 3-clause license
.. image:: https://img.shields.io/pypi/dm/pyoscode?color=indigo 
        :target: https://pypi.org/project/pyoscode/
        :alt: PyPI downloads
.. image:: https://joss.theoj.org/papers/d4c9396ef9b2b595e2f3881a4f8a7cda/status.svg
        :target: https://joss.theoj.org/papers/d4c9396ef9b2b595e2f3881a4f8a7cda
        :alt: JOSS status
.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.4322958.svg
        :target: https://doi.org/10.5281/zenodo.4322958
        :alt: Zenodo doi

|
|

.. contents::
   :local:

|

About
-----

Oscode is a C++ tool with a Python interface that solves **osc**\illatory
**o**\rdinary **d**\ifferential **e**\quations efficiently. It is designed to
deal with equations of the form

.. image:: 
    https://github.com/fruzsinaagocs/oscode/raw/master/pyoscode/images/oscillator.png

where |gamma| (friction term) and |omega| (frequency) can be given as arrays.

.. |gamma| image:: https://github.com/fruzsinaagocs/oscode/raw/master/pyoscode/images/gamma.png

.. |omega| image:: https://github.com/fruzsinaagocs/oscode/raw/master/pyoscode/images/omega.png

Oscode makes use of an analytic approximation of x(t) embedded in a
stepping procedure to skip over long regions of oscillations, giving a reduction
in computing time. The approximation is valid when the frequency changes slowly
relative to the timescales of integration, it is therefore worth applying when
this condition holds for at least some part of the integration range. 

For the details of the numerical method used by oscode, see Citation_.


Installation
------------

Dependencies
~~~~~~~~~~~~

Basic requirements for using the C++ interface:

- C++11 or later
- `Eigen <http://eigen.tuxfamily.org/index.php?title=Main_Page>`__ (a header-only library included in this source)

The strictly necessary Python dependencies are automatically installed when you use `pip` or the `setup.py`. They are:

- `numpy <https://pypi.org/project/numpy/>`__

The *optional* dependencies are: 

- for tests
    - `scipy <https://pypi.org/project/scipy/>`__ 
    - `pytest <https://docs.pytest.org/en/stable/getting-started.html>`__ 
- for examples/plotting
    - `matplotlib <https://pypi.org/project/matplotlib/>`__
    - `scipy <https://pypi.org/project/scipy/>`__ 
- for generating offline documentation
    - `sphinx <https://pypi.org/project/Sphinx/>`__ 
    - `doxygen <https://www.doxygen.nl/index.html>`__
    - `breathe <https://pypi.org/project/breathe/>`__
    - `exhale <https://pypi.org/project/exhale/>`__


Python
~~~~~~

``pyoscode`` can be installed via pip 

.. code:: bash
   
   pip install pyoscode

or via the setup.py

.. code:: bash

   git clone https://github.com/fruzsinaagocs/oscode
   cd oscode
   python setup.py install --user

or

.. code:: bash

   git clone https://github.com/fruzsinaagocs/oscode
   cd oscode
   pip install .

You can then import ``pyoscode`` from anywhere. Omit the ``--user`` option if
you wish to install globally or in a virtual environment. If you have any
difficulties, check out the `FAQs - Installation
<https://github.com/fruzsinaagocs/oscode#installation-1>`__ section below. 

You can check that things are working by running `tests/` (also ran by Travis continuous integration):

.. code:: bash

   pytest tests/

C++
~~~

``oscode`` is a header-only C++ package, it requires no installation.

.. code:: bash

   git clone https://github.com/fruzsinaagocs/oscode

and then include the relevant header files in your C++ code:

.. code:: c

    #include "solver.hpp"
    #include "system.hpp"


Quick start
-----------

Try the following quick examples. They are available in the `examples
<https://github.com/fruzsinaagocs/oscode/tree/master/examples/>`__.

Python
~~~~~~

:Introduction to pyoscode: |intro_binder|
:Cosmology examples: |cosmology_binder|
:Scipy 2020 lecture notebook: |scipy_binder|

.. |intro_binder| image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/fruzsinaagocs/oscode/master?filepath=examples/introduction_to_pyoscode.ipynb

.. |cosmology_binder| image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/fruzsinaagocs/oscode/master?filepath=examples/cosmology.ipynb

.. |scipy_binder| image:: https://mybinder.org/badge_logo.svg
 :target: https://mybinder.org/v2/gh/fruzsinaagocs/oscode/master?filepath=examples/pyoscode_scipy.ipynb


.. image::
    https://github.com/fruzsinaagocs/oscode/raw/master/pyoscode/images/spectra.gif
    :width: 800

C++
~~~

:Introduction to oscode: `examples/burst.cpp`
:To plot results from `burst.cpp`: `examples/plot_burst.py`

To compile and run:

.. code:: bash

    g++ -g -Wall -std=c++11 -c -o burst.o burst.cpp
    g++ -g -Wall -std=c++11 -o burst burst.o
    ./burst


Documentation
-------------

Documentation is hosted at `readthedocs <https://oscode.readthedocs.io>`__.

To build your own local copy of the documentation you can run:

.. code:: bash

   cd pyoscode/docs
   make html

Citation
--------

If you use ``oscode`` to solve equations for a publication, please cite:

- `Efficient method for solving highly oscillatory ordinary differential equations with applications to physical systems <https://doi.org/10.1103/PhysRevResearch.2.013030>`__,
- `Dense output for highly oscillatory numerical solutions  <https://arxiv.org/abs/2007.05013>`__

Contributing
------------

Any comments and improvements to this project are welcome. You can contribute
by:

- Opening and `issue <https://www.github.com/fruzsinaagocs/oscode/issues/>`__ to report bugs and propose new features.
- Making a pull request.

Further help
------------

You can get help by submitting an issue or posting a message on `Gitter <https://gitter.im/oscode-help/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge>`__.

FAQs
----

Installation
~~~~~~~~~~~~

1. Eigen import errors:
    .. code:: bash

       pyoscode/_pyoscode.hpp:6:10: fatal error: Eigen/Dense: No such file or directory
        #include <Eigen/Dense>
                  ^~~~~~~~~~~~~

    Try explicitly including the location of your Eigen library via the
    ``CPLUS_INCLUDE_PATH`` environment variable, for example:

    .. code:: bash

       CPLUS_INCLUDE_PATH=/usr/include/eigen3 python setup.py install --user
       # or 
       CPLUS_INCLUDE_PATH=/usr/include/eigen3 pip install pyoscode

    where  ``/usr/include/eigen3`` should be replaced with your system-specific
    eigen location.

Thanks
------

Many thanks to **Will Handley**, **Lukas Hergt**, **Anthony Lasenby**, and **Mike Hobson** for
their support and advice regarding the algorithm behind `oscode`.
There are many packages without which some part of `oscode` (e.g. testing and
examples) wouldn't run as nicely and smoothly, thank you all developers for
making and maintaining these open-source projects. A special thanks goes to the
devs of `exhale <https://pypi.org/project/exhale/>`__ for making the beautiful C++ documentation possible. 


Changelog
---------

- 1.0.4:
    - set minimally required numpy version: numpy>=1.20.0
    - drop Python 2.7 support, instead support 3.8 and 3.9 in addition to 3.7
- 1.0.3: 
    - paper accepted to JOSS
- 1.0.2:
    - Fixed getting correct numpy include directories
- 1.0.1:
    - Added `pyproject.toml` to handle build dependencies (numpy)
- 1.0.0:
    - Dense output
    - Arrays for frequency and damping term need not be evenly spaced
    - Automatic C++ documentation on readthedocs
    - Eigen included in source for pip installability
    - First pip release :)
- 0.1.2:
    - Bug that occurred when beginning and end of integration coincided
      corrected
- 0.1.1:
    - Automatic detection of direction of integration
- 0.1.0:
    - Memory leaks at python interface fixed
    - C++ documentation added 
.. title:: oscode (C++ interface)

================================
Using the C++ interface (oscode)
================================

.. sectnum:: 


Overview
--------

This documentation illustrates how one can use ``oscode`` via its C++ interface.
Usage of ``oscode`` involves

- defining an equation to solve,
- solving the equation,
- and extracting the solution and other statistics about the run.

The next sections will cover each of these. For a complete reference, see the
:doc:`C++ interface reference<cpp-docs/oscode-reference>` page, and for examples
see the `examples
<https://github.com/fruzsinaagocs/oscode/tree/master/examples/>`__ directory on
GitHub.

Defining an equation
--------------------

The equations ``oscode`` can be used to solve are of the form 

.. math::

   \ddot{x}(t) + 2\gamma(t)\dot{x}(t) + \omega^2(t)x(t) = 0,

where :math:`x(t)`, :math:`\gamma(t)`, :math:`\omega(t)` can be complex. We will
call :math:`t` the independent variable, :math:`x` the dependent variable,
:math:`\omega(t)` the frequency term, and :math:`\gamma(t)` the friction or
first-derivative term. 

Defining an equation is via

- giving the frequency :math:`\omega(t)`,
- giving the first-derivative term :math:`\gamma(t)`,

Defining the frequency and the first-derivative term can either be done by
giving them as **functions explicitly**, or by giving them as **sequences** evaluated on a grid of :math:`t`.

:math:`\omega` and :math:`\gamma` as explicit functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If :math:`\omega` and :math:`\gamma` are closed-form functions of time, then
define them as

.. code:: c
    
    #include "solver.hpp" // de_system, Solution defined in here 

    std::complex<double> g(double t){
        return 0.0;
    };
    
    std::complex<double> w(double t){
        return std::pow(9999,0.5)/(1.0 + t*t);
    };

Then feed them to the solver via the de_system class:

.. code:: c
    
    de_system sys(&w, &g);   
    Solution solution(sys, ...) // other arguments left out

:math:`\omega` and :math:`\gamma` as time series
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sometimes :math:`\omega` and :math:`\gamma` will be results of numerical
integration, and they will have no closed-form functional form. In this case,
they can be specified on a grid, and ``oscode`` will perform linear
interpolation on the given grid to find their values at any timepoint. Because
of this, some important things to **note** are:

- ``oscode`` will assume the grid of timepoints :math:`\omega` and :math:`\gamma` are **not evenly spaced**. If the grids are evenly sampled, set ``even=true`` in the call for ``de_system()``, this will speed linear interpolation up significantly.
- The timepoints grid needs to be **monotonically increasing**.
- The timepoints grid needs to **include the range of integration** (:math:`t_i`,:math:`t_f`). 
- The grids for the timepoints, frequencies, and first-derivative terms have to be the **same size**.
- The speed/efficiency of the solver depends on how accurately it can carry out numerical integrals of the frequency and the first-derivative terms, therefore the **grid fineness** needs to be high enough. (Typically this means that linear interpolation gives a :math:`\omega(t)` value that is accurate to 1 part in :math:`10^{6}` or so.) If you want `oscode` to check whether the grids were sampled finely enough, set ``check_grid=true`` in the call for ``de_system()``.

To define the grids, use any array-like container which is **contiguous in
memory**, e.g. an ``Eigen::Vector``, ``std::array``, ``std::vector``:

.. code:: c
    
    #include "solver.hpp" // de_system, Solution defined in here 

    // Create a fine grid of timepoints and 
    // a grid of values for w, g
    N = 10000; 
    std::vector<double> ts(N);
    std::vector<std::complex<double>> ws(N), gs(N);
    
    // Fill up the grids
    for(int i=0; i<N; i++){
        ts[i] = i;
        ws[i] = std::sqrt(i);
        gs[i] = 0.0;
    }   

They can then be given to the solver again by feeding a pointer to their underlying
data to the ``de_system`` class:

.. code:: c
    
    de_system sys(ts.data(), ws.data(), gs.data());   
    Solution solution(sys, ...) // other arguments left out


Often :math:`\omega` and :math:`\gamma` are much easier to perform linear
interpolation on once taken natural log of. This is what the optional ``islogw``
and ``islogg`` arguments of the overloaded ``de_system::de_system()``
constructor are for:

.. code:: c
    
    #include "solver.hpp" // de_system, Solution defined in here 

    // Create a fine grid of timepoints and 
    // a grid of values for w, g
    N = 10000; 
    std::vector<double> ts(N);
    std::vector<std::complex<double> logws(N), gs(N); // Note the log!
    
    // Fill up the grids
    for(int i=0; i<N; i++){
        ts[i] = i;
        logws[i] = 0.5*i;
        gs[i] = 0.0; // Will not be logged
    }   
    
    // We want to tell de_system that w has been taken natural log of, but g
    // hasn't. Therefore islogw=true, islogg=false:
    de_system sys(ts.data(), logws.data(), gs.data(), true, false);
    Solution solution(sys, ... ) // other arguments left out


DIY interpolation
=================

For some problems, linear interpolation of :math:`\omega` and :math:`\gamma` (or
their natural logs) might simply not be enough.

For example, the user could carry out cubic spline interpolation and feed
:math:`\omega` and :math:`\gamma` as functions to ``de_system``. 

Another example for wanting to do (linear) interpolation outside of ``oscode`` is
when ``Solution.solve()`` is ran in a loop, and for each iteration a large grid
of :math:`\omega` and :math:`\gamma` is required, depending on some parameter.
Instead of generating them over and over again, one could define them as
functions, making use of some underlying vectors that are independent of the
parameter we iterate over:

.. code:: c

    // A, B, and C are large std::vectors, same for each run
    // k is a parameter, different for each run
    // the grid of timepoints w, g are defined on starts at tstart, and is
    // evenly spaced with a spacing tinc.

    // tstart, tinc, A, B, C defined here

    std::complex<double> g(double t){
        int i;
        i=int((t-tstart)/tinc);
        std::complex<double> g0 = 0.5*(k*k*A[i] + 3.0 - B[i] + C[i]*k;
        std::complex<double> g1 = 0.5*(k*k*A[i+1] + 3.0 - B[i+1] + C[i+1]*k);
        return (g0+(g1-g0)*(t-tstart-tinc*i)/tinc);
    };


Solving an equation
-------------------

Once the equation to be solver has been defined as an instance of the
``de_system`` class, the following additional information is necessary to solve
it: 

- initial conditions, :math:`x(t_i)` and :math:`\dot{x}(t_f)`,
- the range of integration, from :math:`t_i` and :math:`t_f`,
- (optional) set of timepoints at which dense output is required,
- (optional) order of WKB approximation to use, ``order=3``,
- (optional) relative tolerance, ``rtol=1e-4``,
- (optional) absolute tolerance ``atol=0.0``,
- (optional) initial step ``h_0=1``,
- (optional) output file name ``full_output=""``,

**Note** the following about the optional arguments:

- ``rtol``, ``atol`` are tolerances on the local error. The global error in the solution is not guaranteed to stay below these values, but the error per step is. In the RK regime (not oscillatory solution), the global error will rise above the tolerance limits, but in the WKB regime, the global error usually stagnates.
- The initial step should be thought of as an initial estimate of what the first stepsize should be. The solver will determine the largest possible step within the given tolerance limit, and change ``h_0`` if necessary.
- The full output of ``solve()`` will be written to the filename contained in ``full_output``, if specified.  

Here's an example to illustrate usage of all of the above variables:

.. code:: c
    
    #include "solver.hpp" // de_system, Solution defined in here 

    // Define the system
    de_system sys(...) // For args see previous examples

    // Necessary parameters:
    // initial conditions
    std::complex<double> x0=std::complex<double>(1.0,1.0), dx0=0.0;
    // range of integration
    double ti=1.0, tf=100.0;
    
    // Optional parameters:
    // dense output will be required at the following points:
    int n = 1000;
    std::vector t_eval(n);
    for(int i=0; i<n; i++){
        t_eval[i] = i/10.0;
    }
    // order of WKB approximation to use
    int order=2;
    // tolerances
    double rtol=2e-4, atol=0.0;
    // initial step
    double h0 = 0.5;
    // write the solution to a file
    std::string outfile="output.txt";

    Solution solution(sys, x0, dx0, ti, tf, t_eval.data(), order, rtol, atol, h0, outfile);
    // Solve the equation:
    solution.solve()

Here, we've also called the ``solve()`` method of the ``Solution`` class, to
carry out the integration. Now all information about the solution is in
``solution`` (and written to ``output.txt``).

Using the solution
------------------

Let's break down what ``solution`` contains (what ``Solution.solve()`` returns).
An instance of a ``Solution`` object is returned with the following attributes:

- ``times`` [std::list of double]: timepoints at which the solution was determined. These are **not** supplied by the user, rather they are internal steps that the solver has takes. The list starts with :math:`t_i` and ends with :math:`t_f`, these points are always guaranteed to be included.
- ``sol`` [std::list of std::complex<double>]: the solution at the timepoints specified in ``times``.
- ``dsol`` [std::list of std::complex<double>]: first derivative of the solution at timepoints specified in ``times``. 
- ``wkbs`` [std::list of int/bool]: types of steps takes at each timepoint in ``times``. **1** if the step was WKB, **0** if it was RK.  
- ``ssteps`` [int]: total number of accepted steps.  
- ``totsteps`` [int]: total number of attempted steps (accepted + rejected).  
- ``wkbsteps`` [int]: total number of successful WKB steps. 
- ``x_eval`` [std::list of std::complex<double>]: dense output, i.e. the solution evaluated at the points specified in the ``t_eval`` optional argument 
- ``dx_eval`` [std::list of std::complex<double>]: dense output of the derivative of the solution, evaluted at the points specified in ``t_eval`` optional argument.



pyoscode
========

.. automodule:: pyoscode
   :members:
   :undoc-members:
   :show-inheritance:
(py)oscode's documentation
==========================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

    Introduction <introduction>
    Python interface reference <pyoscode>
    Using the C++ interface <oscode>
    C++ interface reference <cpp-docs/oscode-reference>

.. title:: Introduction

========================================================================
(py)oscode: Oscillatory ordinary differential equation solver
========================================================================

.. image:: https://codecov.io/gh/fruzsinaagocs/oscode/branch/joss-paper/graph/badge.svg
    :target: https://codecov.io/gh/fruzsinaagocs/oscode
.. image:: https://travis-ci.org/fruzsinaagocs/oscode.svg?branch=master
    :target: https://travis-ci.org/fruzsinaagocs/oscode
    :alt: Travis CI build status
.. image:: https://readthedocs.org/projects/oscode/badge/?version=latest
    :target: https://oscode.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status
.. image:: https://badges.gitter.im/oscode-help/community.svg
    :target: https://gitter.im/oscode-help/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge
    :alt: Chat on gitter
.. image:: https://img.shields.io/badge/License-BSD%203--Clause-blue.svg
    :target: https://opensource.org/licenses/BSD-3-Clause
    :alt: BSD 3-clause license

|
|

.. contents::
   :local:

|

About
-----

``oscode`` is a C++ tool with a Python interface that solves **osc**\illatory
**o**\rdinary **d**\ifferential **e**\quations efficiently. It is designed to
deal with equations of the form

.. math:: 

	\ddot{x}(t) + 2\gamma(t)\dot{x}(t) + \omega^2(t)x(t) = 0,

where :math:`\gamma(t)` and :math:`\omega(t)` can be given as arrays.


``oscode`` makes use of an analytic approximation of :math:`x(t)` embedded in a
stepping procedure to skip over long regions of oscillations, giving a reduction
in computing time. The approximation is valid when the frequency
:math:`\omega(t)` changes slowly relative to the timescales of integration, it
is therefore worth applying when this condition holds for at least some part of
the integration range.

For the details of the numerical method used by ``oscode``, see the Citations
section.


Installation
------------

Dependencies
~~~~~~~~~~~~

Basic requirements for using the C++ interface:

- C++11 or later
- `Eigen <http://eigen.tuxfamily.org/index.php?title=Main_Page>`__ (a header-only library included in this source)

The strictly necessary Python dependencies are automatically installed when you use `pip` or the `setup.py`. They are:

- `numpy <https://pypi.org/project/numpy/>`__

The *optional* dependencies are: 

- for tests
    - `scipy <https://pypi.org/project/scipy/>`__ 
    - `pytest <https://docs.pytest.org/en/stable/getting-started.html>`__ 
- for examples/plotting
    - `matplotlib <https://pypi.org/project/matplotlib/>`__
    - `scipy <https://pypi.org/project/scipy/>`__ 
- for generating offline documentation
    - `sphinx <https://pypi.org/project/Sphinx/>`__ 
    - `doxygen <https://www.doxygen.nl/index.html>`__
    - `breathe <https://pypi.org/project/breathe/>`__
    - `exhale <https://pypi.org/project/exhale/>`__


Python
~~~~~~

``pyoscode`` can be installed via pip 

.. code:: bash

   pip install pyoscode

or via the setup.py

.. code:: bash

   git clone https://github.com/fruzsinaagocs/oscode
   cd oscode
   python setup.py install --user

You can then import ``pyoscode`` from anywhere. Omit the ``--user`` option if
you wish to install globally or in a virtual environment. If you have any
difficulties, check out the FAQs_ section below. 

You can check that things are working by running `tests/` (also ran by Travis continuous integration):

.. code:: bash

   pytest tests/


C++
~~~

``oscode`` is a header-only C++ package, it requires no installation.

.. code:: bash

   git clone https://github.com/fruzsinaagocs/oscode

and then include the relevant header files in your C++ code:

.. code:: c

    #include "solver.hpp"
    #include "system.hpp"


Quick start
-----------

Try the following quick examples. They are available in the `examples
<https://github.com/fruzsinaagocs/oscode/tree/master/examples/>`__.

Python
~~~~~~

:Introduction to pyoscode: |intro_binder|
:Cosmology examples: |cosmology_binder|

.. |intro_binder| image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/fruzsinaagocs/oscode/joss-paper?filepath=examples/introduction_to_pyoscode.ipynb

.. |cosmology_binder| image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/fruzsinaagocs/oscode/joss-paper?filepath=examples/cosmology.ipynb

.. image::
    https://github.com/fruzsinaagocs/oscode/raw/master/pyoscode/images/spectra.gif
    :width: 800


C++
~~~

:Introduction to oscode: `examples/burst.cpp`
:To plot results from `burst.cpp`: `examples/plot_burst.py`

To compile and run:

.. code:: bash

    g++ -g -Wall -std=c++11 -c -o burst.o burst.cpp
    g++ -g -Wall -std=c++11 -o burst burst.o
    ./burst


Documentation
-------------

Documentation is hosted at `readthedocs <https://oscode.readthedocs.io>`__.

To build your own local copy of the documentation you can run:

.. code:: bash

   cd pyoscode/docs
   make html

Citation
--------

If you use ``oscode`` to solve equations for a publication, please cite:

- `Efficient method for solving highly oscillatory ordinary differential equations with applications to physical systems <https://doi.org/10.1103/PhysRevResearch.2.013030>`__,
- `Dense output for highly oscillatory numerical solutions  <https://arxiv.org/abs/2007.05013>`__


Contributing
------------

Any comments and improvements to this project are welcome. You can contribute
by:

- Opening and `issue <https://www.github.com/fruzsinaagocs/oscode/issues/>`__ to report bugs and propose new features.
- Making a pull request.

Further help
------------

You can get help by submitting an issue or posting a message on `Gitter <https://gitter.im/oscode-help/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge>`__.


FAQs
----

Installation
~~~~~~~~~~~~

1. Eigen import errors:
    .. code:: bash

       pyoscode/_pyoscode.hpp:6:10: fatal error: Eigen/Dense: No such file or directory
        #include <Eigen/Dense>
                  ^~~~~~~~~~~~~

    Try explicitly including the location of your Eigen library via the
    ``CPLUS_INCLUDE_PATH`` environment variable, for example:

    .. code:: bash

       CPLUS_INCLUDE_PATH=/usr/include/eigen3 python setup.py install --user
       # or 
       CPLUS_INCLUDE_PATH=/usr/include/eigen3 pip install pyoscode

    where  ``/usr/include/eigen3`` should be replaced with your system-specific
    eigen location.


Thanks
------

Many thanks to **Will Handley**, **Lukas Hergt**, **Anthony Lasenby**, and **Mike Hobson** for
their support and advice regarding the algorithm behind `oscode`.
There are many packages without which some part of `oscode` (e.g. testing and
examples) wouldn't run as nicely and smoothly, thank you all developers for
making and maintaining these open-source projects. A special thanks goes to the
devs of `exhale <https://pypi.org/project/exhale/>`__ for making the beautiful C++ documentation possible. 


Changelog
---------

- 1.0.0: current version
    - Dense output
    - Arrays for frequency and damping term need not be evenly spaced
    - Automatic C++ documentation on readthedocs
    - Eigen included in source for pip installability
    - First pip release :)
- 0.1.2:
    - Bug that occurred when beginning and end of integration coincided
      corrected
- 0.1.1:
    - Automatic detection of direction of integration
- 0.1.0:
    - Memory leaks at python interface fixed
    - C++ documentation added 

