# Eigentools

Eigentools is a set of tools for studying linear eigenvalue problems. The underlying eigenproblems are solved using [Dedalus](http://dedalus-project.org), which provides a domain-specific language for partial differential equations. Eigentools extends Dedalus's `EigenvalueProblem` object and provides

* automatic rejection of unresolved eigenvalues
* simple plotting of specified eigenmodes
* simple plotting of spectra
* computation of pseudospectra for any Differential-Algebraic Equations with **user-specifiable norms**
* tools to find critical parameters for linear stability analysis
* ability to project eigenmode onto 2- or 3-D domain for visualization
* ability to output projected eigenmodes as Dedalus-formatted HDF5 file to be used as initial conditions for Initial Value Problems
* simple plotting of drift ratios (both ordinal and nearest) to evaluate tolerance for eigenvalue rejection

## Installation

Eigentools can be `pip` installed, though it requires [Dedalus](http://dedalus-project.org/), which has non-`pip` installable dependencies. See the [installation instructions](https://eigentools.readthedocs.io/en/latest/pages/installation.html) for details.

## Documentation

Documentation (including detailed API documentation) can be found at [Read the Docs](https://eigentools.readthedocs.io/).

If you are upgrading from version 1 to version 2, you can find a guide to API changes [here](https://eigentools.readthedocs.io/en/latest/pages/upgrading.html)

## Contributing

Eigentools welcomes community contributions from issue reports to code contributions. For details, please see [our contribution policy](CONTRIBUTING.md).

## Developers
The core development team consists of 

* Jeff Oishi (<jsoishi@gmail.com>)
* Keaton Burns (<keaton.burns@gmail.com>)
* Susan Clark (<susanclark19@gmail.com>)
* Evan Anders (<evan.anders@northwestern.edu>)
* Ben Brown (<bpbrown@gmail.com>)
* Geoff Vasil (<geoffrey.m.vasil@gmail.com>)
* Daniel Lecoanet (<daniel.lecoanet@northwestern.edu>)

## Support 
Eigentools was developed with support from the Research Corporation under award Scialog Collaborative Award (TDA) ID# 24231.


<!--  LocalWords:  Eigentools eigenproblems Dedalus EigenvalueProblem
 -->
<!--  LocalWords:  eigenmodes pseudospectra eigenmode HDF conda Oishi
 -->
<!--  LocalWords:  eigentools Anders Geoff Vasil Lecoanet Scialog TDA
 -->
# Contributing to Eigentools #

We welcome contributions, including issue reports, bug fixes, and feature implementations.
Contributions are reviewed on Github via pull request; to get started, fork the repository, make changes, and issue a pull request.
You can also contribute by submitting an issue.

## Reporting issues ##

If you find a bug or unexpected behavior, please file an issue report on [github](https://github.com/DedalusProject/eigentools/issues).
Please provide as much detail as possible, including version of both eigentools and dedalus, platform (Mac/Linux), and a stand alone `.py` file that demonstrates the problem in as simple a manner as possible.

## Proposing features ##

You can propose new features on the issue tracker using the "enhancement" tag.

## Contributing code ##

Code contributions ranging from fixing typos to implementing additional features are most welcome! We use [pull requests](https://github.com/DedalusProject/eigentools/pulls) to integrate contributions into the main codebase. If you would like to contribute and are looking for a place to start, please don't hesitate to contact the authors (emails are in [readme.md](readme.md)).

Eigentools
**********

Eigentools is a set of tools for studying linear eigenvalue problems. The underlying eigenproblems are solved using `Dedalus <http://dedalus-project.org>`_, which provides a domain-specific language for partial differential equations. Each entry in the following list of features links to a Jupyter notebook giving an example of its use.

* :ref:`automatic rejection of unresolved eigenvalues </notebooks/Orr Somerfeld pseudospectra.ipynb#Eigenmode-rejection>`
* :ref:`simple plotting of drift ratios (both ordinal and nearest) to evaluate tolerance for eigenvalue rejection </pages/getting_started.rst#Mode-rejection>`

* :ref:`simple plotting of specified eigenmodes </notebooks/Mixed Layer Instability.ipynb#Plotting-eigenmodes>`
* :ref:`simple plotting of spectra </notebooks/Orr Somerfeld pseudospectra.ipynb#Plotting-Spectra>`
* :ref:`computation of pseudospectra for any Differential-Algebraic Equations </notebooks/Orr Somerfeld pseudospectra.ipynb#Pseudospectra>` with :ref:`user-specifiable norms </notebooks/Orr Somerfeld pseudospectra.ipynb#Choosing-an-inner-product-and-norm>`
* :ref:`tools to find critical parameters for linear stability analysis </notebooks/Mixed Layer Instability.ipynb#Finding-critical-parameters>` with :ref:`user-specifiable definitions of growth and stability </notebooks/Mixed Layer Instability.ipynb#Specifying-a-definition-of-stability>`
* :ref:`ability to project eigenmode onto 2- or 3-D domain for visualization </notebooks/Mixed Layer Instability.ipynb#Projection-onto-higher-dimensional-domains>`
* :ref:`ability to output projected eigenmodes as Dedalus-formatted HDF5 file to be used as initial conditions for Initial Value Problems </notebooks/Mixed Layer Instability.ipynb#Writing-Dedalus-HDF5-files>`

Contents
========

.. toctree::
   :maxdepth: 2

   pages/installation
   pages/getting_started
   pages/upgrading

Example notebooks
-----------------

.. toctree::
   :maxdepth: 1

   Example 1: Orr-Somerfield pseudospectra </notebooks/Orr Somerfeld pseudospectra.ipynb>
   Example 2: Mixed Layer instability </notebooks/Mixed Layer Instability.ipynb>

API reference
-------------

.. toctree::
   :maxdepth: 2
              
   Eigentools API reference <autoapi/eigentools/index>

Developers
==========
The core development team consists of 

* Jeff Oishi (<jsoishi@gmail.com>)
* Keaton Burns (<keaton.burns@gmail.com>)
* Susan Clark (<susanclark19@gmail.com>)
* Evan Anders (<evan.anders@northwestern.edu>)
* Ben Brown (<bpbrown@gmail.com>)
* Geoff Vasil (<geoffrey.m.vasil@gmail.com>)
* Daniel Lecoanet (<daniel.lecoanet@northwestern.edu>)

Support
=======
Eigentools was developed with support from the Research Corporation under award Scialog Collaborative Award (TDA) ID# 24231.

Installing eigentools
*********************

eigentools requires Dedalus, which you can install via any of the methods found in `the Dedalus installation instructions <https://dedalus-project.readthedocs.io/en/latest/pages/installation.html>`_.

Once Dedalus is installed, eigentools is `pip` installable::

  pip install eigentools

If you would like the development version, you can clone the repository and install locally::

  git clone https://github.com/DedalusProject/eigentools.git
  pip install -e eigentools



Getting Started
***************

eigentools comes with several `examples <https://github.com/DedalusProject/eigentools/tree/master/examples>`_ to get you started, but let's outline some basics in a very simple problem, the 1-D wave equation with :math:`u = 0` at both ends.
This is not quite as trivial a problem as it might seem, because we are expanding the solution in Chebyshev polynomials, but the eigenmodes are sines and cosines. 


.. code-block:: python
                
    from eigentools import Eigenproblem
    import dedalus.public as de
    
    Nx = 128
    x = de.Chebyshev('x',Nx, interval=(-1, 1))
    d = de.Domain([x])
    
    string = de.EVP(d, ['u','u_x'], eigenvalue='omega')
    string.add_equation("omega*u + dx(u_x) = 0")
    string.add_equation("u_x - dx(u) = 0")
    string.add_bc("left(u) = 0")
    string.add_bc("right(u) = 0")
    
    EP = Eigenproblem(string)
    EP.solve(sparse=False)
    ax = EP.plot_spectrum()
    print("there are {} good eigenvalues.".format(len(EP.evalues)))
    ax.set_ylim(-1,1)
    ax.figure.savefig('waves_spectrum.png')

    ax = EP.plot_drift_ratios()
    ax.figure.savefig('waves_drift_ratios.png')

That code takes about 10 seconds to run on a 2020 Core-i7 laptop, produces about 68 "good" eigenvalues, and produces the following output:

.. image:: ../images/waves_spectrum.png
           :width: 400
           :alt: A spectrum for waves on a string

eigentools has taken a Dedalus eigenvalue problem, automatically run it at 1.5 times the specified resolution, rejected any eigenvalues that do not agree to a default precision of one part in :math:`10^{-6}` and plotted a spectrum in six extra lines of code!           

Most of the plotting functions in eigentools return a `matplotlib` `axes` object, making it easy to modify the plot defaults.
Here, we set the y-limits manually, because the eigenvalues of a string are real.
Try removing the `ax.set_ylim(-1,1)` line and see what happens.

Mode Rejection
--------------
One of the most important tasks eigentools performs is spurious mode rejection. It does so by computing the "drift ratio" [Boyd2000]_ between the eigenvalues at the given resolution and a higher resolution problem that eigentools automatically assembles. By default, the "high" resolution case is 1.5 times the given resolution, though this is user configurable via the `factor` keyword option to `Eigenproblem()`.

The drift ratio :math:`\delta` is calculated using either the **ordinal** (e.g. first mode of low resolution to first mode of high resolution) or **nearest** (mode with smallest difference between a given high mode and all low modes). In order to visualize this, `EP.plot_drift_ratios()` in the above code returns an `axes` object making a plot of the *inverse drift ratio* (:math:`1/\delta`),

.. image:: ../images/waves_drift_ratios.png
           :width: 400
           :alt: Plot of inverse drift ratios vs. mode number for waves on a string.

Good modes are those *above* the horizontal line at :math:`10^{6}`; bad modes are also grayed out. In this case, the **nearest** and **ordinal** methods produce identical results. If the problem contains more than one wave *family*, **nearest** typically fails. For an example, see the `MRI example script <https://github.com/DedalusProject/eigentools/blob/master/examples/mri.py>`_. Note that **nearest** is the default criterion used by eigentools.


.. [Boyd2000] Boyd, J (2000). "Chebyshev and Fourier Spectral Methods." Dover. `<http://www-personal.umich.edu/~jpboyd/aaabook_9500toc.pdf>`_
Upgrading eigentools scripts
****************************

Version 2 of eigentools has made significant changes to the API and will necessitate some changes (for the better, we hope) to the user experience. The guiding principle behind the new API is that one should no longer need to touch the Dedalus :code:`EVP` object that defines the eigenvalue problem at hand. 

**Most importantly**, no changes need to be made to the underlying Dedalus :code:`EVP` object.

Basic :code:`eigenproblem` usage
--------------------------------

Choosing a sparse or dense solve is no longer done when instantiating :code:`Eigenproblem` objects. Instead, this is a choice at *solve* time:

.. code-block:: python
                
    EP = Eigenproblem(string, reject=True)
    EP.solve(sparse=False)

Also, notice that rejection of spurious modes is now done automatically with :code:`EP.solve` if :code:`reject=True` is selected at instantiation time. Note that although in the above code, we explicitly set :code:`reject=True`, this is **unnecessary**, as it is the default. The :code:`EP.reject_spurious()` function has been removed

In addition, solving again with different parameters has been greatly simplified from the previous version. You now simply *pass a dictionary* with the parameters you wish to change to solve itself. Let's revisit the simple waves-on-a-string problem from :ref:`the getting started page </pages/getting_started.rst>`,  but add a parameter, :code:`c2`, the wave speed squared.

Here, we solve twice, once with :code:`c1 = 1` and once with :code:`c2 = 2`. Given the dispersion relation for this problem is :math:`\omega^2 = c^2 k` and our eigenvalue :code:`omega` is really :math:`\omega^2`, we expect the eigenvalues for the second solve to be twice those for the first.

.. code-block:: python
                
    import numpy as np
    from eigentools import Eigenproblem
    import dedalus.public as de
    
    Nx = 128
    x = de.Chebyshev('x',Nx, interval=(-1, 1))
    d = de.Domain([x])
    
    string = de.EVP(d, ['u','u_x'], eigenvalue='omega')
    string.parameters['c2'] = 1
    string.add_equation("omega*u + c2*dx(u_x) = 0")
    string.add_equation("u_x - dx(u) = 0")
    string.add_bc("left(u) = 0")
    string.add_bc("right(u) = 0")
    EP = Eigenproblem(string)
    EP.solve(sparse=False)
    evals_c1 = EP.evalues
    EP.solve(sparse=False, parameters={'c2':2})
    evals_c2 = EP.evalues
    
    print(np.allclose(evals_c2, 2*evals_c1))

Getting eigenmodes
==================

Getting eigenmodes has also been simplified and significantly extended. Previously, getting an eigenmode corresponding to an eigenvalue required using the :code:`set_state()` method on the underlying :code:`EVP` object. In keeping with the principle of not needing to manipulate the :code:`EVP`, we provide a new :code:`.eigenmode(index)`, where :code:`index` is the mode number corresponding to the eigenvalue index in :code:`EP.evalues`. By default, with mode rejection on, these are the "good" eigenmodes.
`.eigenmode(index)` returns a Dedalus :code:`FieldSystem` object, with a Dedalus :code:`Field` for each field in the eigenmode:

.. code-block:: python
                
    emode = EP.eigenmode(0)
    print([f.name for f in emode.fields])
    u = emode.fields[0]
    u_x = emode.fields[1]


Finding critical parameters
---------------------------

This has been considerably cleaned up. The two major things to note are that

1. one no longer needs to create a shim function to translate between an x-y grid and the parameter names within the :code:`EVP`.
2. The parameter grid is no longer defined inside :code:`CriticalFinder`, but is instead created by the user and passed in

For example, here are the relevant changes necessary for the `MRI test problem <https://github.com/DedalusProject/eigentools/tree/master/examples/mri.py>`_.

First, replace

.. code-block:: python
                
    EP = Eigenproblem(mri, sparse=True)
    
    # create a shim function to translate (x, y) to the parameters for the eigenvalue problem:
    def shim(x,y):
        iRm = 1/x
        iRe = (iRm*Pm)
        print("Rm = {}; Re = {}; Pm = {}".format(1/iRm, 1/iRe, Pm))
        gr, indx, freq = EP.growth_rate({"Q":y,"iRm":iRm,"iR":iRe})
        ret = gr+1j*freq
        return ret
     
    cf = CriticalFinder(shim, comm)

with

.. code-block:: python
   
    EP = Eigenproblem(mri)

    cf = CriticalFinder(EP, ("Q", "Rm"), comm, find_freq=False)

**Important:** note that :code:`find_freq` is specified at instantiation rather than when calling :code:`cf.crit_finder` later.

Once this is done, the grid generation changes from

.. code-block:: python
                
    mins = np.array((4.6, 0.5))
    maxs = np.array((5.5, 1.5))
    ns   = np.array((10,10))
    logs = np.array((False, False))
    
    cf.grid_generator(mins, maxs, ns, logs=logs)

to

.. code-block:: python
                
    nx = 20
    ny = 20
    xpoints = np.linspace(0.5, 1.5, nx)
    ypoints = np.linspace(4.6, 5.5, ny)
    
    cf.grid_generator((xpoints, ypoints), sparse=True)


