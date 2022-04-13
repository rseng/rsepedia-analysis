# Contributing to `NGSTrefftz`

🎉**Thanks for taking the time to contribute!** 🎉

To get an overview of the project, check out the [README](README.md).

The [issue tracker](https://github.com/paulst/ngstrefftz/issues)
is the preferred channel for [bug reports](#bugs), [features
requests](#features) and [submitting pull requests](#pull-requests). 

For personal support requests you can reach me at [p.stocker@math.uni-goettingen.de](mailto:p.stocker@math.uni-goettingen.de).

## Bug reports

A bug is a _demonstrable problem_ that is caused by the code in the repository.
Bug reports are extremely helpful - thank you!

Guidelines for bug reports:

1. **Check if the issue has been fixed**: try to reproduce it using the latest `main` or development branch in the repository.

2. **Use the GitHub issue search**: check if the issue has already been reported.

3. **Isolate the problem**: Create a minimal example showing the problem.

4. **Open an issue**: Using the [issue tracker](https://github.com/paulst/ngstrefftz/issues), describe the expected outcome and report the OS, the compiler, and NGSolve version you are using.

## Pull requests

Pull requests - patches, improvements, new features - are a fantastic
help. They should remain focused in scope and avoid containing unrelated commits.
**Please ask first** before embarking on any significant pull request.

Tips on opening a pull request:

1. [Fork](http://help.github.com/fork-a-repo/) the project.

2. Create a branch and implement your feature.
   ```bash
   git checkout -b <your-feature-name>
   ```

3. Run the test suite by calling 
   ```bash
   make test
   ```
   in your build directory. Consider adding new tests for you feature - have a look in the test folder.

4. Once the implementation is done, use Git's
   [interactive rebase](https://help.github.com/articles/interactive-rebase)
   feature to tidy up your commits.
   ```bash
   git rebase --interactive --fork-point main <your-feature-name> 
   ```

5. Push your topic branch up to your fork and [open a Pull Request](https://help.github.com/articles/using-pull-requests/).

**IMPORTANT**: By submitting a patch, you agree to allow the project owners to license your work under the terms of the *LGPL License*.

# NGSTrefftz
**an add-on to NGSolve for Trefftz methods**

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/PaulSt/NGSTrefftz/HEAD?filepath=doc%2Fnotebooks%2Findex.ipynb)
[![Docker Image Version (latest semver)](https://img.shields.io/docker/v/paulstdocker/ngstrefftz?label=docker&logo=docker&sort=semver)](https://hub.docker.com/r/paulstdocker/ngstrefftz)
[![PyPI](https://img.shields.io/pypi/v/ngstrefftz?color=blue&logo=pypi)](https://pypi.org/project/ngstrefftz/)
[![GitHub Workflow Status](https://img.shields.io/github/workflow/status/PaulSt/NGSTrefftz/build?logo=github)](https://github.com/PaulSt/NGSTrefftz/actions)
[![status](https://joss.theoj.org/papers/c2f4e85b118c22b81aa27d7799265409/status.svg)](https://joss.theoj.org/papers/c2f4e85b118c22b81aa27d7799265409)
[![docs](https://img.shields.io/badge/docs-NGSTrefftz-blue?logo=readthedocs)](https://paulst.github.io/NGSTrefftz/)

NGSTrefftz provides a framework to implement Trefftz finite element spaces for [NGSolve](https://www.ngsolve.org), with several Trefftz spaces already implemented. Additionally, Trefftz-DG on tent-pitched meshes for the acoustic wave equation is implemented using meshes provided by [ngstents](https://github.com/jayggg/ngstents). Furthermore, the package includes an implementation of the embedded Trefftz method.

## Try it out!
You can try out some jupyter notebooks:
* Launch the Binder here:   
  [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/PaulSt/NGSTrefftz/HEAD?filepath=doc%2Fnotebooks%2Findex.ipynb)
* Or run the docker locally (you need to have docker installed):

```bash
git clone https://github.com/PaulSt/NGSTrefftz
cd NGSTrefftz && docker build -t ngstrefftz_jupyter .
docker run -p 8888:8888 ngstrefftz_jupyter
```

## Installing the package
You can either:
 * install using pip

```bash
pip install ngstrefftz
```

 * or build from source

```bash
git clone --recursive https://github.com/PaulSt/NGSTrefftz
mkdir ./NGSTrefftz/make && cd ./NGSTrefftz/make
cmake ../src && make install
```

### Dependencies
To build the package the following software needs to be installed
 * [cmake](https://cmake.org/) >= 3.1
 * [gcc](https://gcc.gnu.org/) >= 9 or [clang](https://clang.llvm.org/) >= 10
 * [ngsolve](https://www.ngsolve.org) >= 6.2

To access the newest features the nightly version of NGSolve works best and [lapack](http://www.netlib.org/lapack/) >= 3.9 is required.

## News
Mar, 2022: NGSTrefftz now has a [website](https://paulst.github.io/NGSTrefftz/)! 🚀

Feb, 2022: If you are using NGSolve nightly releases: [NGSolve@eda758d](https://github.com/NGSolve/ngsolve/commit/eda758d99483888851913d8a5c9aff4d0cbc9cc2) breaks a dependency and [NGSolve@3d52ecd](https://github.com/NGSolve/ngsolve/commit/3d52ecd615f2b7c409219eebaba99288ea19c1bc) produces import issue. Make sure to update ngstrefftz submodules and move to newest ngsolve version, at least [NGSolve@5839a09](https://github.com/NGSolve/ngsolve/commit/5839a09810235a938bd85807d8e638d3a0b6c69d).

Jan, 2022: NGSTrefftz is now available via pip! 🚀

Nov, 2021: NGSTrefftz now comes in a docker and with binder notebooks! 🚀

## Papers using the code
* Tent pitching and Trefftz-DG method for the acoustic wave equation  
[![arXiv](https://img.shields.io/badge/arXiv-1907.02367-b31b1b.svg)](https://arxiv.org/abs/1907.02367)
* A space-time quasi-Trefftz DG method for the wave equation with piecewise-smooth coefficients  
[![arXiv](https://img.shields.io/badge/arXiv-2011.04617-b31b1b.svg)](https://arxiv.org/abs/2011.04617)
* Embedded Trefftz discontinuous Galerkin methods  
[![arXiv](https://img.shields.io/badge/arXiv-2201.07041-b31b1b.svg)](https://arxiv.org/abs/2201.07041)


![](.github/wave.gif)

---
title: '`NGSTrefftz`: Add-on to NGSolve for Trefftz methods'
tags:
  - numerical methods for PDEs
  - finite elements 
  - Trefftz methods
authors:
  - name: Paul Stocker^[corresponding author] 
    orcid: 0000-0001-5073-3366 
    affiliation: "1" # (Multiple affiliations must be quoted)
affiliations:
 - name: Georg-August-Universit&#0228;t, G&#0246;ttingen, Germany
   index: 1
date: 29 January 2021
bibliography: paper.bib

---

# Summary
`NGSTrefftz` is an add-on to [`Netgen/NGSolve`](https://ngsolve.org/), a finite element software for the numerical treatment of partial differential equations (PDEs).
The package implements Trefftz based discontinuous Galerkin (DG) methods in `NGSolve`.
Trefftz methods reduce the number of unknowns in the discretization of PDE problems by injecting knowledge of the PDE into the approximation functions.
Like `NGSolve`, `NGSTrefftz` is written in C++ and integrates seamlessly with the easy-to-use Python interface of `NGSolve`.

Trefftz methods originate from @trefftz1926 and have since been developed for a wide range of problems, for an overview see @cbe1997trefftz, @TrefftzSurvey, @Qin05, @LLHC08, @kk95.
The central principle of Trefftz methods is the construction of a discrete basis of solutions to the differential operator under consideration, making the space of Trefftz functions problem dependent. 
In combination with finite elements the Trefftz basis is constructed locally, on each mesh element, continuity and boundary conditions are then enforced in the variational formulation.

# Statement of need

`NGSTrefftz` provides a framework to implement Trefftz finite element spaces for `NGSolve`, with spaces for Laplace equation, Helmholtz equation, and acoustic wave equation already implemented. 
`NGSolve` provides a flexible framework to implement variational formulations for arbitrary physical models. 
The Trefftz finite element spaces can then be used with the (bi-)linear forms generated in `NGSolve`.
The focus lies on the combination of Trefftz functions with DG methods or least-squares formulations,
this approach is also often referred to as frameless Trefftz elements.

On top of that, the package provides unique features that are:

* A quasi-Trefftz space that provides Trefftz-like properties for the acoustic wave equation with smooth coefficients. The method and results are presented in @qtrefftz.

* A space--time Trefftz method for the acoustic wave equation on tent-pitched meshes. 
Tent-pitching is a space--time meshing strategy that provides mesh elements that conform to the causality constraint of a hyperbolic system. 
The meshes are generated using [`ngstents`](https://github.com/jayggg/ngstents) and can be used with the Trefftz and quasi-Trefftz space for the acoustic wave equation, results are shown in @StockerSchoeberl and @qtrefftz.

* A general framework to produce Trefftz spaces implicitly is provided by an implementation of the embedded Trefftz method, see @embtrefftz.
The approach produces a Galerkin projection of an underlying discontinuous Galerkin method onto a subspace of Trefftz-type. 
It can be applied to very general cases, including inhomogeneous sources and non-constant coefficient differential operators.

The aim of this package is to facilitate research into Trefftz methods and to make them more accessible to a broad audience.
To the best of our knowledge, the only other open source software package that provides Trefftz finite element methods is [`FreeHyTe`](https://www.sites.google.com/site/ionutdmoldovan/freehyte) [@freehyte], implemented in `MATLAB`. 
Examples for the usage of all the features in `NGSTrefftz` are provided in the form of jupyter notebooks.

# Acknowledgements
P. Stocker has been supported by the Austrian Science Fund (FWF) through the projects F 65, W 1245, and by the German Research Foundation (DFG) through grant 432680300 - SFB 1456.
Special thanks to the `NGSolve` team for their continuous support.

# References
.. NGSTrefftz documentation master file, created by
   sphinx-quickstart on Wed Mar  2 15:24:26 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. raw:: html

    <style>
        .p-Widget {
            height: 400px;
        }
        .dg.main {
            margin-left: 0px;
        }
        div.p-Widget div div div div.dg ul li {
            list-style: none;
            margin-left: 0px;
        }
        div.p-Widget div div div div.dg ul li div.dg {
            margin-bottom: 0px;
        }
    </style>

.. only:: html
    .. role:: raw-html(raw)
        :format: html


.. mdinclude:: ../README.md


Introduction to Trefftz-DG
========================================

We give a short introduction to Trefftz-DG methods for the model problem of the Laplace equation.

.. math::

    \newcommand{\Th}{{\mathcal{T}_h}} 
    \newcommand{\Fh}{\mathcal{F}_h} 
    \newcommand{\dom}{\Omega} 
    \newcommand{\jump}[1]{[\![ #1 ]\!]}
    \newcommand{\tjump}[1]{[\![{#1} ]\!]_\tau}
    \newcommand{\avg}[1]{\{\!\!\{#1\}\!\!\}}
    \newcommand{\nx}{n_\mathbf{x}} 
    \begin{align*} \begin{split}
        \begin{cases}
        -\Delta u = 0 &\text{ in } \Omega, \\
        u=g &\text{ on } \partial \Omega,
        \end{cases}
    \end{split} \end{align*}


Harmonic polynomials
----------------------------------------
For an element :math:`K` in out mesh :math:`\Th` we define the local Trefftz space as

.. math::

    \begin{align*} \begin{split}
    \mathbb{T}^p(K):=\big\{
    f\in\mathbb{P}^p(K) \mid \Delta f = 0
    \big\},
    \qquad p\in \mathbb{N}.
    \end{split} \end{align*}


For the Laplace problem the Trefftz space is given by harmonic polynomials, e.g. for polynomials of order 3 in two dimensions the space is given by

.. math::

    \begin{align*} \begin{split}
    \mathbb{T}^3(K)=\text{span}\{1,x,y,xy,x^2-y^2,x^3-3y^2x,y^3-3x^2y\}
    \end{split} \end{align*}

The local degrees of freedom are reduced from the dimension of the full polynomials, given by :math:`\big(\begin{smallmatrix}p+n\\p\end{smallmatrix}\big)`, to

.. math::

    \begin{align*} \begin{split}
    \dim\mathbb{T}^p(K)=\begin{cases} 2 & n=1\\ 2p+1 & n=2\\ (p+1)^2 & n=3\end{cases}.
    \end{split} \end{align*}

Numerical treatment
----------------------------------------

The IP-DG method for the Laplace equation is given by

.. math::

    \begin{align}\label{eq:dglap}
        \begin{split}
        a_h(u,v) &= \int_\dom \nabla u\nabla v\ dV
        -\int_{\Fh^\text{int}}\left(\avg{\nabla u}\jump{v}+\avg{\nabla v}\jump{u} 
        - \frac{\alpha p^2}{h}\jump{u}\jump{v} \right) dS \\
               &\qquad -\int_{\Fh^\text{bnd}}\left(\nx\cdot\nabla u v+\nx\cdot\nabla v u-\frac{\alpha p^2}{h} u v \right) dS\\
        \ell(v) &= \int_{\Fh^\text{bnd}}\left(\frac{\alpha p^2}{h} gv -\nx\cdot\nabla vg\right) dS.
        \end{split}
    \end{align}

where :math:`h` is the mesh size, :math:`p` is the polynomial degree and :math:`\alpha>0`.

.. jupyter-execute::
    :hide-output:
    :hide-code:

    from ngsolve import *
    from ngstrefftz import *
    from netgen.geom2d import unit_square
    p=5

For a given mesh we construct the Trefftz finite element space as 

.. math::

    \begin{align*} \begin{split}
    \mathbb{T}^p(\Th):=\prod_{K\in\Th} \mathbb{T}^p(K)
    \end{split} \end{align*}

in NGSolve this is done by

.. jupyter-execute::
    :hide-output:

    mesh = Mesh(unit_square.GenerateMesh(maxh=.5))
    trefftzfespace(mesh,order=p,eq="laplace",dgjumps=True)

We will compare the approximation properties of the Trefftz space with the space of (all) piecewise polynomials. In NGSolve

.. jupyter-execute::
    :hide-output:

    L2(mesh,order=p,dgjumps=True)

.. jupyter-execute::
    :hide-code:
    :hide-output:

    exact = exp(x)*sin(y)
    bndc = exact

    def lap_problem(fes):
        mesh = fes.mesh
        order = fes.globalorder
        alpha = 4
        n = specialcf.normal(mesh.dim)
        h = specialcf.mesh_size
        u = fes.TrialFunction()
        v = fes.TestFunction()

        jump_u = (u-u.Other())*n
        jump_v = (v-v.Other())*n
        mean_dudn = 0.5 * (grad(u)+grad(u.Other()))
        mean_dvdn = 0.5 * (grad(v)+grad(v.Other()))

        a = BilinearForm(fes,symmetric=True)
        a += grad(u)*grad(v) * dx \
            +alpha*order**2/h*jump_u*jump_v * dx(skeleton=True) \
            +(-mean_dudn*jump_v-mean_dvdn*jump_u) * dx(skeleton=True) \
            +alpha*order**2/h*u*v * ds(skeleton=True) \
            +(-n*grad(u)*v-n*grad(v)*u)* ds(skeleton=True)

        f = LinearForm(fes)
        f += alpha*order**2/h*bndc*v * ds(skeleton=True) \
             +(-n*grad(v)*bndc)* ds(skeleton=True)

        with TaskManager():
            a.Assemble()
            f.Assemble()
        return (a,f)

    import time
    fes = L2
    tfes = trefftzfespace

    terrors = []
    tndofs = []
    ttime = []
    l2errors = []
    l2ndofs = []
    l2time = []
    def eoc(FES,errors,ndofs,timers):
        for p in range(1,6):
            errors.append([])
            ndofs.append([])
            timers.append([])
            mesh = Mesh(unit_square.GenerateMesh(maxh=1))
            for h in range(6):
                with TaskManager():
                    fes = FES(mesh,order=p,eq="laplace",dgjumps=True)
                    (a,f) = lap_problem(fes)
                    gfu = GridFunction(fes)
                    start = time.time()
                    gfu.vec.data = a.mat.Inverse(inverse='sparsecholesky') * f.vec
                    timers[p-1].append(time.time()-start)
                ndofs[p-1].append(fes.ndof)
                errors[p-1].append(sqrt(Integrate((gfu-exact)**2, mesh)))
                mesh.Refine()
    eoc(fes,l2errors,l2ndofs,l2time)
    eoc(tfes,terrors,tndofs,ttime)

The numerical results for the Trefftz space are plotted in solid lines, while the results for the full polynomial space are the dashed lines.
We show the convergence rates with respect to polynomial degree :math:`p` and mesh size :math:`h`. 
In the Figure on the left we show the error compared to the global number of degrees of freedom for varying :math:`p`.
In the Figure on the right we show the error with respect to :math:`h`.

.. jupyter-execute::
    :hide-code:

    import matplotlib.pyplot as plt
    #colors = ['b','g','r','c','m','y']
    colors = ['C'+str(i) for i in range(6)]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
    fig.tight_layout() 
    fig.subplots_adjust(wspace=0.2)

    for c,h in zip(colors,range(3,6)):
        ax1.loglog([tndofs[i][h] for i in range(5)], [terrors[i][h] for i in range(5)],c+'-o',label='$h = 2^{-'+str(h)+'}$')
        ax1.loglog([l2ndofs[i][h] for i in range(5)], [l2errors[i][h] for i in range(5)],c+'--o')
    ax1.set(xlabel="#dofs", ylabel="$L^2$-error", title="$p$-convergence")
    ax1.grid(True)
    ax1.legend()

    for c,p in zip(colors,range(1,6)):
        ax2.loglog([0.5**i for i in range(6)], terrors[p-1],c+'-o',label=f'$p = {p}$')
        ax2.loglog([0.5**i for i in range(6)], l2errors[p-1],c+'--o')
    ax2.set(xlabel="h", ylabel="$L^2$-error", title="$h$-convergence")
    ax2.grid(True)
    ax2.legend()
    ax2.invert_xaxis()

The Trefftz space shows optimal rate of convergence, using fewer degrees of freedom.
As an exact solution we used

.. math::

    u(x,y) = \exp(x)\sin(y) \quad \text{in} \quad [0,1]^2.

.. jupyter-execute::

    from ngsolve.webgui import Draw
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.3))
    fes = trefftzfespace(mesh,order=5,eq="laplace",dgjumps=True)
    (a,f) = lap_problem(fes)
    gfu = GridFunction(fes)
    gfu.vec.data = a.mat.Inverse(inverse='sparsecholesky') * f.vec
    Draw(gfu)


How to get started?
----------------------------------------

- If you want to see more implementations using Trefftz-DG methods or use an already implemented Trefftz space, have a look at the `notebooks`_ and the `documentation`_.

- If you are looking to implement a (polynomial) Trefftz space, a good starting point is to have a look at `trefftzspace.hpp <https://github.com/PaulSt/NGSTrefftz/blob/main/src/trefftzfespace.hpp>`_. For Trefftz spaces based on plane waves check out `planewavefe.hpp <https://github.com/PaulSt/NGSTrefftz/blob/main/src/planewavefe.hpp>`_.

- Before undertaking the implementation of a new Trefftz space in C++, or especially if your PDE does not have a easy-to-construct Trefftz space, consider checking out the `embedded Trefftz method`_.

If you are implementing a new method using NGSTrefftz consider `contributing`_.

.. _notebooks: notebooks/index.html
.. _embedded Trefftz method: notebooks/embTrefftz.html


.. _documentation: 
Documentation
========================================

Trefftz spaces
----------------------------------------
Trefftz finite element spaces extend the ngsolve class of :code:`FESpace`. They provide test and trial functions that are elementwise solutions to chosen PDEs. They are ment to be used in a DG setting with :code:`dgjumps=True`.

.. autoclass:: ngstrefftz.trefftzfespace(mesh,**kwargs)


Embedded Trefftz method
----------------------------------------
The embedded Trefftz method produces a Galerkin projection of an underlying discontinuous Galerkin method onto a subspace of Trefftz-type. 
It can be applied to very general cases, including inhomogeneous sources and non-constant coefficient differential operators.
:math:`\newcommand{\Th}{\mathcal{T}_h}
\newcommand{\Vhp}{V^p(\Th)}
\newcommand{\bT}{\mathbf{T}}
\newcommand{\bW}{\mathbf{W}}
\newcommand{\bl}{\mathbf{l}}
\newcommand{\bM}{\mathbf{M}}
\newcommand{\bL}{\mathbf{L}}
\newcommand{\bA}{\mathbf{A}}
\newcommand{\bU}{\mathbf{U}}
\newcommand{\bV}{\mathbf{V}}
\newcommand{\calL}{\mathcal{L}}
\newcommand{\bu}{\mathbf{u}}
\newcommand{\IT}{\mathbb{T}}
\newcommand{\calG}{\mathcal{G}}
\newcommand{\be}{\mathbf{e}}
\newcommand{\bx}{{\mathbf x}}
\newcommand{\inner}[1]{\langle #1 \rangle}
\DeclareMathOperator\Ker{ker}`

Specifically, for an operator :math:`\calL` we are looking to find a basis for its kernel over the space of polynomails. We call this projection :math:`\bT`
Then we can solve the reduced problem: Find :math:`\bu_\IT` so that 

.. math::

    \begin{equation} \label{eq:trefftzlinearsys}
       \bT^T\bA\bT ~ \bu_\IT = \bT^T \bl.
    \end{equation}

The solution in the full polynomial space is then given by :math:`\bu=\bT\bu_\IT`


.. autofunction:: ngstrefftz.TrefftzEmbedding

Trefftz + tent-pitching 
----------------------------------------
Given a tent-pitched mesh produced by ngstents the following function returns a solver using Trefftz, or quas-Trefftz, finite elements to solve the acoustic wave equation.

.. autofunction:: ngstrefftz.TWave

The solver class returned by :code:`TWave` is used to set initial and boundary conditions and to propagate the solution in the given tent slab.

.. autoclass:: ngstrefftz.TrefftzTents

   .. automethod:: Propagate
   .. automethod:: SetInitial
   .. automethod:: SetBoundaryCF


.. _contributing: 
.. mdinclude:: ../CONTRIBUTING.md


NGSTrefftz Notebooks
========================================

Take a look at the `notebooks`_! 
You can run them online `here`_. 

.. _here: https://mybinder.org/v2/gh/PaulSt/NGSTrefftz/HEAD?filepath=doc%2Fnotebooks%2Findex.ipynb

.. toctree::
   :maxdepth: 1

   notebooks/laplace.ipynb
   notebooks/helmholtz.ipynb
   notebooks/twave.ipynb
   notebooks/qtwave.ipynb
   notebooks/twavetents.ipynb
   notebooks/embTrefftz.ipynb
   notebooks/embTrefftz-poi.ipynb
   notebooks/embTrefftz-wave.ipynb
   notebooks/embTrefftz-helm.ipynb
   notebooks/embTrefftz-adv.ipynb
