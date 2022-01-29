---
title: '``pytreegrav``: A fast Python gravity solver'
tags:
  - Python
  - physics
  - gravity
  - simulations
authors:
  - name: Michael Y. Grudić
    orcid: 0000-0002-1655-5604
    affiliation: "1,2"
  - name: Alexander B. Gurvich
    orcid: 0000-0002-6145-3674
    affiliation: 2
affiliations:
 - name: NASA Hubble Fellow, Carnegie Observatories
   index: 1
 - name: Department of Physics & Astronomy and CIERA, Northwestern University
   index: 2
date: 9 June 2021
bibliography: paper.bib
---

# Summary

Gravity is important in a wide variety of science problems. In particular, questions in astrophysics nearly all involve gravity, and can have large ($\gg10^4$) numbers of gravitating masses, such as the stars in a cluster or galaxy, or the discrete fluid elements in a hydrodynamics simulation. Often the gravitational field of such a large number of masses can be too computationally expensive to compute by directly summing the contribution of every single element at every point of interest.

``pytreegrav`` is a multi-method Python package for computing gravitational fields and potentials. It includes an exact direct-summation ("brute force") solver and a fast, approximate tree-based method that can be orders of magnitude faster than the naïve method. It can compute fields and potentials from arbitrary particle distributions at arbitrary points, with arbitrary softening/smoothing lengths, and is parallelized with OpenMP.

# Statement of need

The problem addressed by ``pytreegrav`` is the following: given an arbitrary set of "source" masses $m_i$ with 3D coordinates $\mathbf{x}_i$, and optionally each having a finite spatial extent $h_i$ (the _softening radius_), one would like to compute the gravitational potential $\Phi$ and/or the gravitational field $\mathbf{g}$ at an arbitrary set of "target" points in space $\mathbf{y}_i$. A common application for this is N-body simulations (wherein $\mathbf{y}_i=\mathbf{x}_i$). It is also often useful for _analyzing_ simulation results after the fact -- $\Phi$ and $\mathbf{g}$ are sometimes not saved in simulation outputs, and even when they are it is often useful to analyze the gravitational interactions between specific _subsets_ of the mass elements in the simulation. Computing $\mathbf{g}$ is also important for generating equilibrium _initial conditions_ for N-body simulations [@makedisk;@galic], and for identifying interesting gravitationally-bound structures such as halos, star clusters, and giant molecular clouds [@rockstar;@grudic2018;@guszejnov2020].

Many gravity simulation codes (or multi-physics simulation codes _including_ gravity) have been written that address the problem of gravity computation in a variety of ways for their own internal purposes [@aarseth_nbody;@dehnen]. However, ``pykdgrav`` (the precursor of ``pytreegrav``) was the first Python package to offer a generic, modular, trivially-installable gravity solver that could be easily integrated into any other Python code, using the fast, approximate tree-based @barneshut method to be practical for large particle numbers. ``pykdgrav`` used a KD-tree implementation accelerated with ``numba`` [@numba] to achieve high performance in the potential/field evaluation, however the prerequisite tree-building step had relatively high overhead and a very large memory footprint, because the entire dataset was redundantly stored at every level in the tree hierarchy. This made it difficult to scale to various practical research problems, such as analyzing high-resolution galaxy simulations [@fire_pressurebalance]. ``pytreegrav`` is a full refactor of ``pykdgrav`` that addresses these shortcomings with a new octree implementation, with drastically reduced tree-build time and memory footprint, and a more efficient non-recursive tree traversal for field summation. This makes it suitable for post-processing datasets from state-of-the-art astrophysics simulations, with upwards of $10^8$ particles in the region of interest. 

# Methods

``pytreegrav`` can compute $\Phi$ and $\mathbf{g}$ using one of two methods: by "brute force" (explcitly summing the field of every particle, which is exact to machine precision), or using the fast, approximate @barneshut tree-based method (which is approximate, but much faster for large particle numbers). In $N$-body problems where the fields at all particle positions must be known, the cost of the brute-force method scales as $\propto N^2$, while the cost of the tree-based method scales less steeply, as $\propto N \log N$.

![Wall-clock time per particle running ``pytreegrav`` on a sample of $N$ particles from a @plummer distribution for various $N$. Test was run on an Intel i9 9900K workstation on a single core (_left_) and in parallel on 16 logical cores (_right_).\label{fig:cputime}](images/CPU_Time_both.png)

The brute-force methods are often fastest for small ($<10^3$ particle) point sets because they lack the overheads of tree construction and traversal, while the tree-based methods will typically be faster for larger datasets because they reduce the number of floating-point operations required. Both methods are optimized with the ``numba`` LLVM JIT compiler [@numba], and the basic ``Accel`` and ``Potential`` front-end functions will automatically choose the method is likely to be faster, based on this heuristic crossover point of $10^3$ particles. Both methods can also optionally be parallelized with OpenMP, via the ``numba`` ``@njit(parallel=True)`` interface.

The implementation of the tree build and tree-based field summation largely follows that of ``GADGET-2`` [@gadget2]. Starting with an initial cube enclosing all particles, particles are inserted into the tree one at a time. Nodes are divided into 8 subnodes until each subnode contains at most one particle. The indices of the 8 subnodes of each node are stored for an initial recursive traversal of the completed tree, but an optimized tree traversal only needs to know the _first_ subnode (if the node is to be refined) and the index of the next branch of the tree (if the field due to the node is summed directly), so these indices are recorded in the initial recursive tree traversal, and the 8 explicit subnode indices are then deleted, saving memory and removing any empty nodes from consideration. Once these "next branch" and "first subnode" indices are known, the tree field summations can be done in a single ``while`` loop with no recursive function calls, which generally improves performance and memory usage.

The field summation itself uses the @barneshut geometric opening criterion, with improvements suggested by @dubinski: for a node of side length $L$ with centre of mass located at distance $r$ from the target point, its contribution is summed using the monopole approximation (treating the whole node as a point mass) only if $r > L/\Theta + \delta$, where $\Theta=0.7$ by default (giving $\sim 1\%$ RMS error in $\mathbf{g}$), $\delta$ is the distance from the node's geometric center to its center of mass. If the conditions for approximation are not satisfied, the node's subnodes are considered in turn, until the field contribution of all mass within the node is summed.

``pytreegrav`` supports gravitational softening by assuming the mass distribution of each particle takes the form of a standard M4 cubic spline kernel, which is zero beyond the softening radius $h$ (outside which the field reduces to that of a point mass). Explicit expressions for this form of the softened gravitational potential and field are given in @gizmo. $h$ is allowed to vary from particle to particle, and when summing the field the larger of the source or the target softening is used (symmetrizing the force between overlapping particles). When softenings are nonzero, the largest softening $h_{\rm max}$ of all particles in a node is stored, and a node is always opened in the field summation if $r < 0.6L + \max\left(h_{\rm target}, h_{\rm max}\right) + \delta$, where $h_{\rm target}$ is the softening of the target particle where the field is being summed. This ensures that any interactions between physically-overlapping particles are summed directly with the softening kernel.

# Acknowledgements

We acknowledge code contributions from Ben Keller and Martin Beroiz, and helpful feedback from Elisa Bortolas, Thorsten García, and GitHub user ``herkesg`` during the development of ``pykdgrav``, which were incorporated into ``pytreegrav``.

# References
[![PyPI](https://img.shields.io/pypi/v/pytreegrav)](https://pypi.org/project/pytreegrav)[![Documentation Status](https://readthedocs.org/projects/pytreegrav/badge/?version=latest)](https://pytreegrav.readthedocs.io/en/latest/?badge=latest)

# Introduction
pytreegrav is a package for computing the gravitational potential and/or field of a set of particles. It includes methods for brute-force direction summation and for the fast, approximate Barnes-Hut treecode method. For the Barnes-Hut method we implement an oct-tree as a numba jitclass to achieve much higher peformance than the equivalent pure Python implementation, without writing a single line of C or Cython. Full documentation is available [here](http://pytreegrav.readthedocs.io).

# Installation

```pip install pytreegrav``` or clone the repo and run ```python setup.py install``` from the repo directory.

# Walkthrough
First let's import the stuff we want and generate some particle positions and masses - these would be your particle data for whatever your problem is.


```python
import numpy as np
from pytreegrav import Accel, Potential
```


```python
N = 10**5 # number of particles
x = np.random.rand(N,3) # positions randomly sampled in the unit cube
m = np.repeat(1./N,N) # masses - let the system have unit mass
h = np.repeat(0.01,N) # softening radii - these are optional, assumed 0 if not provided to the frontend functions
```

Now we can use the ``Accel`` and ``Potential`` functions to compute the gravitational field and potential at each particle position:


```python
print(Accel(x,m,h))
print(Potential(x,m,h))
```

    [[-0.1521787   0.2958852  -0.30109005]
     [-0.50678204 -0.37489886 -1.0558666 ]
     [-0.24650087  0.95423467 -0.175074  ]
     ...
     [ 0.87868472 -1.28332176 -0.22718531]
     [-0.41962742  0.32372245 -1.31829084]
     [ 2.45127054  0.38292881  0.05820412]]
    [-2.35518057 -2.19299372 -2.28494218 ... -2.11783337 -2.1653377
     -1.80464695]


By default, pytreegrav will try to make the optimal choice between brute-force and tree methods for speed, but we can also force it to use one method or another. Let's try both and compare their runtimes:


```python
from time import time
t = time()
# tree gravitational acceleration
accel_tree = Accel(x,m,h,method='tree')
print("Tree accel runtime: %gs"%(time() - t)); t = time()

accel_bruteforce = Accel(x,m,h,method='bruteforce')
print("Brute force accel runtime: %gs"%(time() - t)); t = time()

phi_tree = Potential(x,m,h,method='tree')
print("Tree potential runtime: %gs"%(time() - t)); t = time()

phi_bruteforce = Potential(x,m,h,method='bruteforce')
print("Brute force potential runtime: %gs"%(time() - t)); t = time()
```

    Tree accel runtime: 0.927745s
    Brute force accel runtime: 44.1175s
    Tree potential runtime: 0.802386s
    Brute force potential runtime: 20.0234s


As you can see, the tree-based methods can be much faster than the brute-force methods, especially for particle counts exceeding 10^4. Here's an example of how much faster the treecode is when run on a Plummer sphere with a variable number of particles, on a single core of an Intel i9 9900k workstation:
![Benchmark](images/CPU_Time_serial.png)


But there's no free lunch here: the tree methods are approximate. Let's quantify the RMS errors of the stuff we just computed, compared to the exact brute-force solutions:


```python
acc_error = np.sqrt(np.mean(np.sum((accel_tree-accel_bruteforce)**2,axis=1))) # RMS force error
print("RMS force error: ", acc_error)
phi_error = np.std(phi_tree - phi_bruteforce)
print("RMS potential error: ", phi_error)
```

    RMS force error:  0.006739311224338851
    RMS potential error:  0.0003888328578588027


The above errors are typical for default settings: ~1% force error and ~0.1\% potential error. The error in the tree approximation is controlled by the Barnes-Hut opening angle ``theta``, set to 0.7 by default. Smaller ``theta`` gives higher accuracy, but also runs slower:


```python
thetas = 0.1,0.2,0.4,0.8 # different thetas to try
for theta in thetas:
    t = time()    
    accel_tree = Accel(x,m,h,method='tree',theta=theta)
    acc_error = np.sqrt(np.mean(np.sum((accel_tree-accel_bruteforce)**2,axis=1)))
    print("theta=%g Runtime: %gs RMS force error: %g"%(theta, time()-t, acc_error))
```

    theta=0.1 Runtime: 63.1738s RMS force error: 3.78978e-05
    theta=0.2 Runtime: 14.3356s RMS force error: 0.000258755
    theta=0.4 Runtime: 2.91292s RMS force error: 0.00148698
    theta=0.8 Runtime: 0.724668s RMS force error: 0.0105937


Both brute-force and tree-based calculations can be parallelized across all available logical cores via OpenMP, by specifying ``parallel=True``. This can speed things up considerably, with parallel scaling that will vary with your core and particle number:


```python
from time import time
t = time()
# tree gravitational acceleration
accel_tree = Accel(x,m,h,method='tree',parallel=True)
print("Tree accel runtime in parallel: %gs"%(time() - t)); t = time()

accel_bruteforce = Accel(x,m,h,method='bruteforce',parallel=True)
print("Brute force accel runtime in parallel: %gs"%(time() - t)); t = time()

phi_tree = Potential(x,m,h,method='tree',parallel=True)
print("Tree potential runtime in parallel: %gs"%(time() - t)); t = time()

phi_bruteforce = Potential(x,m,h,method='bruteforce',parallel=True)
print("Brute force potential runtime in parallel: %gs"%(time() - t)); t = time()
```

    Tree accel runtime in parallel: 0.222271s
    Brute force accel runtime in parallel: 7.25576s
    Tree potential runtime in parallel: 0.181393s
    Brute force potential runtime in parallel: 5.72611s


# What if I want to evaluate the fields at different points than where the particles are?

We got you covered. The ``Target`` methods do exactly this: you specify separate sets of points for the particle positions and the field evaluation, and everything otherwise works exactly the same (including optional parallelization and choice of solver):


```python
from pytreegrav import AccelTarget, PotentialTarget

# generate a separate set of "target" positions where we want to know the potential and field
N_target = 10**4
x_target = np.random.rand(N_target,3)
h_target = np.repeat(0.01,N_target) # optional "target" softening: this sets a floor on the softening length of all forces/potentials computed

accel_tree = AccelTarget(x_target, x,m, h_target=h_target, h_source=h,method='tree') # we provide the points/masses/softenings we generated before as the "source" particles
accel_bruteforce = AccelTarget(x_target,x,m,h_source=h,method='bruteforce')

acc_error = np.sqrt(np.mean(np.sum((accel_tree-accel_bruteforce)**2,axis=1))) # RMS force error
print("RMS force error: ", acc_error)

phi_tree = PotentialTarget(x_target, x,m, h_target=h_target, h_source=h,method='tree') # we provide the points/masses/softenings we generated before as the "source" particles
phi_bruteforce = PotentialTarget(x_target,x,m,h_target=h_target, h_source=h,method='bruteforce')

phi_error = np.std(phi_tree - phi_bruteforce)
print("RMS potential error: ", phi_error)
```

    RMS force error:  0.006719983300560105
    RMS potential error:  0.0003873676304955059

# Quickstart
pytreegrav is a package for computing the gravitational potential and/or field of a set of particles. It includes methods for brute-force direction summation and for the fast, approximate Barnes-Hut treecode method. For the Barnes-Hut method we implement an oct-tree as a numba jitclass to achieve much higher peformance than the equivalent pure Python implementation.

First let's import the stuff we want and generate some particle positions and masses - these would be your particle data for whatever your problem is.


```python
import numpy as np
from pytreegrav import Accel, Potential
```


```python
N = 10**5 # number of particles
x = np.random.rand(N,3) # positions randomly sampled in the unit cube
m = np.repeat(1./N,N) # masses - let the system have unit mass
h = np.repeat(0.01,N) # softening radii - these are optional, assumed 0 if not provided to the frontend functions
```

Now we can use the ``Accel`` and ``Potential`` functions to compute the gravitational field and potential at each particle position:


```python
print(Accel(x,m,h))
print(Potential(x,m,h))
```

    [[-0.1521787   0.2958852  -0.30109005]
     [-0.50678204 -0.37489886 -1.0558666 ]
     [-0.24650087  0.95423467 -0.175074  ]
     ...
     [ 0.87868472 -1.28332176 -0.22718531]
     [-0.41962742  0.32372245 -1.31829084]
     [ 2.45127054  0.38292881  0.05820412]]
    [-2.35518057 -2.19299372 -2.28494218 ... -2.11783337 -2.1653377
     -1.80464695]


By default, pytreegrav will try to make the optimal choice between brute-force and tree methods for speed, but we can also force it to use one method or another. Let's try both and compare their runtimes:


```python
from time import time
t = time()
# tree gravitational acceleration
accel_tree = Accel(x,m,h,method='tree')
print("Tree accel runtime: %gs"%(time() - t)); t = time()

accel_bruteforce = Accel(x,m,h,method='bruteforce')
print("Brute force accel runtime: %gs"%(time() - t)); t = time()

phi_tree = Potential(x,m,h,method='tree')
print("Tree potential runtime: %gs"%(time() - t)); t = time()

phi_bruteforce = Potential(x,m,h,method='bruteforce')
print("Brute force potential runtime: %gs"%(time() - t)); t = time()
```

    Tree accel runtime: 0.927745s
    Brute force accel runtime: 44.1175s
    Tree potential runtime: 0.802386s
    Brute force potential runtime: 20.0234s


As you can see, the tree-based methods can be much faster than the brute-force methods, especially for particle counts exceeding 10^4. Here's an example of how much faster the treecode is when run on a Plummer sphere with a variable number of particles, on a single core of an Intel i9 9900k workstation:
![Benchmark](./CPU_Time_serial.png)


But there's no free lunch here: the tree methods are approximate. Let's quantify the RMS errors of the stuff we just computed, compared to the exact brute-force solutions:


```python
acc_error = np.sqrt(np.mean(np.sum((accel_tree-accel_bruteforce)**2,axis=1))) # RMS force error
print("RMS force error: ", acc_error)
phi_error = np.std(phi_tree - phi_bruteforce)
print("RMS potential error: ", phi_error)
```

    RMS force error:  0.006739311224338851
    RMS potential error:  0.0003888328578588027


The above errors are typical for default settings: ~1% force error and ~0.1\% potential error. The error in the tree approximation is controlled by the Barnes-Hut opening angle ``theta``, set to 0.7 by default. Smaller ``theta`` gives higher accuracy, but also runs slower:


```python
thetas = 0.1,0.2,0.4,0.8 # different thetas to try
for theta in thetas:
    t = time()    
    accel_tree = Accel(x,m,h,method='tree',theta=theta)
    acc_error = np.sqrt(np.mean(np.sum((accel_tree-accel_bruteforce)**2,axis=1)))
    print("theta=%g Runtime: %gs RMS force error: %g"%(theta, time()-t, acc_error))
```

    theta=0.1 Runtime: 63.1738s RMS force error: 3.78978e-05
    theta=0.2 Runtime: 14.3356s RMS force error: 0.000258755
    theta=0.4 Runtime: 2.91292s RMS force error: 0.00148698
    theta=0.8 Runtime: 0.724668s RMS force error: 0.0105937


Both brute-force and tree-based calculations can be parallelized across all available logical cores via OpenMP, by specifying ``parallel=True``. This can speed things up considerably, with parallel scaling that will vary with your core and particle number:


```python
from time import time
t = time()
# tree gravitational acceleration
accel_tree = Accel(x,m,h,method='tree',parallel=True)
print("Tree accel runtime in parallel: %gs"%(time() - t)); t = time()

accel_bruteforce = Accel(x,m,h,method='bruteforce',parallel=True)
print("Brute force accel runtime in parallel: %gs"%(time() - t)); t = time()

phi_tree = Potential(x,m,h,method='tree',parallel=True)
print("Tree potential runtime in parallel: %gs"%(time() - t)); t = time()

phi_bruteforce = Potential(x,m,h,method='bruteforce',parallel=True)
print("Brute force potential runtime in parallel: %gs"%(time() - t)); t = time()
```

    Tree accel runtime in parallel: 0.222271s
    Brute force accel runtime in parallel: 7.25576s
    Tree potential runtime in parallel: 0.181393s
    Brute force potential runtime in parallel: 5.72611s

 
## What if I want to evaluate the fields at different points than where the particles are?

We got you covered. The ``Target`` methods do exactly this: you specify separate sets of points for the particle positions and the field evaluation, and everything otherwise works exactly the same (including optional parallelization and choice of solver):


```python
from pytreegrav import AccelTarget, PotentialTarget

# generate a separate set of "target" positions where we want to know the potential and field
N_target = 10**4
x_target = np.random.rand(N_target,3)
h_target = np.repeat(0.01,N_target) # optional "target" softening: this sets a floor on the softening length of all forces/potentials computed

accel_tree = AccelTarget(x_target, x,m, h_target=h_target, h_source=h,method='tree') # we provide the points/masses/softenings we generated before as the "source" particles
accel_bruteforce = AccelTarget(x_target,x,m,h_source=h,method='bruteforce')

acc_error = np.sqrt(np.mean(np.sum((accel_tree-accel_bruteforce)**2,axis=1))) # RMS force error
print("RMS force error: ", acc_error)

phi_tree = PotentialTarget(x_target, x,m, h_target=h_target, h_source=h,method='tree') # we provide the points/masses/softenings we generated before as the "source" particles
phi_bruteforce = PotentialTarget(x_target,x,m,h_target=h_target, h_source=h,method='bruteforce')

phi_error = np.std(phi_tree - phi_bruteforce)
print("RMS potential error: ", phi_error)
```

    RMS force error:  0.006719983300560105
    RMS potential error:  0.0003873676304955059

