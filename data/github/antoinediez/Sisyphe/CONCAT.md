---
title: '`SiSyPHE`: A Python package for the Simulation of Systems of interacting mean-field Particles with High Efficiency'
tags:
  - Python
  - GPU
  - particles
  - mean-field
  - self-organization
  - swarming
authors:
  - name: Antoine Diez
    affiliation: 1
affiliations:
 - name: Department of Mathematics, Imperial College London, South Kensington Campus, London, SW7 2AZ, UK
   index: 1
date: 27 June 2021
bibliography: doc/biblio_sisyphe.bib

---

# Summary

Over the past decades, the study of systems of particles has become an important 
part of many research areas, from theoretical physics to applied biology and 
computational mathematics. One of the main motivations in mathematical biology is the modelling of large
animal societies and the emergence of complex patterns from simple behavioral rules, e.g., flocks of birds, fish schools, ant colonies, etc. In the microscopic
world, particle systems are used to model a wide range of phenomena, 
from the collective motion of spermatozoa to the anarchical development of 
cancer cells. Within this perspective, there are at least three important reasons to conduct large scale computer simulations of particle systems. First, numerical experiments are essential to 
calibrate the models and test the influence of each parameter in a controlled 
environment. For instance, the renowned Vicsek model [@vicsek_novel_1995]
is a minimal model of *flocking*, which exhibits a complex behavior, studied numerically in particular in [@chate_collective_2008]. Secondly, particle simulations are used to check the validity of *macroscopic* models that describe the statistical behavior of particle systems.
These models are usually based on partial differential equations (PDE) derived 
using phenomenological considerations that are often difficult to justify 
mathematically [@degond_continuum_2008; @dimarco_self-alignment_2016; @degond_bulk_2021]. Finally, inspired by models in biology, there is an ever growing
literature on the design of algorithms based on the simulation of *artificial*
particle systems to solve tough optimization problems [@kennedy_particle_1995; @pinnau_consensus-based_2017; @totzeck_trends_2021; @grassi_particle_2020] and to construct new more efficient Markov Chain Monte Carlo methods [@del_moral_measure-valued_1998; @del_moral_mean_2013; @doucet_sequential_2001; @cappe_population_2004; @clarte_collective_2021]. The simulation of systems of particles is also at the core of molecular dynamics [@leihmkuhler], although the present library is not specifically written for this purpose. The `SiSyPHE` library builds on recent advances in hardware and software 
for the efficient simulation of large scale interacting *mean-field* particle systems, 
both on the GPU and on the CPU. The versatile object-oriented Python interface of the library is designed for the simulation and comparison of new and classical many-particle models of collective dynamics in mathematics and active matter physics, enabling ambitious numerical experiments and leading to novel conjectures and results.

# Statement of need

A major difficulty in the simulation of systems of particles is the high computational cost, typically quadratic in the number of particles, which prevents large scale experiments. The implementation of `SiSyPHE` is based on recent libraries 
originally developed for machine learning purposes to significantly accelerate 
tensor (array) computations, namely the `PyTorch` package [@paszke_pytorch_2019] and the `KeOps` library [@charlier_kernel_2021]. On a GPU, the `SiSyPHE` library speeds up both traditional Python and low-level implementations by one to three orders 
of magnitude for systems with up to several millions of particles. 

In addition, to the best of our knowledge, only model-specific packages such as @motsch_vicsek_microflat_2016 are available. The `SiSyPHE` library includes, within a common framework, the implementation of many classical models and their variants as well as recent models for which no implementation was previously available. All the models detailed in the Example gallery of the documentation are directly taken from the literature on collective dynamics in mathematics and active matter physics. Moreover, the `SiSyPHE` library is designed in such a way that new custom models can easily be added in order to facilite the study and comparison of models from a research perspective. 

The development of the `SiSyPHE` library was initially motivated by the study of *body-oriented particles* [@giacomin_alignment_2019]. 
The (formal) derivation of a macroscopic PDE model from the particle system has lead to a novel conjecture 
which postulates the existence of a class of so-called *bulk topological states* in [@degond_bulk_2021]. The quantitative comparison
between this theoretical prediction and the numerical simulation of the particle system in a suitable regime (with more than
$10^6$ particles) has confirmed the existence of these new states of matter. The study of their physical properties
which are observed in the numerical experiments but not readily explained by the PDE model is an ongoing work.

# A typical example

A typical model that is implemented in the `SiSyPHE` library is the variant of the Vicsek model
introduced by @degond_continuum_2008 and defined by the system of $2N$ Stratonovich Stochastic Differential Equations
\begin{equation}\label{eq:sde}
\mathrm{d}X^i_t = c_0 V^i_t \mathrm{d}t, \quad \mathrm{d}V^i_t = \sigma \mathsf{P}(V^i_t)\circ(J^i_t\mathrm{d}t+\mathrm{d}B^i_t),
\end{equation}
where the position at time $t$ of a particle indexed by $i\in\{1,\ldots,N\}$ is a vector $X^i_t\in\mathbb{R}^d$ and its orientation (or velocity) is a unit vector $V^i_t\in\mathbb{R}^d$ with $|V^i_t|=1$. The coefficient $c_0>0$ is the speed
of the particles (assumed to be constant), the matrix $\mathsf{P}(V^i_t)= I_d - V^i_t\otimes V^i_t$ is the orthogonal projection matrix on the plane orthgonal to $V^i_t$, $(B^i_t)^{}_t$ is an independent Brownian motion, and $\sigma>0$ is a diffusion coefficient which models the level of noise. 
The quantity $J^i_t\in\mathbb{R}^d$ is called a *target*; it is the orientation that particle $i$ is trying to adopt. 
In the Vicsek model introduced by @degond_continuum_2008, 
\begin{equation}\label{eq:target}
J^i_t = \frac{\sum_{j=1}^N K(|X^j_t-X^i_t|)V^j_t}{\big|\sum_{j=1}^N K(|X^j_t-X^i_t|)V^j_t\big|}, 
\end{equation}
where the *kernel* $K:[0,+\infty)\to[0,+\infty)$ is a smooth nonnegative
function vanishing at infinity which models the visual perception of the particles; 
in the Vicsek model, the vision of the particles depends on the distance between them. 
With the target given by \autoref{eq:target}, each particle tries to adopt the average orientation of its neighbors, which is a typical *flocking* behavior. 

On a computer, the time-continuous system given by \autoref{eq:sde} needs to be discretized first. For the Vicsek model, a natural discretization method is the (geometric) Euler-Maruyama scheme [@kloeden; @piggott]. In general, the discretization method depends on the model considered as illustrated in the Example gallery. Then, at each time step, 
the most expensive operation is the computation of the target given by \autoref{eq:target}, which requires $\mathcal{O}(N)$
operations for each of the $N$ particles. The total simulation cost is thus $\mathcal{O}(N^2T)$ where $T$ is the
total number of iterations. Within the framework of the `KeOps` library on which `SiSyPHE` is based, 
the computation of the target \autoref{eq:target} is called a *kernel operation*, which is efficiently carried out
using a *symbolic* definition of the $N\times N$ interaction matrix whose $(i,j)$-entry is $K(|X^j_t-X^i_t|)$. The computation of the target is then understood as a symbolic matrix-vector product between the interaction matrix and the vector of orientations.  


# Acknowledgements

The development of this library would not have been possible without the help of Jean Feydy, 
his constant support and precious advice. This project was initiated by Pierre Degond and 
has grown out of many discussions with him.

# References
# Simulation of Systems of interacting mean-field Particles with High Efficiency

<p align="center">
<img src="./doc/_static/ball.png" alt="logo" width="224"/>
</p>

Please visit the [website](https://sisyphe.readthedocs.io/en/latest/) for a full documentation.

------------------------------------------------------------------------------------------------


The SiSyPHE library builds on recent advances in hardware and software for the efficient simulation of **large scale interacting particle systems**, both on the **GPU** and on the CPU. The implementation is based on recent libraries originally developed for machine learning purposes to significantly accelerate tensor (array) computations, namely the [PyTorch](https://github.com/pytorch/pytorch) package and the [KeOps](https://www.kernel-operations.io/keops/index.html) library. The **versatile object-oriented Python interface** is well suited to the comparison of new and classical many-particle models, enabling ambitious numerical experiments and leading to novel conjectures. The SiSyPHE library speeds up both traditional Python and low-level implementations by **one to three orders of magnitude** for systems with up to **several millions** of particles. 

<p align="center">
<img src="./doc/_static/mill.gif" alt="mill">
</p>

## Citation 

If you use SiSyPHE in a research paper, please cite the [JOSS publication](https://joss.theoj.org/papers/10.21105/joss.03653#) :

```
@article{Diez2021,
  doi = {10.21105/joss.03653},
  url = {https://doi.org/10.21105/joss.03653},
  year = {2021},
  publisher = {The Open Journal},
  volume = {6},
  number = {65},
  pages = {3653},
  author = {Antoine Diez},
  title = {`SiSyPHE`: A Python package for the Simulation of Systems of interacting mean-field Particles with High Efficiency},
  journal = {Journal of Open Source Software}
}
```

Diez, A., (2021). SiSyPHE: A Python package for the Simulation of Systems of interacting mean-field Particles with High Efficiency. Journal of Open Source Software, 6(65), 3653, https://doi.org/10.21105/joss.03653

[![DOI](https://joss.theoj.org/papers/10.21105/joss.03653/status.svg)](https://doi.org/10.21105/joss.03653)

## Installation 

### Requirements

- **Python 3** with packages **NumPy** and **SciPy** 
- **PyTorch** : version>= 1.5
- **PyKeops** : version>= 1.5

### Using pip

In a terminal, type:

```
pip install sisyphe
```
    
### On Google Colab

The easiest way to get a working version of SiSyPHE is to use the free virtual machines provided by [Google Colab](https://colab.research.google.com).

1. On a new Colab notebook, navigate to Edit→Notebook Settings and select GPU from the Hardware Accelerator drop-down.

2. Install PyKeops with the Colab specifications **first** by typing

```    
!pip install pykeops[colab]
```

3. Install SiSyPHE by typing 

```
!pip install sisyphe
```    

### Testing the installation

In a Python terminal, type 

```python
import sisyphe
sisyphe.test_sisyphe()
```    

<p align="center">
<img src="./doc/_static/band.gif" alt="band">
</p>   

## Contributing

Contributions to make SiSyPHE grow are warmly welcome! Examples of possible (and ongoing) developments include the following. 

* The implementation of new models.

* The implementation of more complex boundary conditions and of models on *non-flat* manifolds. 

* An improved visualization method (currently only basic visualization functions relying on [Matplotlib](https://matplotlib.org/) are implemented). 

Contributions can be made by opening an issue on the GitHub repository, via a pull request or by contacting directly the author.  


## Author

- [Antoine Diez](https://antoinediez.gitlab.io), Imperial College London 

### Acknowledgments

The development of this library would not have been possible without the help of [Jean Feydy](https://www.jeanfeydy.com/), his constant support and precious advice. This project was initiated by [Pierre Degond](https://sites.google.com/site/degond/) and has grown out of many discussions with him. 
# How to contribute

Thank you for your interest in SiSyPHE, contributions to make the library grow are warmly welcome! 

## Examples of contributions

Examples of possible (and ongoing) developments to which you are welcome to contribute include the following. 

* The implementation of new models.

* The implementation of more complex boundary conditions and of models on *non-flat* manifolds. 

* An improved visualization method (currently only basic visualization functions relying on [Matplotlib](https://matplotlib.org/) are implemented). 

## Submitting changes 

If you find a bug or an unexpected behaviour, please open an issue on Github. To add a new feature to the library, please send a pull request with a clear description of what you have done. If you are working on a model which you would like to be included in SiSyPHE, I will be happy to discuss it beforehand either directly on Github or by email first. ========
Examples
========

Some examples for various models. =========
Tutorials
=========

Basic features. ===============================
Installation
===============================


Requirements
================

- **Python 3** with packages **NumPy** and **SciPy** 
- **PyTorch** : version>= 1.5
- **PyKeops** : version>= 1.5

.. note::
    In order to have a working version of PyKeOps, it may be necessary to install and upgrade Cmake to Cmake>=3.18 as discussed `here <https://github.com/getkeops/keops/issues/142>`_.

Installation
============= 

Using pip
------------

In a terminal, type:

.. code-block:: bash

    pip install sisyphe
    
On Google Colab
-------------------

The easiest way to get a working version of SiSyPHE is to use the free virtual machines provided by `Google Colab <https://colab.research.google.com>`_.

1. On a new Colab notebook, navigate to Edit→Notebook Settings and select GPU from the Hardware Accelerator drop-down.

2. Install PyKeops with the Colab specifications **first** by typing
    
    .. code-block:: bash

        !pip install pykeops[colab]
    
3. Install SiSyPHE by typing 
    
    .. code-block:: bash

        !pip install sisyphe
    
From source
-------------------

Alternatively, you can clone the `git repository <https://github.com/antoinediez/Sisyphe>`_ at a location of your choice. 


Testing the installation
============================

The following test function will check the configuration and run the simulation of a system of body-oriented particles (see the example gallery). This simulation uses the main features of the SiSyPHE library: a complex system, nontrivial boundary conditions, a sampling-based interaction mechanism, the blocksparse reduction method and a nontrivial initial condition. Moreover, this model is `theoretically well-understood <https://arxiv.org/abs/2101.10864>`_ which provides a theoretical baseline to check the accuracy of the output of the simulation. 

.. warning::
    This function is mainly intended to be runned on a GPU. Running this function on a CPU will take a long time! See below for a quick testing procedure. 

In a Python terminal, type 

.. code-block:: python

    import sisyphe
    sisyphe.test_sisyphe()
    
On a fresh environment, it should return

.. code-block:: console

    Welcome! This test function will create a system of body-oriented particles in a ``milling configuration'' (cf. the example gallery). The test will be considered as successful if the computed milling speed is within a 5% relative error range around the theoretical value.

     Running test, this may take a few minutes...

     Check configuration... 
    [pyKeOps] Initializing build folder for dtype=float32 and lang=torch in /root/.cache/pykeops-1.5-cpython-37 ... done.
    [pyKeOps] Compiling libKeOpstorch180bebcc11 in /root/.cache/pykeops-1.5-cpython-37:
           formula: Sum_Reduction(SqNorm2(x - y),1)
           aliases: x = Vi(0,3); y = Vj(1,3); 
           dtype  : float32
    ... 
    [pyKeOps] Compiling pybind11 template libKeOps_template_574e4b20be in /root/.cache/pykeops-1.5-cpython-37 ... done.
    Done.

    pyKeOps with torch bindings is working!

    Done.

     Sample an initial condition... 
    Done.

     Create a model... 
    Done.

     Run the simulation... 
    [pyKeOps] Compiling libKeOpstorch269aaf150e in /root/.cache/pykeops-1.5-cpython-37:
           formula: Sum_Reduction((Step((Var(5,1,2) - Sum(Square((((Var(0,3,1) - Var(1,3,0)) + (Step(((Minus(Var(2,3,2)) / Var(3,1,2)) - (Var(0,3,1) - Var(1,3,0)))) * Var(2,3,2))) - (Step(((Var(0,3,1) - Var(1,3,0)) - (Var(2,3,2) / Var(4,1,2)))) * Var(2,3,2))))))) * Var(6,16,1)),0)
           aliases: Var(0,3,1); Var(1,3,0); Var(2,3,2); Var(3,1,2); Var(4,1,2); Var(5,1,2); Var(6,16,1); 
           dtype  : float32
    ... 
    Done.
    Progress:100%Done.

     Check the result... 
    Done.

     SiSyPHE is working!    

    
The core functionalities of the library are automatically and continuously tested through a GitHub workflow based on the module :mod:`sisyphe.test.quick_test`. The testing functions include basic computations on simple examples (computation of simple local averages in various situations) and small scales simulations. Note that unlike the function :meth:`sisyphe.test_sisyphe()`, these testing functions do not check the accuracy of the output of the simulations but only check that the code runs without errors. It is possible to use the `Pytest package <https://docs.pytest.org/en/6.2.x/>`_ to run these tests manually: on a Python terminal, type

.. code-block:: python

    import pytest
    from sisyphe.test import quick_test
    retcode = pytest.main([quick_test.__file__,])


    Benchmarks
###########

We compare 

* a `Fortran implementation <https://github.com/smotsch/Vicsek_microFlat>`_ due to Sébastien Motsch using the Verlet list method with double precision float numbers, run on the `NextGen <http://sysnews.ma.ic.ac.uk/compute-cluster/>`_ compute cluster at Imperial College London. 
* the CPU version of SiSyPHE with double precision tensors (float64) and the blocksparse-reduction method (BSR), run on an Intel MacBook Pro (2GHz Intel Core i5 processor with 8 Go of memory) ; 
* the GPU version of SiSyPHE with single precision tensors (float32) with and without the :ref:`blocksparse-reduction method <tutobsr>` (BSR), run on a `GPU cluster <http://sysnews.ma.ic.ac.uk/GPU-computing/>`_ at Imperial College London using an nVidia GTX 2080 Ti GPU ;
* the GPU version of SiSyPHE with double precision tensors (float64) with and without the :ref:`blocksparse-reduction method <tutobsr>` (BSR), run on a `GPU cluster <http://sysnews.ma.ic.ac.uk/GPU-computing/>`_ at Imperial College London using an nVidia GTX 2080 Ti GPU.

.. note::

    The choice of single or double precision tensors is automatically made according to the floating point precision of the input tensors. 

We run a :class:`Vicsek <sisyphe.models.Vicsek>` model in a square periodic box with fixed parameters :math:`L=100`, :math:`\nu=5`, :math:`\sigma=1`, :math:`c=R`, :math:`dt=0.01` and various choices of :math:`N` and :math:`R`. The simulation is run for 10 units of time (i.e. 1000 iterations). 

For each value of :math:`N`, three values of :math:`R` are tested which correspond to a dilute regime, a moderate regime and a dense mean-field regime. When the particles are uniformly scattered in the box, the average number of neighbours in the three regimes is respectively :math:`\sim3`, :math:`\sim30` and :math:`\sim300`. The regime has a strong effect on the efficiency of the Verlet list and blocksparse reduction methods. 

In the table below, the total computation times are given in seconds except when they exceed 12 hours. The best times are indicated in **bold** and the worst times in *italic*. 

+----------+---------+-----------+-------------+-------------+-------------+-------------+-------------+
|          |         | | Fortran | | sisyphe64 | | sisyphe32 | | sisyphe32 | | sisyphe64 | | sisyphe64 |
|          |         | |         | | CPU BSR   | | GPU       | | GPU BSR   | | GPU       | | GPU BSR   |
+==========+=========+===========+=============+=============+=============+=============+=============+
|          | R = 1   |    19s    |   *26s*     |             |     3.3s    |             |     3.6s    |
|          +---------+-----------+-------------+             +-------------+             +-------------+
| N = 10k  | R = 3   |   *59s*   |    29s      |             |     3.4s    |             |     3.9s    |
|          +---------+-----------+-------------+  **2.1s**   +-------------+     13s     +-------------+
|          | R = 10  |  *494s*   |    69s      |             |     3.4s    |             |     7.9s    |
+----------+---------+-----------+-------------+-------------+-------------+-------------+-------------+
|          | R = 0.3 |   309s    |  *323s*     |             |   **4.3s**  |             |     9.3s    |
|          +---------+-----------+-------------+             +-------------+             +-------------+
| N = 100k | R = 1   |  *1522s*  |   384s      |     29s     |   **4.5s**  |    973s     |      11s    |
|          +---------+-----------+-------------+             +-------------+             +-------------+
|          | R = 3   |  *3286s*  |   796s      |             |   **4.9s**  |             |      28s    |
+----------+---------+-----------+-------------+-------------+-------------+-------------+-------------+
|          | R = 0.1 |  *>12h*   |   6711s     |             |    **22s**  |             |     120s    |
|          +---------+-----------+-------------+             +-------------+             +-------------+
| N = 1M   | R = 0.3 |  *>12h*   |   6992s     |    2738s    |    **23s**  |  *>12h*     |     135s    |
|          +---------+-----------+-------------+             +-------------+             +-------------+
|          | R = 1   |  *>12h*   |   9245s     |             |    **26s**  |             |     194s    |
+----------+---------+-----------+-------------+-------------+-------------+-------------+-------------+

The GPU implementation is at least 5 times faster in the dilute regime and outperform the other methods by three orders of magnitude in the mean-field regime with large values of :math:`N`. Without the block-sparse reduction method, the GPU implementation does suffer from the quadratic complexity. The block-sparse reduction method is less sensitive to the density of particles than the traditional Verlet list method. It comes at the price of a higher memory cost (see :ref:`tutobsr`) but allows to run large scale simulations even on a CPU. On a CPU, the performances are not critically affected by the precision of the floating-point format so only simulations with double precision tensors are shown.

.. note::

    The parameters of the blocksparse reduction method are chosen using the method :meth:`best_blocksparse_parameters() <sisyphe.particles.Particles.best_blocksparse_parameters>` (see :ref:`tutobsr`). As a guideline, in this example and on the GPU, the best number of cells is respectively :math:`\sim 15^2`, :math:`\sim 42^2` and :math:`\sim 70^2` for :math:`N=10^4`, :math:`N=10^5` and :math:`N=10^6` with sisyphe32 and  :math:`\sim 25^2`, :math:`\sim 60^2` and :math:`\sim 110^2` for sisyphe64. Note that the number of cells will be automatically capped at :math:`(L/R)^2`. For the last case (sisyphe64 with :math:`R=1`), the maximal number of cells depends on the memory available. 
    
As an example, the case sisyphe32 GPU BSR with :math:`N=10^5` and :math:`R=1` is obtained with the following script. 

.. code-block:: python

    import time
    import math
    import torch
    import sisyphe.models as models
    from sisyphe.display import save

    dtype = torch.cuda.FloatTensor

    # Replace by the following line for double precision tensors.
    # dtype = torch.cuda.DoubleTensor

    N = 100000
    L = 100.
    dt = .01

    nu = 5.
    sigma = 1.

    R = 1.
    c = R

    # The type of the initial condition will determine the floating point precision.

    pos = L*torch.rand((N,2)).type(dtype)
    vel = torch.randn(N,2).type(dtype)
    vel = vel/torch.norm(vel,dim=1).reshape((N,1))

    simu = models.Vicsek(
        pos = pos,
        vel = vel,
        v = R,
        sigma = sigma,
        nu = nu,
        interaction_radius = R,
        box_size = L,
        dt = dt,
        block_sparse_reduction = True,
        number_of_cells = 42**2)

    simu.__next__() # GPU warmup

    s = time.time()
    data = save(simu, [10.], [], [])
    e = time.time()

    print(e-s)



===============================
Background and motivation
===============================

.. image:: _static/corridor.gif
    :scale: 100% 
    :alt: Corridor
    :align: center

Scope of application
=======================


Over the past decades, the study of systems of particles has become an important part of many research areas, from theoretical physics to applied biology and computational mathematics. Some current research trends include the following. 

* Since the 80's, there is a joint effort from biologists, physicists, computer graphics scientists and more recently mathematicians to propose accurate models for large **self-organized** animal societies. While the motivations of the different communities are largely independent, a commmon goal is to be able to explain and reproduce the **emergence of complex patterns from simple interaction rules**. Typical examples of complex systems include `flocks of birds <https://en.wikipedia.org/wiki/Flock_(birds)>`_, `fish schools <https://en.wikipedia.org/wiki/Shoaling_and_schooling>`_, `herds of mammals <https://www.youtube.com/watch?v=5IFLz4CETj4>`_ or `ant colonies <https://en.wikipedia.org/wiki/Ant_colony>`_. It seems reasonable to consider that none of these animals has an accute consciousness of the whole system organization but rather follows simple behavioral patterns: stay close to the group, do not collide with another individual, copy the attitude of the close neighbours etc. These models are often called **swarming models** or **collective dynamics** models :cite:`vicsek_collective_2012, degond_mathematical_2018, albi_vehicular_2019, naldi_particle_2010`.

* The microscopic world if full of complex systems made of many simple entities. Colonies of bacteria :cite:`wensink_meso-scale_2012` are a natural example reminiscent of the macroscopic animal societies described above. But complex patterns can also emerge from **systems of non-living particles**. Our body is actually a very complex self-organized assembly of cells which dictate every aspect of our life: our brain works thanks to the communications between billions of neurons :cite:`fournier_toy_2016, andreis_mckeanvlasov_2018`, the reproduction is based on the competition between millions of spermatozoa :cite:`creppy_turbulence_2015` and sometimes, death is sadly caused by the relentless division of cancer cells. Unlike the previous models, the communication between the particles cannot be based on their visual perception, but rather on chemical agents, on physical (geometrical) constraints etc. 

* On a completely different side, there is an ever growing number of methods in computational mathematics which are based on the simulation of systems of particles. The pioneering `Particle Swarm Optimization method <https://en.wikipedia.org/wiki/Particle_swarm_optimization>`_ :cite:`kennedy_particle_1995` has shown how biologically-inspired artificial systems of particles can be used to solve tough optimization problems. Following these ideas, other **Swarm Intelligence** optimization algorithms have been proposed up to very recently :cite:`pinnau_consensus-based_2017, totzeck_numerical_2018, grassi_particle_2020, totzeck_trends_2021`. Since the beginning of the 2000's, particle systems are also at the core of `filtering <https://en.wikipedia.org/wiki/Particle_filter>`_ :cite:`crisan_survey_2002, doucet_sequential_2001, kantas_overview_2009, del_moral_feynman-kac_2004, del_moral_mean_2013` and `sampling <https://en.wikipedia.org/wiki/Monte_Carlo_method>`_ methods :cite:`cappe_population_2004, clarte_collective_2021`. Recently, a particle-based interpretation :cite:`mei_mean_2018, chizat_global_2018, rotskoff_trainability_2019, de_bortoli_quantitative_2020, sirignano_mean_2020` of the training task of neural networks has lead to new theoretical convergence results.

Why do we need to simulate particle systems?
==================================================

Beyond the self-explanatory applications for particle-based algorithms in computational mathematics, the simulation of systems of particles is also a crucial **modelling tool** in Physics and Biology. 

* On the one hand, field experiments are useful to collect data and trigger new modelling ideas. On the other hand, numerical simulations become necessary to test these ideas and to calibrate the models. **Numerical experiments** can be conducted to identify which mechanisms are able to produce a specific phenomena in a **controlled environment** :cite:`chate_collective_2008`. 

* On a more theoretical side, the direct mathematical study of particle systems can rapidly become incredibly difficult. Inspired by the `kinetic theory of gases <https://en.wikipedia.org/wiki/Kinetic_theory_of_gases>`_, many **mesoscopic** and **macroscopic** models have been proposed to model the average statistical behavior of particle systems rather than the individual motion of each particle :cite:`toner_flocks_1998, degond_continuum_2008, naldi_particle_2010`. Within this framework, a particle system is rather seen as a **fluid** and it is decribed by `Partial Differential Equations <https://en.wikipedia.org/wiki/Partial_differential_equation>`_ (PDE) reminiscent from the theory of `fluid dynamics <https://en.wikipedia.org/wiki/Fluid_dynamics>`_. PDE models are more easily theoretically and numerically tractable. However, **cheking the validity** of these models is not always easy and is sometimes only postulated based on phenomenological considerations. In order to design good models, it is often necessary to go back-and-forth between the PDE models and the numerical simulation of the underlying particle systems. 

The development of the SiSyPHE library was initially motivated by the study of :class:`body-oriented particles <sisyphe.particles.BOParticles>` :cite:`giacomin_alignment_2019`. The (formal) derivation of a macroscopic PDE model from the particle system has lead to a novel conjecture which postulates the existence of a class of so-called **bulk topological states** :cite:`degond_bulk_2021`. The quantitative comparison between this theoretical prediction and the `numerical simulation of the particle system <https://figshare.com/projects/Bulk_topological_states_in_a_new_collective_dynamics_model/96491>`_ in a suitable regime (with more than :math:`10^6` particles) has confirmed the existence of these new states of matter. The study of their physical properties which are observed in the numerical experiments but not readily explained by the PDE model is an ongoing work.


Mean-field particle systems
==============================================

Currently, the models implemented in the SiSyPHE library belong to the family of **mean-field models**. It means that the motion of each particle is influenced by the average behavior of the whole system. The `Vicsek model <https://en.wikipedia.org/wiki/Vicsek_model>`_ :cite:`vicsek_novel_1995, degond_continuum_2008` and the Cucker-Smale model :cite:`cucker_mathematics_2007, ha_emergence_2009, naldi_particle_2010` are two popular examples of mean-field models where each particle tries to move in the average direction of motion of its neighbours (it produces a so-called *flocking* behavior). 


The mathematical point of view   
--------------------------------

From a mathematical point of view, a mean-field particle system with :math:`N` particles is defined as a Markov process in :math:`(\mathbb{R}^d)^N` with a generator :math:`\mathcal{L}_N` whose action on a test function :math:`\varphi_N` is of the form:

.. math::

    \mathcal{L}_N\varphi_N(x^1,\ldots,x^N) = \sum_{i=1}^N L_{\mu_{\mathbf{x}^N}}\diamond_i \varphi_N (x^1,\ldots,x^N),
    
where :math:`\mathbf{x}^N = (x^1,\ldots,x^N)\in (\mathbb{R}^d)^N` and

.. math::

    \mu_{\mathbf{x}^N} := \frac{1}{N}\sum_{i=1}^N \delta_{x^i},

is the so-called **empirical measure**. The operator :math:`L_{\mu_{\mathbf{x}_N}}` depends on the empirical measure and acts on one-variable test functions on :math:`\mathbb{R}^d`. Given an operator :math:`L` and a test function :math:`\varphi_N`, the notation :math:`L\diamond_i \varphi_N` denotes the function 

.. math::

    L\diamond_i \varphi_N : (x^1,\ldots,x^N) \in (\mathbb{R}^d)^N \mapsto L[x\mapsto\varphi_N(x_1,\ldots,x^{i-1},x,x^{i+1},\ldots,x^N)](x^i)\in \mathbb{R}.
    
    
The dynamics of the mean-field system depends on the operator :math:`L_{\mu_{\mathbf{x}_N}}` which is typically either a `diffusion operator <https://en.wikipedia.org/wiki/Diffusion_process>`_ :cite:`degond_continuum_2008, ha_emergence_2009` or a `jump operator <https://en.wikipedia.org/wiki/Jump_process>`_ :cite:`dimarco_self-alignment_2016, andreis_mckeanvlasov_2018`. In the first case, the mean-field particle systems can be alternatively defined by a system of :math:`N` coupled `Stochastic Differential Equations <https://en.wikipedia.org/wiki/Stochastic_differential_equation>`_ of the form: 

.. math::

    \mathrm{d}X^i_t = b(X^i_t,\mu_{\mathcal{X}^N_t})\mathrm{d}t + \sigma(X^i_t,\mu_{\mathcal{X}^N_t})\mathrm{d} B^i_t,\quad i\in\{1,\ldots,N\},

.. math::

    \mu_{\mathcal{X}^N_t} := \frac{1}{N}\sum_{i=1}^N \delta_{X^i_t},
    
where :math:`(B^i_t)_t` are :math:`N` independent Brownian motions and the coefficients :math:`b` and :math:`\sigma` are called the *drift* and *diffusion* coefficients. In most cases, these coefficients are *linear* or *quasi-linear* which means that they are of the form 

.. math::

    b(X^i_t,\mu_{\mathcal{X}^N_t}) = \tilde{b}{\left(X^i_t, \frac{1}{N}\sum_{j=1}^N K(X^i_t,X^j_t)\right)}, 

for two given functions :math:`\tilde{b}:\mathbb{R}^d\times\mathbb{R}^n\to\mathbb{R}^d` and :math:`K:\mathbb{R}^d\times\mathbb{R}^d\to \mathbb{R}^n`. 

When the particles are initially statistically independent, it can be shown that when :math:`N\to+\infty`, each particle :math:`X^i_t` converges towards an independent copy of the solution of the so-called **McKean-Vlasov** diffusion process :cite:`mckean_propagation_1969, sznitman_topics_1991, meleard_asymptotic_1996` defined by the Stochastic Differential Equation 

.. math::

    \mathrm{d}\overline{X}_t = b(\overline{X}_t,f_t)\mathrm{d}t + \sigma(\overline{X}_t,f_t)\mathrm{d} B_t,
    
where :math:`(B_t)_t` is a Brownian motion and :math:`f_t` is the law of the process :math:`\overline{X}_t`. It satisfies the **Fokker-Planck** Partial Differential Equation 

.. math::

    \partial_t f_t = - \nabla\cdot(b(x,f_t)f_t) + \frac{1}{2}\sum_{i,j=1}^N \partial_{x_i}\partial_{x_j}(a_{ij}(x,f_t)f_t),
    
with :math:`x=(x_1,\ldots,x_d)\in \mathbb{R}^d` and :math:`a=\sigma\sigma^\mathrm{T}`. This phenomenon is called **propagation of chaos**, following the terminology introduced by Kac in the 50's :cite:`kac_foundations_1956, mckean_propagation_1969, sznitman_topics_1991, meleard_asymptotic_1996, hauray_kacs_2014`. 


.. note::
    
    A popular example of mean-field particle system is the (stochastic) Cucker-Smale model :cite:`cucker_mathematics_2007, ha_emergence_2009, naldi_particle_2010`. Each particle is defined by its position :math:`X^i_t \in\mathbb{R}^d` and its velocity :math:`V^i_t\in\mathbb{R}^d` which evolve according to the system of :math:`2N` Stochastic Differential Equations: 
    
    .. math::
    
        \mathrm{d}X^i_t = V^i_t\mathrm{d}t,\quad \mathrm{d}V^i_t = \frac{1}{N}\sum_{i=1}^N K(|X^j_t-X^i_t|)(V^j_t-V^i_t)\mathrm{d}t + \mathrm{d}B^i_t,
        
    where :math:`K:[0,+\infty)\to[0,+\infty)` is an **observation kernel** which models the visual perception of the particles. In this example, the function :math:`K` is a smooth function vanishing at infinity and the communication between the particles is based on the distance between them. The motion of the particles follows the Newton's laws of motion with an additional stochastic term. The term :math:`V^j_t-V^i_t` is a relaxation force (for a quadratic potential) which tends to align the velocities of particle :math:`i` and particle :math:`j`. 
        


Simulating mean-field particle systems
------------------------------------------------------

On a computer, the above examples which are time-continuous needs to be discretized, using for instance `one of the classical numerical schemes <https://en.wikipedia.org/wiki/Euler%E2%80%93Maruyama_method>`_ for Stochastic Differential Equations. Then the difficulty lies in the evaluation of the empirical measure **at each time-step** of the numerical scheme. 

In the example of the Cucker-Smale model, the force exerted on a particle is the sum of :math:`N` small relaxation forces of order :math:`1/N`. The total number of operations required is thus of order :math:`\mathcal{O}(N)` **for each particle**. Since there are :math:`N` particles, the total time complexity of the algorithm is thus :math:`\mathcal{O}(N^2T)` where :math:`T` is the total number of iterations of the numerical scheme. 

This **quadratic cost** is the main bottleneck in the simulation of mean-field particle systems. As explained in `the documentation of the KeOps library <https://www.kernel-operations.io/keops/introduction/why_using_keops.html>`_, the evaluation of the :math:`N` forces at time :math:`t`

.. math::
    
    F_i(t) = \frac{1}{N}\sum_{i=1}^N K(|X^j_t-X^i_t|)(V^j_t-V^i_t), \quad i\in\{1,\ldots,N\},
    
is called a **kernel operation** and can be understood as a discrete convolution operation (or matrix-vector product) between the matrix of distances :math:`(K_{ij})_{i,j\in\{1\,\ldots,N\}}` where :math:`K_{ij} = K(|X^j_t-X^i_t|)` and the vector of velocities :math:`(V^j_t-V^i_t)_{j\in\{1,\ldots, N\}}`. When :math:`N` is large (say :math:`N>10^4`), such operation is too costly even for array-based programming languages such as Matlab: the :math:`N\times N` kernel matrix :math:`(K_{ij})_{i,j}` would simply not fit into the memory. With lower level languages (Fortran, C), this operation can be implemented more efficiently with two nested loops but with a significantly higher global coding effort and with less versatility. 

Over the past decades, several workarounds have been proposed. Popular methods include 

* the `low-rank decomposition <https://en.wikipedia.org/wiki/Low-rank_matrix_approximations>`_ of the kernel matrix, 

* the fast-multipole methods :cite:`greengard_fast_1987` which to treat differently short- and long-range interactions, 

* the `Verlet list method <https://en.wikipedia.org/wiki/Verlet_list>`_ which is based on a grid decompostion of the spatial domain to reduce the problem to only short-range interactions between subsets of the particle system, 

* the Random Batch Method :cite:`jin_random_2019` which is based on a stochastic approximation where only interactions between randomly sampled subsets (*batches*) of the particle system are computed. 

All these methods require an significant amount of work, either in terms of code or to justify the approximation procedures. 


The SiSyPHE library
========================

The present implementation is based on recent libraries originally developed for machine learning purposes to significantly accelerate such tensor (array) computations, namely the `PyTorch <https://github.com/pytorch/pytorch>`_ package and the `KeOps <https://www.kernel-operations.io/keops/index.html>`_ library :cite:`charlier_kernel_2021`. Using the KeOps framework, the kernel matrix is a **symbolic matrix** defined by a mathematical formula and no approximation is required in the computation of the interactions (up to the time discretization). The SiSyPHE library speeds up both traditional Python and low-level CPU implementations by **one to three orders of magnitude** for systems with up to several millions of particles. Although the library is mainly intended to be used on a GPU, the implementation is fully functional on the CPU with a significant gain in efficiency.   

Moreover, the **versatile object-oriented Python interface** is well suited to the comparison and study of new and classical many-particle models. This aspect is fundamental in applications in order to conduct ambitious numerical experiments in a systematic framework, even for particles with a complex structure and with a significantly reduced computational cost

References
===============

.. bibliography::



.. Sisyphe documentation master file, created by
   sphinx-quickstart on Thu Apr 29 16:21:29 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Simulation of Systems of interacting mean-field Particles with High Efficiency
************************************************************************************

.. image:: _static/ball.png
    :scale: 70% 
    :alt: Ball
    :align: center

The SiSyPHE library builds on recent advances in hardware and software for the efficient simulation of **large scale interacting particle systems**, both on the **GPU** and on the CPU. The implementation is based on recent libraries originally developed for machine learning purposes to significantly accelerate tensor (array) computations, namely the `PyTorch <https://github.com/pytorch/pytorch>`_ package and the `KeOps <https://www.kernel-operations.io/keops/index.html>`_ library. The **versatile object-oriented Python interface** is well suited to the comparison of new and classical many-particle models, enabling ambitious numerical experiments and leading to novel conjectures. The SiSyPHE library speeds up both traditional Python and low-level implementations by **one to three orders of magnitude** for systems with up to **several millions** of particles. 

The project is hosted on `GitHub <https://github.com/antoinediez/Sisyphe>`_, under the permissive `MIT license <https://en.wikipedia.org/wiki/MIT_License>`_. |br|  
|PyPi version| |JOSS paper| |doc badge|


.. image:: _static/mill.gif
    :scale: 100% 
    :alt: Mill
    :align: center
    
Citation
=====================

If you use SiSyPHE in a research paper, please cite the `JOSS publication <https://joss.theoj.org/papers/10.21105/joss.03653#>`_ : ::

  @article{Diez2021,
  doi = {10.21105/joss.03653},
  url = {https://doi.org/10.21105/joss.03653},
  year = {2021},
  publisher = {The Open Journal},
  volume = {6},
  number = {65},
  pages = {3653},
  author = {Antoine Diez},
  title = {`SiSyPHE`: A Python package for the Simulation of Systems of interacting mean-field Particles with High Efficiency},
  journal = {Journal of Open Source Software}}
  
Diez, A., (2021). SiSyPHE: A Python package for the Simulation of Systems of interacting mean-field Particles with High Efficiency. Journal of Open Source Software, 6(65), 3653, https://doi.org/10.21105/joss.03653

Contributing
======================

Contributions to make SiSyPHE grow are warmly welcome! Examples of possible (and ongoing) developments include the following. 

* The implementation of new models.

* The implementation of more complex boundary conditions and of models on *non-flat* manifolds. 

* An improved visualization method (currently only basic visualization functions relying on `Matplotlib <https://matplotlib.org/>`_ are implemented). 

Contributions can be made by opening an issue on the GitHub repository, via a pull request or by contacting directly the author.  

.. image:: _static/band.gif
    :scale: 100% 
    :alt: Band
    :align: center


Author
===========

`Antoine Diez <https://antoinediez.gitlab.io/>`_, Imperial College London 

Acknowledgments
-------------------

The development of this library would not have been possible without the help of `Jean Feydy  <https://www.jeanfeydy.com/>`_, his constant support and precious advice. This project was initiated by `Pierre Degond <https://sites.google.com/site/degond/>`_ and has grown out of many discussions with him. 

Table of contents
========================

.. toctree::
   :maxdepth: 2

   background.rst

.. toctree::
   :maxdepth: 2

   installation.rst

.. toctree::
   :maxdepth: 2

   benchmark.rst

.. toctree::
   :maxdepth: 2
   :caption: Tutorials and examples

   _auto_tutorials/index
   _auto_examples/index

.. toctree::
   :maxdepth: 2
   :caption: API

   api/API_particles
   api/API_models
   api/API_kernels
   api/API_sampling
   api/API_display
   api/API_toolbox
   

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. |PyPi version| image:: https://img.shields.io/pypi/v/sisyphe?color=blue
   :target: https://pypi.org/project/sisyphe/
   
.. |JOSS paper| image:: https://joss.theoj.org/papers/10.21105/joss.03653/status.svg
   :target: https://doi.org/10.21105/joss.03653
   
.. |doc badge| image:: https://readthedocs.org/projects/sisyphe/badge/?version=latest
   :target: https://sisyphe.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status
   
.. |br| raw:: html

  <br/>

:orphan:

.. _sphx_glr__auto_examples_sg_execution_times:

Computation times
=================
**29:34.952** total execution time for **_auto_examples** files:

+----------------------------------------------------------------------------------------+-----------+--------+
| :ref:`sphx_glr__auto_examples_plot_bands.py` (``plot_bands.py``)                       | 29:34.952 | 0.0 MB |
+----------------------------------------------------------------------------------------+-----------+--------+
| :ref:`sphx_glr__auto_examples_plot_bo_milling.py` (``plot_bo_milling.py``)             | 00:00.000 | 0.0 MB |
+----------------------------------------------------------------------------------------+-----------+--------+
| :ref:`sphx_glr__auto_examples_plot_milling.py` (``plot_milling.py``)                   | 00:00.000 | 0.0 MB |
+----------------------------------------------------------------------------------------+-----------+--------+
| :ref:`sphx_glr__auto_examples_plot_sphere_bands.py` (``plot_sphere_bands.py``)         | 00:00.000 | 0.0 MB |
+----------------------------------------------------------------------------------------+-----------+--------+
| :ref:`sphx_glr__auto_examples_plot_volume_exclusion.py` (``plot_volume_exclusion.py``) | 00:00.000 | 0.0 MB |
+----------------------------------------------------------------------------------------+-----------+--------+

.. DO NOT EDIT.
.. THIS FILE WAS AUTOMATICALLY GENERATED BY SPHINX-GALLERY.
.. TO MAKE CHANGES, EDIT THE SOURCE PYTHON FILE:
.. "_auto_examples/plot_milling.py"
.. LINE NUMBERS ARE GIVEN BELOW.

.. only:: html

    .. note::
        :class: sphx-glr-download-link-note

        Click :ref:`here <sphx_glr_download__auto_examples_plot_milling.py>`
        to download the full example code

.. rst-class:: sphx-glr-example-title

.. _sphx_glr__auto_examples_plot_milling.py:


.. _examplemill:

Mills
============================================

Examples of milling behaviours. 

.. GENERATED FROM PYTHON SOURCE LINES 12-13

First of all, some standard imports. 

.. GENERATED FROM PYTHON SOURCE LINES 13-28

.. code-block:: default


    import os
    import sys
    import time
    import math
    import torch
    import numpy as np 
    from matplotlib import pyplot as plt
    import sisyphe.models as models
    from sisyphe.display import display_kinetic_particles

    use_cuda = torch.cuda.is_available()
    dtype = torch.cuda.FloatTensor if use_cuda else torch.FloatTensor









.. GENERATED FROM PYTHON SOURCE LINES 29-53

Milling in the D'Orsogna et al. model
-----------------------------------------

Let us create an instance of the attraction-repulsion model introduced in

M. R. D’Orsogna, Y. L. Chuang, A. L. Bertozzi, L. S. Chayes, Self-Propelled Particles with Soft-Core Interactions: Patterns, Stability, and Collapse, *Phys. Rev. Lett.*, Vol. 96, No. 10 (2006).

The particle system satisfies the ODE:

.. math::  

        \frac{\mathrm{d}X^i_t}{\mathrm{d}t} = V^i_t

.. math::

         \frac{\mathrm{d}V^i_t}{\mathrm{d}t} = (\alpha-\beta|V^i_t|^2)V^i_t - \frac{m}{N}\nabla_{x^i}\sum_{j\ne i} U(|X^i_t-X^j_t|)

where :math:`U` is the Morse potential 

.. math::

       U(r) := -C_a\mathrm{e}^{-r/\ell_a}+C_r\mathrm{e}^{-r/\ell_r}



.. GENERATED FROM PYTHON SOURCE LINES 53-89

.. code-block:: default



    N = 10000
    mass = 1000.
    L = 10. 

    Ca = .5
    la = 2.
    Cr = 1.
    lr = .5

    alpha = 1.6
    beta = .5
    v0 = math.sqrt(alpha/beta)

    pos = L*torch.rand((N,2)).type(dtype)
    vel = torch.randn(N,2).type(dtype)
    vel = vel/torch.norm(vel,dim=1).reshape((N,1))
    vel = v0*vel

    dt = .01

    simu = models.AttractionRepulsion(pos=pos,
                     vel=vel,
                     interaction_radius=math.sqrt(mass),
                     box_size=L,
                     propulsion = alpha,
                     friction = beta,                        
                     Ca = Ca,
                     la = la,
                     Cr = Cr,
                     lr = lr,                        
                     dt=dt,
                     p=1,                        
                     isaverage=True)








.. GENERATED FROM PYTHON SOURCE LINES 90-92

Run the simulation over 100 units of time and plot 10 frames. The ODE system is solved using the Runge-Kutta 4 numerical scheme. 


.. GENERATED FROM PYTHON SOURCE LINES 92-99

.. code-block:: default


    frames = [0,1,2,3,4,5,10,40,70,100]

    s = time.time()
    it, op = display_kinetic_particles(simu,frames)
    e = time.time()




.. rst-class:: sphx-glr-horizontal


    *

      .. image:: /_auto_examples/images/sphx_glr_plot_milling_001.png
          :alt: Self-propulsion and Attraction-Repulsion  Parameters: N=10000 ; alpha=1.6 ; beta=0.5 ; Ca=0.5 ; la=2.0 ; Cr=1.0 ; lr=0.5  Time=0.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_milling_002.png
          :alt: Self-propulsion and Attraction-Repulsion  Parameters: N=10000 ; alpha=1.6 ; beta=0.5 ; Ca=0.5 ; la=2.0 ; Cr=1.0 ; lr=0.5  Time=1.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_milling_003.png
          :alt: Self-propulsion and Attraction-Repulsion  Parameters: N=10000 ; alpha=1.6 ; beta=0.5 ; Ca=0.5 ; la=2.0 ; Cr=1.0 ; lr=0.5  Time=2.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_milling_004.png
          :alt: Self-propulsion and Attraction-Repulsion  Parameters: N=10000 ; alpha=1.6 ; beta=0.5 ; Ca=0.5 ; la=2.0 ; Cr=1.0 ; lr=0.5  Time=3.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_milling_005.png
          :alt: Self-propulsion and Attraction-Repulsion  Parameters: N=10000 ; alpha=1.6 ; beta=0.5 ; Ca=0.5 ; la=2.0 ; Cr=1.0 ; lr=0.5  Time=4.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_milling_006.png
          :alt: Self-propulsion and Attraction-Repulsion  Parameters: N=10000 ; alpha=1.6 ; beta=0.5 ; Ca=0.5 ; la=2.0 ; Cr=1.0 ; lr=0.5  Time=5.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_milling_007.png
          :alt: Self-propulsion and Attraction-Repulsion  Parameters: N=10000 ; alpha=1.6 ; beta=0.5 ; Ca=0.5 ; la=2.0 ; Cr=1.0 ; lr=0.5  Time=10.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_milling_008.png
          :alt: Self-propulsion and Attraction-Repulsion  Parameters: N=10000 ; alpha=1.6 ; beta=0.5 ; Ca=0.5 ; la=2.0 ; Cr=1.0 ; lr=0.5  Time=40.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_milling_009.png
          :alt: Self-propulsion and Attraction-Repulsion  Parameters: N=10000 ; alpha=1.6 ; beta=0.5 ; Ca=0.5 ; la=2.0 ; Cr=1.0 ; lr=0.5  Time=70.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_milling_010.png
          :alt: Self-propulsion and Attraction-Repulsion  Parameters: N=10000 ; alpha=1.6 ; beta=0.5 ; Ca=0.5 ; la=2.0 ; Cr=1.0 ; lr=0.5  Time=100.0
          :class: sphx-glr-multi-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none


    Progress:0%
    Progress:1%
    Progress:2%
    Progress:3%
    Progress:4%
    Progress:5%
    Progress:6%
    Progress:7%
    Progress:8%
    Progress:9%
    Progress:10%
    Progress:11%
    Progress:12%
    Progress:13%
    Progress:14%
    Progress:15%
    Progress:16%
    Progress:17%
    Progress:18%
    Progress:19%
    Progress:20%
    Progress:21%
    Progress:22%
    Progress:23%
    Progress:24%
    Progress:25%
    Progress:26%
    Progress:27%
    Progress:28%
    Progress:29%
    Progress:30%
    Progress:31%
    Progress:32%
    Progress:33%
    Progress:34%
    Progress:35%
    Progress:36%
    Progress:37%
    Progress:38%
    Progress:39%
    Progress:40%
    Progress:41%
    Progress:42%
    Progress:43%
    Progress:44%
    Progress:45%
    Progress:46%
    Progress:47%
    Progress:48%
    Progress:49%
    Progress:50%
    Progress:51%
    Progress:52%
    Progress:53%
    Progress:54%
    Progress:55%
    Progress:56%
    Progress:57%
    Progress:58%
    Progress:59%
    Progress:60%
    Progress:61%
    Progress:62%
    Progress:63%
    Progress:64%
    Progress:65%
    Progress:66%
    Progress:67%
    Progress:68%
    Progress:69%
    Progress:70%
    Progress:71%
    Progress:72%
    Progress:73%
    Progress:74%
    Progress:75%
    Progress:76%
    Progress:77%
    Progress:78%
    Progress:79%
    Progress:80%
    Progress:81%
    Progress:82%
    Progress:83%
    Progress:84%
    Progress:85%
    Progress:86%
    Progress:87%
    Progress:88%
    Progress:89%
    Progress:90%
    Progress:91%
    Progress:92%
    Progress:93%
    Progress:94%
    Progress:95%
    Progress:96%
    Progress:97%
    Progress:98%
    Progress:99%



.. GENERATED FROM PYTHON SOURCE LINES 100-101

Print the total simulation time and the average time per iteration. 

.. GENERATED FROM PYTHON SOURCE LINES 101-107

.. code-block:: default


    print('Total time: '+str(e-s)+' seconds')
    print('Average time per iteration: '+str((e-s)/simu.iteration)+' seconds')







.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Total time: 86.00608777999878 seconds
    Average time per iteration: 0.008600608777999877 seconds




.. GENERATED FROM PYTHON SOURCE LINES 108-116

Milling in the Vicsek model 
----------------------------------------

Let us create an instance of the Asynchronuos Vicsek model with a bounded cone of vision and a bounded angular velocity, as introduced in: 

A. Costanzo, C. K. Hemelrijk, Spontaneous emergence of milling (vortex state) in a Vicsek-like model, *J. Phys. D: Appl. Phys.*, 51, 134004



.. GENERATED FROM PYTHON SOURCE LINES 116-134

.. code-block:: default


    N = 10000
    R = 1.
    L = 20.

    nu = 1
    sigma = .02
    kappa = nu/sigma

    c = .175
    angvel_max = .175/nu

    pos = L*torch.rand((N,2)).type(dtype)
    vel = torch.randn(N,2).type(dtype)
    vel = vel/torch.norm(vel,dim=1).reshape((N,1))

    dt = .01








.. GENERATED FROM PYTHON SOURCE LINES 135-136

We add an option to the target

.. GENERATED FROM PYTHON SOURCE LINES 136-152

.. code-block:: default


    target = {"name" : "normalised", "parameters" : {}}
    option = {"bounded_angular_velocity" : {"angvel_max" : angvel_max, "dt" : 1./nu}}

    simu=models.AsynchronousVicsek(pos=pos,vel=vel,
                     v=c,
                     jump_rate=nu,kappa=kappa,
                     interaction_radius=R,
                     box_size=L,
                     vision_angle=math.pi, axis = None,
                     boundary_conditions='periodic',
                     variant=target,
                     options=option,
                     sampling_method='projected_normal')









.. GENERATED FROM PYTHON SOURCE LINES 153-155

Run the simulation over 200 units of time and plot 10 frames. 


.. GENERATED FROM PYTHON SOURCE LINES 155-164

.. code-block:: default


    # sphinx_gallery_thumbnail_number = -1

    frames = [0, 10, 30, 50, 75, 100, 125, 150, 175, 200]

    s = time.time()
    it, op = display_kinetic_particles(simu,frames)
    e = time.time()




.. rst-class:: sphx-glr-horizontal


    *

      .. image:: /_auto_examples/images/sphx_glr_plot_milling_011.png
          :alt: Asynchronous Vicsek (normalised)  Parameters: N=10000 ; R=1.0 ; Jump Rate=1 ; kappa=50.0 ; v=0.17  Time=0.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_milling_012.png
          :alt: Asynchronous Vicsek (normalised)  Parameters: N=10000 ; R=1.0 ; Jump Rate=1 ; kappa=50.0 ; v=0.17  Time=10.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_milling_013.png
          :alt: Asynchronous Vicsek (normalised)  Parameters: N=10000 ; R=1.0 ; Jump Rate=1 ; kappa=50.0 ; v=0.17  Time=30.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_milling_014.png
          :alt: Asynchronous Vicsek (normalised)  Parameters: N=10000 ; R=1.0 ; Jump Rate=1 ; kappa=50.0 ; v=0.17  Time=50.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_milling_015.png
          :alt: Asynchronous Vicsek (normalised)  Parameters: N=10000 ; R=1.0 ; Jump Rate=1 ; kappa=50.0 ; v=0.17  Time=75.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_milling_016.png
          :alt: Asynchronous Vicsek (normalised)  Parameters: N=10000 ; R=1.0 ; Jump Rate=1 ; kappa=50.0 ; v=0.17  Time=100.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_milling_017.png
          :alt: Asynchronous Vicsek (normalised)  Parameters: N=10000 ; R=1.0 ; Jump Rate=1 ; kappa=50.0 ; v=0.17  Time=125.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_milling_018.png
          :alt: Asynchronous Vicsek (normalised)  Parameters: N=10000 ; R=1.0 ; Jump Rate=1 ; kappa=50.0 ; v=0.17  Time=150.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_milling_019.png
          :alt: Asynchronous Vicsek (normalised)  Parameters: N=10000 ; R=1.0 ; Jump Rate=1 ; kappa=50.0 ; v=0.17  Time=175.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_milling_020.png
          :alt: Asynchronous Vicsek (normalised)  Parameters: N=10000 ; R=1.0 ; Jump Rate=1 ; kappa=50.0 ; v=0.17  Time=200.0
          :class: sphx-glr-multi-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none


    Progress:0%
    Progress:1%
    Progress:2%
    Progress:3%
    Progress:4%
    Progress:5%
    Progress:6%
    Progress:7%
    Progress:8%
    Progress:9%
    Progress:10%
    Progress:11%
    Progress:12%
    Progress:13%
    Progress:14%
    Progress:15%
    Progress:16%
    Progress:17%
    Progress:18%
    Progress:19%
    Progress:20%
    Progress:21%
    Progress:22%
    Progress:23%
    Progress:24%
    Progress:25%
    Progress:26%
    Progress:27%
    Progress:28%
    Progress:29%
    Progress:30%
    Progress:31%
    Progress:32%
    Progress:33%
    Progress:34%
    Progress:35%
    Progress:36%
    Progress:37%
    Progress:38%
    Progress:39%
    Progress:40%
    Progress:41%
    Progress:42%
    Progress:43%
    Progress:44%
    Progress:45%
    Progress:46%
    Progress:47%
    Progress:48%
    Progress:49%
    Progress:50%
    Progress:51%
    Progress:52%
    Progress:53%
    Progress:54%
    Progress:55%
    Progress:56%
    Progress:57%
    Progress:58%
    Progress:59%
    Progress:60%
    Progress:61%
    Progress:62%
    Progress:63%
    Progress:64%
    Progress:65%
    Progress:66%
    Progress:67%
    Progress:68%
    Progress:69%
    Progress:70%
    Progress:71%
    Progress:72%
    Progress:73%
    Progress:74%
    Progress:75%
    Progress:76%
    Progress:77%
    Progress:78%
    Progress:79%
    Progress:80%
    Progress:81%
    Progress:82%
    Progress:83%
    Progress:84%
    Progress:85%
    Progress:86%
    Progress:87%
    Progress:88%
    Progress:89%
    Progress:90%
    Progress:91%
    Progress:92%
    Progress:93%
    Progress:94%
    Progress:95%
    Progress:96%
    Progress:97%
    Progress:98%
    Progress:99%



.. GENERATED FROM PYTHON SOURCE LINES 165-166

Print the total simulation time and the average time per iteration. 

.. GENERATED FROM PYTHON SOURCE LINES 166-172

.. code-block:: default


    print('Total time: '+str(e-s)+' seconds')
    print('Average time per iteration: '+str((e-s)/simu.iteration)+' seconds')







.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Total time: 69.85851073265076 seconds
    Average time per iteration: 0.003492925536632538 seconds





.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 2 minutes  38.151 seconds)


.. _sphx_glr_download__auto_examples_plot_milling.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download sphx-glr-download-python

     :download:`Download Python source code: plot_milling.py <plot_milling.py>`



  .. container:: sphx-glr-download sphx-glr-download-jupyter

     :download:`Download Jupyter notebook: plot_milling.ipynb <plot_milling.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_

.. DO NOT EDIT.
.. THIS FILE WAS AUTOMATICALLY GENERATED BY SPHINX-GALLERY.
.. TO MAKE CHANGES, EDIT THE SOURCE PYTHON FILE:
.. "_auto_examples/plot_bands.py"
.. LINE NUMBERS ARE GIVEN BELOW.

.. only:: html

    .. note::
        :class: sphx-glr-download-link-note

        Click :ref:`here <sphx_glr_download__auto_examples_plot_bands.py>`
        to download the full example code

.. rst-class:: sphx-glr-example-title

.. _sphx_glr__auto_examples_plot_bands.py:


Bands
============================================

The classical Vicsek model in a square periodic domain `is known <https://arxiv.org/abs/0712.2062>`_ to produce band-like structures in a very dilute regime. These structures also appears in a mean-field regime. To showcase the efficiency of the SiSyPHE library, we simulate a mean-field :class:`Vicsek <sisyphe.models.Vicsek>` model with the target :meth:`max_kappa() <sisyphe.particles.KineticParticles.max_kappa>` and :math:`10^6` particles. 

.. GENERATED FROM PYTHON SOURCE LINES 10-11

First of all, some standard imports. 

.. GENERATED FROM PYTHON SOURCE LINES 11-24

.. code-block:: default


    import os
    import sys
    import time
    import math
    import torch
    import numpy as np 
    from matplotlib import pyplot as plt

    use_cuda = torch.cuda.is_available()
    dtype = torch.cuda.FloatTensor if use_cuda else torch.FloatTensor









.. GENERATED FROM PYTHON SOURCE LINES 25-26

Set the parameters and create an instance of the Vicsek model. 

.. GENERATED FROM PYTHON SOURCE LINES 26-57

.. code-block:: default


    import sisyphe.models as models

    N = 1000000
    L = 1.
    dt = .01

    nu = 3
    sigma = 1.
    kappa = nu/sigma

    R = .01
    c = .1

    pos = L*torch.rand((N,2)).type(dtype)
    vel = torch.randn(N,2).type(dtype)
    vel = vel/torch.norm(vel,dim=1).reshape((N,1))

    simu=models.Vicsek(pos=pos,vel=vel,
                 v=c,
                 sigma=sigma,nu=nu,
                 interaction_radius=R,
                 box_size=L,
                 boundary_conditions='periodic',
                 variant = {"name" : "max_kappa", "parameters" : {"kappa_max" : 10.}},
                 options = {},
                 numerical_scheme='projection',
                 dt=dt,
                 block_sparse_reduction=True)









.. GENERATED FROM PYTHON SOURCE LINES 58-59

Check that we are in a mean field regime... 

.. GENERATED FROM PYTHON SOURCE LINES 59-66

.. code-block:: default


    Nneigh = simu.number_of_neighbours()

    print("The most isolated particle has " + str(Nneigh.min().item()) + " neighbours.")
    print("The least isolated particle has " + str(Nneigh.max().item()) + " neighbours.")






.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    The most isolated particle has 232.0 neighbours.
    The least isolated particle has 404.0 neighbours.




.. GENERATED FROM PYTHON SOURCE LINES 67-68

Set the block sparse parameters to their optimal value. 

.. GENERATED FROM PYTHON SOURCE LINES 68-74

.. code-block:: default


    fastest, nb_cells, average_simu_time, simulation_time = simu.best_blocksparse_parameters(40,100)

    plt.plot(nb_cells,average_simu_time)
    plt.show()




.. image:: /_auto_examples/images/sphx_glr_plot_bands_001.png
    :alt: plot bands
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none


    Progress:0.0%
    Progress:1.67%
    Progress:3.33%
    Progress:5.0%
    Progress:6.67%
    Progress:8.33%
    Progress:10.0%
    Progress:11.67%
    Progress:13.33%
    Progress:15.0%
    Progress:16.67%
    Progress:18.33%
    Progress:20.0%
    Progress:21.67%
    Progress:23.33%
    Progress:25.0%
    Progress:26.67%
    Progress:28.33%
    Progress:30.0%
    Progress:31.67%
    Progress:33.33%
    Progress:35.0%
    Progress:36.67%
    Progress:38.33%
    Progress:40.0%
    Progress:41.67%
    Progress:43.33%
    Progress:45.0%
    Progress:46.67%
    Progress:48.33%
    Progress:50.0%
    Progress:51.67%
    Progress:53.33%
    Progress:55.0%
    Progress:56.67%
    Progress:58.33%
    Progress:60.0%
    Progress:61.67%
    Progress:63.33%
    Progress:65.0%
    Progress:66.67%
    Progress:68.33%
    Progress:70.0%
    Progress:71.67%
    Progress:73.33%
    Progress:75.0%
    Progress:76.67%
    Progress:78.33%
    Progress:80.0%
    Progress:81.67%
    Progress:83.33%
    Progress:85.0%
    Progress:86.67%
    Progress:88.33%
    Progress:90.0%
    Progress:91.67%
    Progress:93.33%
    Progress:95.0%
    Progress:96.67%
    Progress:98.33%



.. GENERATED FROM PYTHON SOURCE LINES 75-76

Create the function which compute the center of mass of the system (on the torus).

.. GENERATED FROM PYTHON SOURCE LINES 76-87

.. code-block:: default


    def center_of_mass(particles):
        cos_pos = torch.cos((2*math.pi / L) * particles.pos)
        sin_pos = torch.sin((2*math.pi / L) * particles.pos)
        average_cos = cos_pos.sum(0)
        average_sin = sin_pos.sum(0)
        center = torch.atan2(average_sin, average_cos)
        center = (L / (2*math.pi)) * torch.remainder(center, 2*math.pi)
        return center









.. GENERATED FROM PYTHON SOURCE LINES 88-89

Let us save the positions and velocities of 100k particles and the center of mass of the system during 300 units of time. 

.. GENERATED FROM PYTHON SOURCE LINES 89-98

.. code-block:: default


    from sisyphe.display import save

    frames = [50., 100., 300.]

    s = time.time()
    data = save(simu,frames,["pos", "vel"],[center_of_mass], Nsaved=100000, save_file=False)
    e = time.time()





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none


    Progress:0%
    Progress:1%
    Progress:2%
    Progress:3%
    Progress:4%
    Progress:5%
    Progress:6%
    Progress:7%
    Progress:8%
    Progress:9%
    Progress:10%
    Progress:11%
    Progress:12%
    Progress:13%
    Progress:14%
    Progress:15%
    Progress:16%
    Progress:17%
    Progress:18%
    Progress:19%
    Progress:20%
    Progress:21%
    Progress:22%
    Progress:23%
    Progress:24%
    Progress:25%
    Progress:26%
    Progress:27%
    Progress:28%
    Progress:29%
    Progress:30%
    Progress:31%
    Progress:32%
    Progress:33%
    Progress:34%
    Progress:35%
    Progress:36%
    Progress:37%
    Progress:38%
    Progress:39%
    Progress:40%
    Progress:41%
    Progress:42%
    Progress:43%
    Progress:44%
    Progress:45%
    Progress:46%
    Progress:47%
    Progress:48%
    Progress:49%
    Progress:50%
    Progress:51%
    Progress:52%
    Progress:53%
    Progress:54%
    Progress:55%
    Progress:56%
    Progress:57%
    Progress:58%
    Progress:59%
    Progress:60%
    Progress:61%
    Progress:62%
    Progress:63%
    Progress:64%
    Progress:65%
    Progress:66%
    Progress:67%
    Progress:68%
    Progress:69%
    Progress:70%
    Progress:71%
    Progress:72%
    Progress:73%
    Progress:74%
    Progress:75%
    Progress:76%
    Progress:77%
    Progress:78%
    Progress:79%
    Progress:80%
    Progress:81%
    Progress:82%
    Progress:83%
    Progress:84%
    Progress:85%
    Progress:86%
    Progress:87%
    Progress:88%
    Progress:89%
    Progress:90%
    Progress:91%
    Progress:92%
    Progress:93%
    Progress:94%
    Progress:95%
    Progress:96%
    Progress:97%
    Progress:98%
    Progress:99%
    Progress:100%



.. GENERATED FROM PYTHON SOURCE LINES 99-100

Print the total simulation time and the average time per iteration. 

.. GENERATED FROM PYTHON SOURCE LINES 100-105

.. code-block:: default


    print('Total time: '+str(e-s)+' seconds')
    print('Average time per iteration: '+str((e-s)/simu.iteration)+' seconds')






.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Total time: 1619.642599105835 seconds
    Average time per iteration: 0.05398628709395804 seconds




.. GENERATED FROM PYTHON SOURCE LINES 106-107

At the end of the simulation, we plot the particles and the evolution of the center of mass. 

.. GENERATED FROM PYTHON SOURCE LINES 107-142

.. code-block:: default


    # sphinx_gallery_thumbnail_number = 2
    f = plt.figure(0, figsize=(12, 12))
    for frame in range(len(data["frames"])):
        x = data["pos"][frame][:,0].cpu()
        y = data["pos"][frame][:,1].cpu()
        u = data["vel"][frame][:,0].cpu()
        v = data["vel"][frame][:,1].cpu()
        ax = f.add_subplot(2,2,frame+1)
        plt.quiver(x,y,u,v)
        ax.set_xlim(xmin=0, xmax=simu.L[0].cpu())
        ax.set_ylim(ymin=0, ymax=simu.L[1].cpu())
        ax.set_title("time="+str(data["frames"][frame]))

    center = data["center_of_mass"]

    center_x = []
    center_y = []

    for c in center:
        center_x.append(c[0])
        center_y.append(c[1])

    f = plt.figure(1)
    plt.plot(data["time"],center_x)
    plt.ylabel("x-coordinate of the center of mass")
    plt.xlabel("time")


    f = plt.figure(2)
    plt.plot(data["time"],center_y)
    plt.ylabel("y-coordinate of the center of mass")
    plt.xlabel("time")
    plt.show()




.. rst-class:: sphx-glr-horizontal


    *

      .. image:: /_auto_examples/images/sphx_glr_plot_bands_002.png
          :alt: time=0, time=50.00999999999862, time=100.00000000001425, time=300.00999999987215
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_bands_003.png
          :alt: plot bands
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_bands_004.png
          :alt: plot bands
          :class: sphx-glr-multi-img





.. GENERATED FROM PYTHON SOURCE LINES 143-144

We are still in a mean-field regime. 

.. GENERATED FROM PYTHON SOURCE LINES 144-152

.. code-block:: default


    Nneigh = simu.number_of_neighbours()

    print("The most isolated particle has " + str(Nneigh.min().item()) + " neighbours.")
    print("The least isolated particle has " + str(Nneigh.max().item()) + " neighbours.")







.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    The most isolated particle has 34.0 neighbours.
    The least isolated particle has 6125.0 neighbours.





.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 29 minutes  34.952 seconds)


.. _sphx_glr_download__auto_examples_plot_bands.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download sphx-glr-download-python

     :download:`Download Python source code: plot_bands.py <plot_bands.py>`



  .. container:: sphx-glr-download sphx-glr-download-jupyter

     :download:`Download Jupyter notebook: plot_bands.ipynb <plot_bands.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_

.. DO NOT EDIT.
.. THIS FILE WAS AUTOMATICALLY GENERATED BY SPHINX-GALLERY.
.. TO MAKE CHANGES, EDIT THE SOURCE PYTHON FILE:
.. "_auto_examples/plot_volume_exclusion.py"
.. LINE NUMBERS ARE GIVEN BELOW.

.. only:: html

    .. note::
        :class: sphx-glr-download-link-note

        Click :ref:`here <sphx_glr_download__auto_examples_plot_volume_exclusion.py>`
        to download the full example code

.. rst-class:: sphx-glr-example-title

.. _sphx_glr__auto_examples_plot_volume_exclusion.py:


Volume exclusion
============================================

.. GENERATED FROM PYTHON SOURCE LINES 8-11

This model is introduced in

S. Motsch, D. Peurichard, From short-range repulsion to Hele-Shaw problem in a model of tumor growth, *J. Math. Biology*, Vol. 76, No. 1, 2017.

.. GENERATED FROM PYTHON SOURCE LINES 14-15

First of all, some standard imports. 

.. GENERATED FROM PYTHON SOURCE LINES 15-31

.. code-block:: default


    import os
    import sys
    import time
    import math
    import torch
    import numpy as np 
    from matplotlib import pyplot as plt
    import sisyphe.models as models
    from sisyphe.display import scatter_particles


    use_cuda = torch.cuda.is_available()
    dtype = torch.cuda.FloatTensor if use_cuda else torch.FloatTensor









.. GENERATED FROM PYTHON SOURCE LINES 32-49

Repulsion force
-------------------------------------------

Each particle is a disk with a (fixed) random radius. The particles repel each other when they overlap. The force exerted by a particle located at :math:`x_j` with radius :math:`R_j` on a particle located at :math:`x_i` with radius :math:`R_i` is

.. math::

    F = -\frac{\alpha}{R_i} \nabla_{x_i} U\left(\frac{|x_i - x_j|^2}{(R_i + R_j)^2}\right),

where the potential is

.. math::

    U(s) = -\log(s) + s - 1\,\,\text{for}\,\, s<1 \,\,\text{and}\,\, U(s) = 0\,\, \text{for}\,\, s>1.

Initially, the particles are clustered in a small region with a strong overlapping. 


.. GENERATED FROM PYTHON SOURCE LINES 49-69

.. code-block:: default



    N = 10000
    rmin = .1
    rmax = 1.
    R = (rmax-rmin)*torch.rand(N).type(dtype)+rmin
    L = 100.
    D0 = 20.
    pos = (D0*torch.rand((N,2)).type(dtype)-D0/2)+torch.tensor([L/2,L/2]).type(dtype)

    dt = .1

    simu = models.VolumeExclusion(pos=pos,
                     interaction_radius=R,
                     box_size=L,
                     alpha=2.5,
                     division_rate=0., 
                     death_rate=0.,                    
                     dt=dt)








.. GENERATED FROM PYTHON SOURCE LINES 70-72

Run the simulation over 200 units of time using an adaptive time-step which ensures that the energy :attr:`E` of the system decreases.


.. GENERATED FROM PYTHON SOURCE LINES 72-81

.. code-block:: default


    # sphinx_gallery_thumbnail_number = 13

    frames = [0,1,2,3,4,5,10,30,50,75,100,150,200]

    s = time.time()
    scatter_particles(simu,frames)
    e = time.time()




.. rst-class:: sphx-glr-horizontal


    *

      .. image:: /_auto_examples/images/sphx_glr_plot_volume_exclusion_001.png
          :alt: Volume exclusion  Parameters: N=10000 ; alpha=2.5 ; Division rate=0.0 ; Death rate=0.0  Time=0.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_volume_exclusion_002.png
          :alt: Volume exclusion  Parameters: N=10000 ; alpha=2.5 ; Division rate=0.0 ; Death rate=0.0  Time=1.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_volume_exclusion_003.png
          :alt: Volume exclusion  Parameters: N=10000 ; alpha=2.5 ; Division rate=0.0 ; Death rate=0.0  Time=2.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_volume_exclusion_004.png
          :alt: Volume exclusion  Parameters: N=10000 ; alpha=2.5 ; Division rate=0.0 ; Death rate=0.0  Time=3.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_volume_exclusion_005.png
          :alt: Volume exclusion  Parameters: N=10000 ; alpha=2.5 ; Division rate=0.0 ; Death rate=0.0  Time=4.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_volume_exclusion_006.png
          :alt: Volume exclusion  Parameters: N=10000 ; alpha=2.5 ; Division rate=0.0 ; Death rate=0.0  Time=5.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_volume_exclusion_007.png
          :alt: Volume exclusion  Parameters: N=10000 ; alpha=2.5 ; Division rate=0.0 ; Death rate=0.0  Time=10.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_volume_exclusion_008.png
          :alt: Volume exclusion  Parameters: N=10000 ; alpha=2.5 ; Division rate=0.0 ; Death rate=0.0  Time=30.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_volume_exclusion_009.png
          :alt: Volume exclusion  Parameters: N=10000 ; alpha=2.5 ; Division rate=0.0 ; Death rate=0.0  Time=50.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_volume_exclusion_010.png
          :alt: Volume exclusion  Parameters: N=10000 ; alpha=2.5 ; Division rate=0.0 ; Death rate=0.0  Time=75.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_volume_exclusion_011.png
          :alt: Volume exclusion  Parameters: N=10000 ; alpha=2.5 ; Division rate=0.0 ; Death rate=0.0  Time=100.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_volume_exclusion_012.png
          :alt: Volume exclusion  Parameters: N=10000 ; alpha=2.5 ; Division rate=0.0 ; Death rate=0.0  Time=150.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_volume_exclusion_013.png
          :alt: Volume exclusion  Parameters: N=10000 ; alpha=2.5 ; Division rate=0.0 ; Death rate=0.0  Time=200.0
          :class: sphx-glr-multi-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none


    Progress:0%
    Progress:1%
    Progress:2%
    Progress:3%
    Progress:4%
    Progress:5%
    Progress:6%
    Progress:7%
    Progress:8%
    Progress:9%
    Progress:10%
    Progress:11%
    Progress:12%
    Progress:13%
    Progress:14%
    Progress:15%
    Progress:16%
    Progress:17%
    Progress:18%
    Progress:19%
    Progress:20%
    Progress:21%
    Progress:22%
    Progress:23%
    Progress:24%
    Progress:25%
    Progress:26%
    Progress:27%
    Progress:28%
    Progress:29%
    Progress:30%
    Progress:31%
    Progress:32%
    Progress:33%
    Progress:34%
    Progress:35%
    Progress:36%
    Progress:37%
    Progress:38%
    Progress:39%
    Progress:40%
    Progress:41%
    Progress:42%
    Progress:43%
    Progress:44%
    Progress:45%
    Progress:46%
    Progress:47%
    Progress:48%
    Progress:49%
    Progress:50%
    Progress:51%
    Progress:52%
    Progress:53%
    Progress:54%
    Progress:55%
    Progress:56%
    Progress:57%
    Progress:58%
    Progress:59%
    Progress:60%
    Progress:61%
    Progress:62%
    Progress:63%
    Progress:64%
    Progress:65%
    Progress:66%
    Progress:67%
    Progress:68%
    Progress:69%
    Progress:70%
    Progress:71%
    Progress:72%
    Progress:73%
    Progress:74%
    Progress:75%
    Progress:76%
    Progress:77%
    Progress:78%
    Progress:79%
    Progress:80%
    Progress:81%
    Progress:82%
    Progress:83%
    Progress:84%
    Progress:85%
    Progress:86%
    Progress:87%
    Progress:88%
    Progress:89%
    Progress:90%
    Progress:91%
    Progress:92%
    Progress:93%
    Progress:94%
    Progress:95%
    Progress:96%
    Progress:97%
    Progress:98%
    Progress:99%
    Progress:100%



.. GENERATED FROM PYTHON SOURCE LINES 82-83

Print the total simulation time and the average time per iteration. 

.. GENERATED FROM PYTHON SOURCE LINES 83-87

.. code-block:: default


    print('Total time: '+str(e-s)+' seconds')
    print('Average time per iteration: '+str((e-s)/simu.iteration)+' seconds')





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Total time: 326.2268440723419 seconds
    Average time per iteration: 0.01703801347847401 seconds




.. GENERATED FROM PYTHON SOURCE LINES 88-89

Note this funny behaviour: the particles are clustered by size! 

.. GENERATED FROM PYTHON SOURCE LINES 91-96

Repulsion force, random births and random deaths
---------------------------------------------------------------------

Same system but this time, particles die at a constant rate and give birth to new particles at the same rate. A new particle is added next to its parent and has the same radius. 


.. GENERATED FROM PYTHON SOURCE LINES 96-116

.. code-block:: default


    N = 10000
    rmin = .1
    rmax = 1.
    R = (rmax-rmin)*torch.rand(N).type(dtype)+rmin
    L = 100.
    D0 = 20.
    pos = (D0*torch.rand((N,2)).type(dtype)-D0/2)+torch.tensor([L/2,L/2]).type(dtype)

    dt = .1

    simu = models.VolumeExclusion(pos=pos,
                     interaction_radius=R,
                     box_size=L,
                     alpha=2.5,
                     division_rate=.3, 
                     death_rate=.3,                    
                     dt=dt,
                     Nmax = 20000)








.. GENERATED FROM PYTHON SOURCE LINES 117-119

Run the simulation over 200 units of time using an adaptive time-step which ensures that the energy :attr:`E <sisyphe.models.VolumeExclusion.E>` of the system decreases.


.. GENERATED FROM PYTHON SOURCE LINES 119-126

.. code-block:: default


    frames = [0,1,2,3,4,5,10,30,50,75,100,150,200]

    s = time.time()
    scatter_particles(simu,frames)
    e = time.time()




.. rst-class:: sphx-glr-horizontal


    *

      .. image:: /_auto_examples/images/sphx_glr_plot_volume_exclusion_014.png
          :alt: Volume exclusion  Parameters: N=10000 ; alpha=2.5 ; Division rate=0.3 ; Death rate=0.3  Time=0.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_volume_exclusion_015.png
          :alt: Volume exclusion  Parameters: N=9985 ; alpha=2.5 ; Division rate=0.3 ; Death rate=0.3  Time=1.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_volume_exclusion_016.png
          :alt: Volume exclusion  Parameters: N=10007 ; alpha=2.5 ; Division rate=0.3 ; Death rate=0.3  Time=2.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_volume_exclusion_017.png
          :alt: Volume exclusion  Parameters: N=9950 ; alpha=2.5 ; Division rate=0.3 ; Death rate=0.3  Time=3.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_volume_exclusion_018.png
          :alt: Volume exclusion  Parameters: N=10014 ; alpha=2.5 ; Division rate=0.3 ; Death rate=0.3  Time=4.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_volume_exclusion_019.png
          :alt: Volume exclusion  Parameters: N=9972 ; alpha=2.5 ; Division rate=0.3 ; Death rate=0.3  Time=5.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_volume_exclusion_020.png
          :alt: Volume exclusion  Parameters: N=9874 ; alpha=2.5 ; Division rate=0.3 ; Death rate=0.3  Time=10.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_volume_exclusion_021.png
          :alt: Volume exclusion  Parameters: N=9540 ; alpha=2.5 ; Division rate=0.3 ; Death rate=0.3  Time=30.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_volume_exclusion_022.png
          :alt: Volume exclusion  Parameters: N=9245 ; alpha=2.5 ; Division rate=0.3 ; Death rate=0.3  Time=50.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_volume_exclusion_023.png
          :alt: Volume exclusion  Parameters: N=9550 ; alpha=2.5 ; Division rate=0.3 ; Death rate=0.3  Time=75.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_volume_exclusion_024.png
          :alt: Volume exclusion  Parameters: N=8759 ; alpha=2.5 ; Division rate=0.3 ; Death rate=0.3  Time=100.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_volume_exclusion_025.png
          :alt: Volume exclusion  Parameters: N=7871 ; alpha=2.5 ; Division rate=0.3 ; Death rate=0.3  Time=150.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_volume_exclusion_026.png
          :alt: Volume exclusion  Parameters: N=6484 ; alpha=2.5 ; Division rate=0.3 ; Death rate=0.3  Time=200.0
          :class: sphx-glr-multi-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none


    Progress:0%
    Progress:1%
    Progress:2%
    Progress:3%
    Progress:4%
    Progress:5%
    Progress:6%
    Progress:7%
    Progress:8%
    Progress:9%
    Progress:10%
    Progress:11%
    Progress:12%
    Progress:13%
    Progress:14%
    Progress:15%
    Progress:16%
    Progress:17%
    Progress:18%
    Progress:19%
    Progress:20%
    Progress:21%
    Progress:22%
    Progress:23%
    Progress:24%
    Progress:25%
    Progress:26%
    Progress:27%
    Progress:28%
    Progress:29%
    Progress:30%
    Progress:31%
    Progress:32%
    Progress:33%
    Progress:34%
    Progress:35%
    Progress:36%
    Progress:37%
    Progress:38%
    Progress:39%
    Progress:40%
    Progress:41%
    Progress:42%
    Progress:43%
    Progress:44%
    Progress:45%
    Progress:46%
    Progress:47%
    Progress:48%
    Progress:49%
    Progress:50%
    Progress:51%
    Progress:52%
    Progress:53%
    Progress:54%
    Progress:55%
    Progress:56%
    Progress:57%
    Progress:58%
    Progress:59%
    Progress:60%
    Progress:61%
    Progress:62%
    Progress:63%
    Progress:64%
    Progress:65%
    Progress:66%
    Progress:67%
    Progress:68%
    Progress:69%
    Progress:70%
    Progress:71%
    Progress:72%
    Progress:73%
    Progress:74%
    Progress:75%
    Progress:76%
    Progress:77%
    Progress:78%
    Progress:79%
    Progress:80%
    Progress:81%
    Progress:82%
    Progress:83%
    Progress:84%
    Progress:85%
    Progress:86%
    Progress:87%
    Progress:88%
    Progress:89%
    Progress:90%
    Progress:91%
    Progress:92%
    Progress:93%
    Progress:94%
    Progress:95%
    Progress:96%
    Progress:97%
    Progress:98%
    Progress:99%
    Progress:100%



.. GENERATED FROM PYTHON SOURCE LINES 127-128

Print the total simulation time and the average time per iteration. 

.. GENERATED FROM PYTHON SOURCE LINES 128-132

.. code-block:: default


    print('Total time: '+str(e-s)+' seconds')
    print('Average time per iteration: '+str((e-s)/simu.iteration)+' seconds')





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Total time: 200.29545760154724 seconds
    Average time per iteration: 0.007975450250917705 seconds





.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 8 minutes  53.860 seconds)


.. _sphx_glr_download__auto_examples_plot_volume_exclusion.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download sphx-glr-download-python

     :download:`Download Python source code: plot_volume_exclusion.py <plot_volume_exclusion.py>`



  .. container:: sphx-glr-download sphx-glr-download-jupyter

     :download:`Download Jupyter notebook: plot_volume_exclusion.ipynb <plot_volume_exclusion.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_

.. DO NOT EDIT.
.. THIS FILE WAS AUTOMATICALLY GENERATED BY SPHINX-GALLERY.
.. TO MAKE CHANGES, EDIT THE SOURCE PYTHON FILE:
.. "_auto_examples/plot_sphere_bands.py"
.. LINE NUMBERS ARE GIVEN BELOW.

.. only:: html

    .. note::
        :class: sphx-glr-download-link-note

        Click :ref:`here <sphx_glr_download__auto_examples_plot_sphere_bands.py>`
        to download the full example code

.. rst-class:: sphx-glr-example-title

.. _sphx_glr__auto_examples_plot_sphere_bands.py:


Boundary clusters in a disk
============================================

A classical mean-field Vicsek model in a bounded disk domain.  

.. GENERATED FROM PYTHON SOURCE LINES 10-11

First of all, some standard imports. 

.. GENERATED FROM PYTHON SOURCE LINES 11-25

.. code-block:: default


    import os
    import sys
    import time
    import torch
    import numpy as np 
    from matplotlib import pyplot as plt
    from sisyphe.display import display_kinetic_particles


    use_cuda = torch.cuda.is_available()
    dtype = torch.cuda.FloatTensor if use_cuda else torch.FloatTensor









.. GENERATED FROM PYTHON SOURCE LINES 26-27

Set the parameters and create an instance of the Vicsek model. 

.. GENERATED FROM PYTHON SOURCE LINES 27-64

.. code-block:: default


    import sisyphe.models as models

    N = 1000000
    L = 10.
    dt = .01

    nu = 3
    sigma = 1.
    kappa = nu/sigma

    R = .1
    c = 1.

    center = torch.tensor([L/2,L/2]).type(dtype).reshape((1,2))
    radius = L/2
    pos = L*torch.rand((N,2)).type(dtype)
    out = ((pos-center)**2).sum(1) > radius**2
    while out.sum()>0:
        pos[out,:] = L*torch.rand((out.sum(),2)).type(dtype)
        out = ((pos-center)**2).sum(1) > radius**2
    vel = torch.randn(N,2).type(dtype)
    vel = vel/torch.norm(vel,dim=1).reshape((N,1))

    simu=models.Vicsek(pos=pos,vel=vel,
                 v=c,
                 sigma=sigma,nu=nu,
                 interaction_radius=R,
                 box_size=L,
                 boundary_conditions='spherical',
                 variant = {"name" : "max_kappa", "parameters" : {"kappa_max" : 10.}},
                 options = {},
                 numerical_scheme='projection',
                 dt=dt,
                 block_sparse_reduction=True)









.. GENERATED FROM PYTHON SOURCE LINES 65-66

Set the block sparse parameters to their optimal value. 

.. GENERATED FROM PYTHON SOURCE LINES 66-73

.. code-block:: default


    fastest, nb_cells, average_simu_time, simulation_time = simu.best_blocksparse_parameters(40,100)

    plt.plot(nb_cells,average_simu_time)
    plt.show()





.. image:: /_auto_examples/images/sphx_glr_plot_sphere_bands_001.png
    :alt: plot sphere bands
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none


    Progress:0.0%
    Progress:1.67%
    Progress:3.33%
    Progress:5.0%
    Progress:6.67%
    Progress:8.33%
    Progress:10.0%
    Progress:11.67%
    Progress:13.33%
    Progress:15.0%
    Progress:16.67%
    Progress:18.33%
    Progress:20.0%
    Progress:21.67%
    Progress:23.33%
    Progress:25.0%
    Progress:26.67%
    Progress:28.33%
    Progress:30.0%
    Progress:31.67%
    Progress:33.33%
    Progress:35.0%
    Progress:36.67%
    Progress:38.33%
    Progress:40.0%
    Progress:41.67%
    Progress:43.33%
    Progress:45.0%
    Progress:46.67%
    Progress:48.33%
    Progress:50.0%
    Progress:51.67%
    Progress:53.33%
    Progress:55.0%
    Progress:56.67%
    Progress:58.33%
    Progress:60.0%
    Progress:61.67%
    Progress:63.33%
    Progress:65.0%
    Progress:66.67%
    Progress:68.33%
    Progress:70.0%
    Progress:71.67%
    Progress:73.33%
    Progress:75.0%
    Progress:76.67%
    Progress:78.33%
    Progress:80.0%
    Progress:81.67%
    Progress:83.33%
    Progress:85.0%
    Progress:86.67%
    Progress:88.33%
    Progress:90.0%
    Progress:91.67%
    Progress:93.33%
    Progress:95.0%
    Progress:96.67%
    Progress:98.33%



.. GENERATED FROM PYTHON SOURCE LINES 74-75

Run the simulation and plot the particles. 

.. GENERATED FROM PYTHON SOURCE LINES 75-84

.. code-block:: default


    # sphinx_gallery_thumbnail_number = -1

    frames = [0, 10, 40, 70, 100, 150, 200, 250, 300]

    s = time.time()
    it, op = display_kinetic_particles(simu, frames, N_dispmax=100000)
    e = time.time()




.. rst-class:: sphx-glr-horizontal


    *

      .. image:: /_auto_examples/images/sphx_glr_plot_sphere_bands_002.png
          :alt: Vicsek(max_kappa)  Parameters: N=1000000 ; R=0.1 ; nu=3 ; sigma=1.0 ; v=1.0  Time=0.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_sphere_bands_003.png
          :alt: Vicsek(max_kappa)  Parameters: N=1000000 ; R=0.1 ; nu=3 ; sigma=1.0 ; v=1.0  Time=10.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_sphere_bands_004.png
          :alt: Vicsek(max_kappa)  Parameters: N=1000000 ; R=0.1 ; nu=3 ; sigma=1.0 ; v=1.0  Time=40.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_sphere_bands_005.png
          :alt: Vicsek(max_kappa)  Parameters: N=1000000 ; R=0.1 ; nu=3 ; sigma=1.0 ; v=1.0  Time=70.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_sphere_bands_006.png
          :alt: Vicsek(max_kappa)  Parameters: N=1000000 ; R=0.1 ; nu=3 ; sigma=1.0 ; v=1.0  Time=100.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_sphere_bands_007.png
          :alt: Vicsek(max_kappa)  Parameters: N=1000000 ; R=0.1 ; nu=3 ; sigma=1.0 ; v=1.0  Time=150.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_sphere_bands_008.png
          :alt: Vicsek(max_kappa)  Parameters: N=1000000 ; R=0.1 ; nu=3 ; sigma=1.0 ; v=1.0  Time=200.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_sphere_bands_009.png
          :alt: Vicsek(max_kappa)  Parameters: N=1000000 ; R=0.1 ; nu=3 ; sigma=1.0 ; v=1.0  Time=250.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_examples/images/sphx_glr_plot_sphere_bands_010.png
          :alt: Vicsek(max_kappa)  Parameters: N=1000000 ; R=0.1 ; nu=3 ; sigma=1.0 ; v=1.0  Time=300.0
          :class: sphx-glr-multi-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none


    Progress:0%
    Progress:1%
    Progress:2%
    Progress:3%
    Progress:4%
    Progress:5%
    Progress:6%
    Progress:7%
    Progress:8%
    Progress:9%
    Progress:10%
    Progress:11%
    Progress:12%
    Progress:13%
    Progress:14%
    Progress:15%
    Progress:16%
    Progress:17%
    Progress:18%
    Progress:19%
    Progress:20%
    Progress:21%
    Progress:22%
    Progress:23%
    Progress:24%
    Progress:25%
    Progress:26%
    Progress:27%
    Progress:28%
    Progress:29%
    Progress:30%
    Progress:31%
    Progress:32%
    Progress:33%
    Progress:34%
    Progress:35%
    Progress:36%
    Progress:37%
    Progress:38%
    Progress:39%
    Progress:40%
    Progress:41%
    Progress:42%
    Progress:43%
    Progress:44%
    Progress:45%
    Progress:46%
    Progress:47%
    Progress:48%
    Progress:49%
    Progress:50%
    Progress:51%
    Progress:52%
    Progress:53%
    Progress:54%
    Progress:55%
    Progress:56%
    Progress:57%
    Progress:58%
    Progress:59%
    Progress:60%
    Progress:61%
    Progress:62%
    Progress:63%
    Progress:64%
    Progress:65%
    Progress:66%
    Progress:67%
    Progress:68%
    Progress:69%
    Progress:70%
    Progress:71%
    Progress:72%
    Progress:73%
    Progress:74%
    Progress:75%
    Progress:76%
    Progress:77%
    Progress:78%
    Progress:79%
    Progress:80%
    Progress:81%
    Progress:82%
    Progress:83%
    Progress:84%
    Progress:85%
    Progress:86%
    Progress:87%
    Progress:88%
    Progress:89%
    Progress:90%
    Progress:91%
    Progress:92%
    Progress:93%
    Progress:94%
    Progress:95%
    Progress:96%
    Progress:97%
    Progress:98%
    Progress:99%



.. GENERATED FROM PYTHON SOURCE LINES 85-86

Print the total simulation time and the average time per iteration. 

.. GENERATED FROM PYTHON SOURCE LINES 86-92

.. code-block:: default


    print('Total time: '+str(e-s)+' seconds')
    print('Average time per iteration: '+str((e-s)/simu.iteration)+' seconds')







.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Total time: 1702.3468585014343 seconds
    Average time per iteration: 0.05674489528338114 seconds





.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 30 minutes  47.639 seconds)


.. _sphx_glr_download__auto_examples_plot_sphere_bands.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download sphx-glr-download-python

     :download:`Download Python source code: plot_sphere_bands.py <plot_sphere_bands.py>`



  .. container:: sphx-glr-download sphx-glr-download-jupyter

     :download:`Download Jupyter notebook: plot_sphere_bands.ipynb <plot_sphere_bands.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_

.. DO NOT EDIT.
.. THIS FILE WAS AUTOMATICALLY GENERATED BY SPHINX-GALLERY.
.. TO MAKE CHANGES, EDIT THE SOURCE PYTHON FILE:
.. "_auto_examples/plot_bo_milling.py"
.. LINE NUMBERS ARE GIVEN BELOW.

.. only:: html

    .. note::
        :class: sphx-glr-download-link-note

        Click :ref:`here <sphx_glr_download__auto_examples_plot_bo_milling.py>`
        to download the full example code

.. rst-class:: sphx-glr-example-title

.. _sphx_glr__auto_examples_plot_bo_milling.py:


Body-oriented mill
============================================

.. GENERATED FROM PYTHON SOURCE LINES 8-12

This model is introduced in 

P. Degond, A. Diez, M. Na, Bulk topological states in a new collective dynamics model,  arXiv:2101.10864, 2021 


.. GENERATED FROM PYTHON SOURCE LINES 14-15

First of all, some standard imports. 

.. GENERATED FROM PYTHON SOURCE LINES 15-30

.. code-block:: default


    import os
    import sys
    import time
    import math
    import torch
    import numpy as np 
    from matplotlib import pyplot as plt
    import sisyphe.models as models
    from sisyphe.display import save

    use_cuda = torch.cuda.is_available()
    dtype = torch.cuda.FloatTensor if use_cuda else torch.FloatTensor









.. GENERATED FROM PYTHON SOURCE LINES 31-41

Body-oriented particles with initial perpendicular twist
---------------------------------------------------------------

The system is composed of body-oriented particles which are initially uniformly scattered in a periodic box but their body-orientations are "twisted". The body orientation of a particle at position :math:`(x,y,z)` is initially: 

.. math:: 

    \left(\begin{array}{ccc} 1 & 0 & 0 \\ 0 & \cos(2\pi z) & -\sin(2\pi z) \\ 0 & \sin(2\pi z) & \cos(2\pi z)\end{array}\right).



.. GENERATED FROM PYTHON SOURCE LINES 41-66

.. code-block:: default



    from sisyphe.initial import cyclotron_twist_z

    N = 1500000
    L = 1
    R = .025
    nu = 40
    c = 1
    kappa = 10

    pos, bo = cyclotron_twist_z(N,L,1,kappa,dtype)

    simu = models.BOAsynchronousVicsek(pos=pos,bo=bo,
                     v=c,
                     jump_rate=nu,kappa=kappa,
                     interaction_radius=R,
                     box_size=L,
                     boundary_conditions='periodic',
                     variant = {"name" : "normalised", "parameters" : {}},
                     options = {},
                     sampling_method='vonmises',
                     block_sparse_reduction=True,
                     number_of_cells=15**3)








.. GENERATED FROM PYTHON SOURCE LINES 67-79

Run the simulation over 5 units of time and save the azimuthal angle of the mean direction of motion defined by: 

.. math::

    \varphi = \mathrm{arg}(\Omega^1+i\Omega^2) \in [0,2\pi],

where :math:`\Omega = (\Omega^1,\Omega^2,\Omega^3)` is the mean direction of motion of the particles with velocities :math:`(\Omega_i)_{1\leq i\leq N}` : 

.. math::

    \Omega := \frac{\sum_{i=1}^N \Omega_i}{|\sum_{i=1}^N \Omega_i|}


.. GENERATED FROM PYTHON SOURCE LINES 79-86

.. code-block:: default


    frames = [5.]

    s = time.time()
    data = save(simu,frames,[],["phi"],save_file=False)
    e = time.time()





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none


    Progress:0%
    Progress:1%
    Progress:2%
    Progress:3%
    Progress:4%
    Progress:5%
    Progress:6%
    Progress:7%
    Progress:8%
    Progress:9%
    Progress:10%
    Progress:11%
    Progress:12%
    Progress:13%
    Progress:14%
    Progress:15%
    Progress:16%
    Progress:17%
    Progress:18%
    Progress:19%
    Progress:20%
    Progress:21%
    Progress:22%
    Progress:23%
    Progress:24%
    Progress:25%
    Progress:26%
    Progress:27%
    Progress:28%
    Progress:29%
    Progress:30%
    Progress:31%
    Progress:32%
    Progress:33%
    Progress:34%
    Progress:35%
    Progress:36%
    Progress:37%
    Progress:38%
    Progress:39%
    Progress:40%
    Progress:41%
    Progress:42%
    Progress:43%
    Progress:44%
    Progress:45%
    Progress:46%
    Progress:47%
    Progress:48%
    Progress:49%
    Progress:50%
    Progress:51%
    Progress:52%
    Progress:53%
    Progress:54%
    Progress:55%
    Progress:56%
    Progress:57%
    Progress:58%
    Progress:59%
    Progress:60%
    Progress:61%
    Progress:62%
    Progress:63%
    Progress:64%
    Progress:65%
    Progress:66%
    Progress:67%
    Progress:68%
    Progress:69%
    Progress:70%
    Progress:71%
    Progress:72%
    Progress:73%
    Progress:74%
    Progress:75%
    Progress:76%
    Progress:77%
    Progress:78%
    Progress:79%
    Progress:80%
    Progress:81%
    Progress:82%
    Progress:83%
    Progress:84%
    Progress:85%
    Progress:86%
    Progress:87%
    Progress:88%
    Progress:89%
    Progress:90%
    Progress:91%
    Progress:92%
    Progress:93%
    Progress:94%
    Progress:95%
    Progress:96%
    Progress:97%
    Progress:98%
    Progress:99%
    Progress:100%



.. GENERATED FROM PYTHON SOURCE LINES 87-88

Print the total simulation time and the average time per iteration. 

.. GENERATED FROM PYTHON SOURCE LINES 88-94

.. code-block:: default



    print('Total time: '+str(e-s)+' seconds')
    print('Average time per iteration: '+str((e-s)/simu.iteration)+' seconds')






.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Total time: 1655.7933239936829 seconds
    Average time per iteration: 0.08278966619968414 seconds




.. GENERATED FROM PYTHON SOURCE LINES 95-96

Plot the azimuthal angle :math:`\varphi`. 

.. GENERATED FROM PYTHON SOURCE LINES 96-101

.. code-block:: default


    plt.plot(data["time"],data["phi"])
    plt.xlabel("time")
    plt.ylabel("Azimuthal angle")




.. image:: /_auto_examples/images/sphx_glr_plot_bo_milling_001.png
    :alt: plot bo milling
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none


    Text(55.847222222222214, 0.5, 'Azimuthal angle')




.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 27 minutes  35.923 seconds)


.. _sphx_glr_download__auto_examples_plot_bo_milling.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download sphx-glr-download-python

     :download:`Download Python source code: plot_bo_milling.py <plot_bo_milling.py>`



  .. container:: sphx-glr-download sphx-glr-download-jupyter

     :download:`Download Jupyter notebook: plot_bo_milling.ipynb <plot_bo_milling.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
:orphan:



.. _sphx_glr__auto_examples:

========
Examples
========

Some examples for various models. 


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="The classical Vicsek model in a square periodic domain `is known &lt;https://arxiv.org/abs/0712.20...">

.. only:: html

 .. figure:: /_auto_examples/images/thumb/sphx_glr_plot_bands_thumb.png
     :alt: Bands

     :ref:`sphx_glr__auto_examples_plot_bands.py`

.. raw:: html

    </div>


.. toctree::
   :hidden:

   /_auto_examples/plot_bands

.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Body-oriented mill">

.. only:: html

 .. figure:: /_auto_examples/images/thumb/sphx_glr_plot_bo_milling_thumb.png
     :alt: Body-oriented mill

     :ref:`sphx_glr__auto_examples_plot_bo_milling.py`

.. raw:: html

    </div>


.. toctree::
   :hidden:

   /_auto_examples/plot_bo_milling

.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Examples of milling behaviours. ">

.. only:: html

 .. figure:: /_auto_examples/images/thumb/sphx_glr_plot_milling_thumb.png
     :alt: Mills

     :ref:`sphx_glr__auto_examples_plot_milling.py`

.. raw:: html

    </div>


.. toctree::
   :hidden:

   /_auto_examples/plot_milling

.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="A classical mean-field Vicsek model in a bounded disk domain.  ">

.. only:: html

 .. figure:: /_auto_examples/images/thumb/sphx_glr_plot_sphere_bands_thumb.png
     :alt: Boundary clusters in a disk

     :ref:`sphx_glr__auto_examples_plot_sphere_bands.py`

.. raw:: html

    </div>


.. toctree::
   :hidden:

   /_auto_examples/plot_sphere_bands

.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Volume exclusion">

.. only:: html

 .. figure:: /_auto_examples/images/thumb/sphx_glr_plot_volume_exclusion_thumb.png
     :alt: Volume exclusion

     :ref:`sphx_glr__auto_examples_plot_volume_exclusion.py`

.. raw:: html

    </div>


.. toctree::
   :hidden:

   /_auto_examples/plot_volume_exclusion
.. raw:: html

    <div class="sphx-glr-clear"></div>



.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-gallery


  .. container:: sphx-glr-download sphx-glr-download-python

    :download:`Download all examples in Python source code: _auto_examples_python.zip </_auto_examples/_auto_examples_python.zip>`



  .. container:: sphx-glr-download sphx-glr-download-jupyter

    :download:`Download all examples in Jupyter notebooks: _auto_examples_jupyter.zip </_auto_examples/_auto_examples_jupyter.zip>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_

:orphan:

.. _sphx_glr__auto_tutorials_sg_execution_times:

Computation times
=================
**06:25.941** total execution time for **_auto_tutorials** files:

+-------------------------------------------------------------------------------------------------+-----------+--------+
| :ref:`sphx_glr__auto_tutorials_plot_b_targetsoptions.py` (``plot_b_targetsoptions.py``)         | 06:25.941 | 0.0 MB |
+-------------------------------------------------------------------------------------------------+-----------+--------+
| :ref:`sphx_glr__auto_tutorials_plot_a_particlesmodels.py` (``plot_a_particlesmodels.py``)       | 00:00.000 | 0.0 MB |
+-------------------------------------------------------------------------------------------------+-----------+--------+
| :ref:`sphx_glr__auto_tutorials_plot_c_boundaryconditions.py` (``plot_c_boundaryconditions.py``) | 00:00.000 | 0.0 MB |
+-------------------------------------------------------------------------------------------------+-----------+--------+
| :ref:`sphx_glr__auto_tutorials_plot_d_blocksparse.py` (``plot_d_blocksparse.py``)               | 00:00.000 | 0.0 MB |
+-------------------------------------------------------------------------------------------------+-----------+--------+
| :ref:`sphx_glr__auto_tutorials_plot_e_kernels.py` (``plot_e_kernels.py``)                       | 00:00.000 | 0.0 MB |
+-------------------------------------------------------------------------------------------------+-----------+--------+

.. DO NOT EDIT.
.. THIS FILE WAS AUTOMATICALLY GENERATED BY SPHINX-GALLERY.
.. TO MAKE CHANGES, EDIT THE SOURCE PYTHON FILE:
.. "_auto_tutorials/plot_c_boundaryconditions.py"
.. LINE NUMBERS ARE GIVEN BELOW.

.. only:: html

    .. note::
        :class: sphx-glr-download-link-note

        Click :ref:`here <sphx_glr_download__auto_tutorials_plot_c_boundaryconditions.py>`
        to download the full example code

.. rst-class:: sphx-glr-example-title

.. _sphx_glr__auto_tutorials_plot_c_boundaryconditions.py:


.. _tuto_boundaryconditions:

Tutorial 03: Boundary conditions
=======================================

The particle systems are defined in a rectangular box whose dimensions are specified by the attribute :attr:`L <sisyphe.Particles.particles.L>`. 

The boundary conditions are specified by the attribute :attr:`bc <sisyphe.Particles.particles.bc>` which can be one of the following. 

* A list of size :math:`d` containing for each dimension either 0 (periodic) or 1 (wall with reflecting boundary conditions).

* The string ``"open"`` : no boundary conditions.

* The string ``"periodic"`` : periodic boundary conditions.

* The string ``"spherical"`` : reflecting boundary conditions on the sphere of diameter :math:`L` enclosed in the square domain :math:`[0,L]^d`. 

.. GENERATED FROM PYTHON SOURCE LINES 22-26

For instance, let us simulate the Vicsek model in an elongated rectangular domain :math:`[0,L_x]\times[0,L_y]` with periodic boundary conditions in the :math:`x`-dimension and reflecting boundary conditions in the :math:`y`-dimension. 

First, some standard imports...


.. GENERATED FROM PYTHON SOURCE LINES 26-35

.. code-block:: default


    import time 
    import torch
    from sisyphe.models import Vicsek
    from sisyphe.display import display_kinetic_particles

    use_cuda = torch.cuda.is_available()
    dtype = torch.cuda.FloatTensor if use_cuda else torch.FloatTensor








.. GENERATED FROM PYTHON SOURCE LINES 36-37

The parameters of the model 

.. GENERATED FROM PYTHON SOURCE LINES 37-49

.. code-block:: default


    N = 100000

    R = .01
    c = .1
    nu = 3.
    sigma = 1.

    dt = .01

    variant = {"name" : "max_kappa", "parameters" : {"kappa_max" : 10.}}








.. GENERATED FROM PYTHON SOURCE LINES 50-52

The spatial domain, the boundary conditions and the initial conditions... 


.. GENERATED FROM PYTHON SOURCE LINES 52-79

.. code-block:: default


    Lx = 3.
    Ly = 1./3.
    L = [Lx, Ly]
    bc = [0,1]

    pos = torch.rand((N,2)).type(dtype)
    pos[:,0] = L[0]*pos[:,0]
    pos[:,1] = L[1]*pos[:,1]
    vel = torch.randn(N,2).type(dtype)
    vel = vel/torch.norm(vel,dim=1).reshape((N,1))

    simu = Vicsek(
        pos = pos.detach().clone(),
        vel = vel.detach().clone(), 
        v = c, 
        sigma = sigma, 
        nu = nu, 
        interaction_radius = R,
        box_size = L,
        boundary_conditions=bc,
        dt = dt,
        variant = variant,
        block_sparse_reduction = True,
        number_of_cells = 100**2)









.. GENERATED FROM PYTHON SOURCE LINES 80-81

Finally run the simulation over 300 units of time.... 

.. GENERATED FROM PYTHON SOURCE LINES 81-90

.. code-block:: default


    # sphinx_gallery_thumbnail_number = 15

    frames = [0., 2., 5., 10., 30., 42., 71., 100, 123, 141, 182, 203, 256, 272, 300]

    s = time.time()
    it, op = display_kinetic_particles(simu, frames, order=True, figsize=(8,3))
    e = time.time()




.. rst-class:: sphx-glr-horizontal


    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_c_boundaryconditions_001.png
          :alt: Vicsek(max_kappa)  Parameters: N=100000 ; R=0.01 ; nu=3.0 ; sigma=1.0 ; v=0.1  Time=0.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_c_boundaryconditions_002.png
          :alt: Vicsek(max_kappa)  Parameters: N=100000 ; R=0.01 ; nu=3.0 ; sigma=1.0 ; v=0.1  Time=2.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_c_boundaryconditions_003.png
          :alt: Vicsek(max_kappa)  Parameters: N=100000 ; R=0.01 ; nu=3.0 ; sigma=1.0 ; v=0.1  Time=5.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_c_boundaryconditions_004.png
          :alt: Vicsek(max_kappa)  Parameters: N=100000 ; R=0.01 ; nu=3.0 ; sigma=1.0 ; v=0.1  Time=10.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_c_boundaryconditions_005.png
          :alt: Vicsek(max_kappa)  Parameters: N=100000 ; R=0.01 ; nu=3.0 ; sigma=1.0 ; v=0.1  Time=30.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_c_boundaryconditions_006.png
          :alt: Vicsek(max_kappa)  Parameters: N=100000 ; R=0.01 ; nu=3.0 ; sigma=1.0 ; v=0.1  Time=42.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_c_boundaryconditions_007.png
          :alt: Vicsek(max_kappa)  Parameters: N=100000 ; R=0.01 ; nu=3.0 ; sigma=1.0 ; v=0.1  Time=71.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_c_boundaryconditions_008.png
          :alt: Vicsek(max_kappa)  Parameters: N=100000 ; R=0.01 ; nu=3.0 ; sigma=1.0 ; v=0.1  Time=100.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_c_boundaryconditions_009.png
          :alt: Vicsek(max_kappa)  Parameters: N=100000 ; R=0.01 ; nu=3.0 ; sigma=1.0 ; v=0.1  Time=123.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_c_boundaryconditions_010.png
          :alt: Vicsek(max_kappa)  Parameters: N=100000 ; R=0.01 ; nu=3.0 ; sigma=1.0 ; v=0.1  Time=141.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_c_boundaryconditions_011.png
          :alt: Vicsek(max_kappa)  Parameters: N=100000 ; R=0.01 ; nu=3.0 ; sigma=1.0 ; v=0.1  Time=182.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_c_boundaryconditions_012.png
          :alt: Vicsek(max_kappa)  Parameters: N=100000 ; R=0.01 ; nu=3.0 ; sigma=1.0 ; v=0.1  Time=203.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_c_boundaryconditions_013.png
          :alt: Vicsek(max_kappa)  Parameters: N=100000 ; R=0.01 ; nu=3.0 ; sigma=1.0 ; v=0.1  Time=256.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_c_boundaryconditions_014.png
          :alt: Vicsek(max_kappa)  Parameters: N=100000 ; R=0.01 ; nu=3.0 ; sigma=1.0 ; v=0.1  Time=272.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_c_boundaryconditions_015.png
          :alt: Vicsek(max_kappa)  Parameters: N=100000 ; R=0.01 ; nu=3.0 ; sigma=1.0 ; v=0.1  Time=300.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_c_boundaryconditions_016.png
          :alt: plot c boundaryconditions
          :class: sphx-glr-multi-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none


    Progress:0%
    Progress:1%
    Progress:2%
    Progress:3%
    Progress:4%
    Progress:5%
    Progress:6%
    Progress:7%
    Progress:8%
    Progress:9%
    Progress:10%
    Progress:11%
    Progress:12%
    Progress:13%
    Progress:14%
    Progress:15%
    Progress:16%
    Progress:17%
    Progress:18%
    Progress:19%
    Progress:20%
    Progress:21%
    Progress:22%
    Progress:23%
    Progress:24%
    Progress:25%
    Progress:26%
    Progress:27%
    Progress:28%
    Progress:29%
    Progress:30%
    Progress:31%
    Progress:32%
    Progress:33%
    Progress:34%
    Progress:35%
    Progress:36%
    Progress:37%
    Progress:38%
    Progress:39%
    Progress:40%
    Progress:41%
    Progress:42%
    Progress:43%
    Progress:44%
    Progress:45%
    Progress:46%
    Progress:47%
    Progress:48%
    Progress:49%
    Progress:50%
    Progress:51%
    Progress:52%
    Progress:53%
    Progress:54%
    Progress:55%
    Progress:56%
    Progress:57%
    Progress:58%
    Progress:59%
    Progress:60%
    Progress:61%
    Progress:62%
    Progress:63%
    Progress:64%
    Progress:65%
    Progress:66%
    Progress:67%
    Progress:68%
    Progress:69%
    Progress:70%
    Progress:71%
    Progress:72%
    Progress:73%
    Progress:74%
    Progress:75%
    Progress:76%
    Progress:77%
    Progress:78%
    Progress:79%
    Progress:80%
    Progress:81%
    Progress:82%
    Progress:83%
    Progress:84%
    Progress:85%
    Progress:86%
    Progress:87%
    Progress:88%
    Progress:89%
    Progress:90%
    Progress:91%
    Progress:92%
    Progress:93%
    Progress:94%
    Progress:95%
    Progress:96%
    Progress:97%
    Progress:98%
    Progress:99%



.. GENERATED FROM PYTHON SOURCE LINES 91-92

Print the total simulation time and the average time per iteration. 

.. GENERATED FROM PYTHON SOURCE LINES 92-96

.. code-block:: default


    print('Total time: '+str(e-s)+' seconds')
    print('Average time per iteration: '+str((e-s)/simu.iteration)+' seconds')





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Total time: 523.2031211853027 seconds
    Average time per iteration: 0.01744010403951009 seconds




.. GENERATED FROM PYTHON SOURCE LINES 97-98

The simulation produces small clusters moving from left to right or from right to left. Each "step" in the order parameter corresponds to a collision between two clusters moving in opposite directions. 


.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 8 minutes  55.368 seconds)


.. _sphx_glr_download__auto_tutorials_plot_c_boundaryconditions.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download sphx-glr-download-python

     :download:`Download Python source code: plot_c_boundaryconditions.py <plot_c_boundaryconditions.py>`



  .. container:: sphx-glr-download sphx-glr-download-jupyter

     :download:`Download Jupyter notebook: plot_c_boundaryconditions.ipynb <plot_c_boundaryconditions.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_

.. DO NOT EDIT.
.. THIS FILE WAS AUTOMATICALLY GENERATED BY SPHINX-GALLERY.
.. TO MAKE CHANGES, EDIT THE SOURCE PYTHON FILE:
.. "_auto_tutorials/plot_a_particlesmodels.py"
.. LINE NUMBERS ARE GIVEN BELOW.

.. only:: html

    .. note::
        :class: sphx-glr-download-link-note

        Click :ref:`here <sphx_glr_download__auto_tutorials_plot_a_particlesmodels.py>`
        to download the full example code

.. rst-class:: sphx-glr-example-title

.. _sphx_glr__auto_tutorials_plot_a_particlesmodels.py:


Tutorial 01: Particles and models 
============================================

.. GENERATED FROM PYTHON SOURCE LINES 8-23

A particle system is an instance of one of the classes defined in the module :mod:`sisyphe.particles`. 

Particles 
        The basic class :class:`sisyphe.particles.Particles` defines a particle system by the positions. 

Kinetic particles
        The class :class:`sisyphe.particles.KineticParticles` defines a particle system by the positions and the velocities.

Body-oriented particles.
        The class :class:`sisyphe.particles.BOParticles` defines a particle system in 3D by the positions and the body-orientations which are a rotation matrices in :math:`SO(3)` stored as quaternions. 

A model is a subclass of a particle class. Several examples are defined in the module :mod:`sisyphe.models`. For example, let us create an instance of the Vicsek model :class:`sisyphe.models.Vicsek` which is a subclass of :class:`sisyphe.particles.KineticParticles`. 

First, some standard imports...


.. GENERATED FROM PYTHON SOURCE LINES 23-27

.. code-block:: default


    import time 
    import torch








.. GENERATED FROM PYTHON SOURCE LINES 28-29

If CUDA is available, the computations will be done on the GPU and on the CPU otherwise. The type of the tensors (simple or double precision) are defined by the type of the initial conditions. Here and throughout the documentation, we work with single precision tensors. 

.. GENERATED FROM PYTHON SOURCE LINES 29-33

.. code-block:: default


    use_cuda = torch.cuda.is_available()
    dtype = torch.cuda.FloatTensor if use_cuda else torch.FloatTensor








.. GENERATED FROM PYTHON SOURCE LINES 34-35

We take initially :math:`N` particles uniformly scattered in a box of size :math:`L` with uniformly sampled directions of motion. 

.. GENERATED FROM PYTHON SOURCE LINES 35-43

.. code-block:: default


    N = 10000
    L = 100 

    pos = L*torch.rand((N,2)).type(dtype)
    vel = torch.randn(N,2).type(dtype)
    vel = vel/torch.norm(vel,dim=1).reshape((N,1))








.. GENERATED FROM PYTHON SOURCE LINES 44-45

Then we define the interaction radius :math:`R`, the speed of the particles :math:`c` and the drift and diffusion coefficients, respectively :math:`\nu` and :math:`\sigma`. 

.. GENERATED FROM PYTHON SOURCE LINES 45-51

.. code-block:: default


    R = 5.
    c = 1.
    nu = 3.
    sigma = 1.








.. GENERATED FROM PYTHON SOURCE LINES 52-53

We take a small discretisation time step. 

.. GENERATED FROM PYTHON SOURCE LINES 53-56

.. code-block:: default


    dt = .01








.. GENERATED FROM PYTHON SOURCE LINES 57-58

Finally, we define an instance of the Vicsek model with these parameters. 

.. GENERATED FROM PYTHON SOURCE LINES 58-71

.. code-block:: default


    from sisyphe.models import Vicsek

    simu = Vicsek(
        pos = pos,
        vel = vel, 
        v = c, 
        sigma = sigma, 
        nu = nu, 
        interaction_radius = R,
        box_size = L,
        dt = dt)








.. GENERATED FROM PYTHON SOURCE LINES 72-74

.. note::
        The boundary conditions are periodic by default, see :ref:`tuto_boundaryconditions`.

.. GENERATED FROM PYTHON SOURCE LINES 76-77

So far, nothing has been computed. All the particles are implemented as Python iterators: in order to compute the next time step of the algorithm, we can call the method :meth:`__next__`. This method increments the iteration counter by one and updates all the relevant quantities (positions and velocities) by calling the method :meth:`update() <sisyphe.models.Vicsek.update>` which defines the model. 

.. GENERATED FROM PYTHON SOURCE LINES 77-83

.. code-block:: default


    print("Current iteration: "+ str(simu.iteration))
    simu.__next__()
    print("Current iteration: "+ str(simu.iteration))






.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Current iteration: 0
    Current iteration: 1




.. GENERATED FROM PYTHON SOURCE LINES 84-85

On a longer time interval, we can use the methods in the module :mod:`sisyphe.display`. For instance, let us fix a list of time frames. 

.. GENERATED FROM PYTHON SOURCE LINES 85-88

.. code-block:: default


    frames = [5., 10., 30., 50., 75., 100]








.. GENERATED FROM PYTHON SOURCE LINES 89-90

Using the method :meth:`sisyphe.display.display_kinetic_particles`, the simulation will run until the last time in the list :data:`frames`. The method also displays a scatter plot of the particle system at each of the times specified in the list and finally compute and plot the order parameter. 

.. GENERATED FROM PYTHON SOURCE LINES 90-97

.. code-block:: default


    from sisyphe.display import display_kinetic_particles

    s = time.time()
    it, op = display_kinetic_particles(simu, frames, order=True)
    e = time.time()




.. rst-class:: sphx-glr-horizontal


    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_a_particlesmodels_001.png
          :alt: Vicsek(normalised)  Parameters: N=10000 ; R=5.0 ; nu=3.0 ; sigma=1.0 ; v=1.0  Time=5.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_a_particlesmodels_002.png
          :alt: Vicsek(normalised)  Parameters: N=10000 ; R=5.0 ; nu=3.0 ; sigma=1.0 ; v=1.0  Time=10.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_a_particlesmodels_003.png
          :alt: Vicsek(normalised)  Parameters: N=10000 ; R=5.0 ; nu=3.0 ; sigma=1.0 ; v=1.0  Time=30.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_a_particlesmodels_004.png
          :alt: Vicsek(normalised)  Parameters: N=10000 ; R=5.0 ; nu=3.0 ; sigma=1.0 ; v=1.0  Time=50.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_a_particlesmodels_005.png
          :alt: Vicsek(normalised)  Parameters: N=10000 ; R=5.0 ; nu=3.0 ; sigma=1.0 ; v=1.0  Time=75.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_a_particlesmodels_006.png
          :alt: Vicsek(normalised)  Parameters: N=10000 ; R=5.0 ; nu=3.0 ; sigma=1.0 ; v=1.0  Time=100.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_a_particlesmodels_007.png
          :alt: plot a particlesmodels
          :class: sphx-glr-multi-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none


    Progress:0%
    Progress:1%
    Progress:2%
    Progress:3%
    Progress:4%
    Progress:5%
    Progress:6%
    Progress:7%
    Progress:8%
    Progress:9%
    Progress:10%
    Progress:11%
    Progress:12%
    Progress:13%
    Progress:14%
    Progress:15%
    Progress:16%
    Progress:17%
    Progress:18%
    Progress:19%
    Progress:20%
    Progress:21%
    Progress:22%
    Progress:23%
    Progress:24%
    Progress:25%
    Progress:26%
    Progress:27%
    Progress:28%
    Progress:29%
    Progress:30%
    Progress:31%
    Progress:32%
    Progress:33%
    Progress:34%
    Progress:35%
    Progress:36%
    Progress:37%
    Progress:38%
    Progress:39%
    Progress:40%
    Progress:41%
    Progress:42%
    Progress:43%
    Progress:44%
    Progress:45%
    Progress:46%
    Progress:47%
    Progress:48%
    Progress:49%
    Progress:50%
    Progress:51%
    Progress:52%
    Progress:53%
    Progress:54%
    Progress:55%
    Progress:56%
    Progress:57%
    Progress:58%
    Progress:59%
    Progress:60%
    Progress:61%
    Progress:62%
    Progress:63%
    Progress:64%
    Progress:65%
    Progress:66%
    Progress:67%
    Progress:68%
    Progress:69%
    Progress:70%
    Progress:71%
    Progress:72%
    Progress:73%
    Progress:74%
    Progress:75%
    Progress:76%
    Progress:77%
    Progress:78%
    Progress:79%
    Progress:80%
    Progress:81%
    Progress:82%
    Progress:83%
    Progress:84%
    Progress:85%
    Progress:86%
    Progress:87%
    Progress:88%
    Progress:89%
    Progress:90%
    Progress:91%
    Progress:92%
    Progress:93%
    Progress:94%
    Progress:95%
    Progress:96%
    Progress:97%
    Progress:98%
    Progress:99%



.. GENERATED FROM PYTHON SOURCE LINES 98-99

Print the total simulation time and the average time per iteration. 

.. GENERATED FROM PYTHON SOURCE LINES 99-110

.. code-block:: default


    print('Total time: '+str(e-s)+' seconds')
    print('Average time per iteration: '+str((e-s)/simu.iteration)+' seconds')




    
    
    
    
    



.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Total time: 19.33805775642395 seconds
    Average time per iteration: 0.0019336124144009549 seconds





.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 0 minutes  22.882 seconds)


.. _sphx_glr_download__auto_tutorials_plot_a_particlesmodels.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download sphx-glr-download-python

     :download:`Download Python source code: plot_a_particlesmodels.py <plot_a_particlesmodels.py>`



  .. container:: sphx-glr-download sphx-glr-download-jupyter

     :download:`Download Jupyter notebook: plot_a_particlesmodels.ipynb <plot_a_particlesmodels.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_

.. DO NOT EDIT.
.. THIS FILE WAS AUTOMATICALLY GENERATED BY SPHINX-GALLERY.
.. TO MAKE CHANGES, EDIT THE SOURCE PYTHON FILE:
.. "_auto_tutorials/plot_d_blocksparse.py"
.. LINE NUMBERS ARE GIVEN BELOW.

.. only:: html

    .. note::
        :class: sphx-glr-download-link-note

        Click :ref:`here <sphx_glr_download__auto_tutorials_plot_d_blocksparse.py>`
        to download the full example code

.. rst-class:: sphx-glr-example-title

.. _sphx_glr__auto_tutorials_plot_d_blocksparse.py:


.. _tutobsr:

Tutorial 04: Block-sparse reduction 
============================================

In many cases, the interaction radius :math:`R` is much smaller than the size of the domain. Consequently, the sums in the local averages (see :ref:`tuto_averages`) contain only a small fraction of non zero terms. To gain in efficiency, we can follow the classical strategy:

* Subdivide the domain into a fixed number of cells of size at least :math:`R`.
* For a particle in a given cell, only look at the contiguous cells to compute the local averages. In dimension :math:`d`, there are :math:`3^d` contiguous cells (including the cell itself). 

A practical implementation is called the *Verlet list method*. However, the implementation below is different than the classical one. It is adapted from the `block-sparse reduction method <https://www.kernel-operations.io/keops/_auto_examples/pytorch/plot_grid_cluster_pytorch.html>`_ implemented in the `KeOps <https://www.kernel-operations.io/keops/index.html>`_ library. 

We illustrate the gain in efficency for the Vicsek model. 

.. note::
    The method is sub-optimal for moderate numbers of particles. As a rule of thumb, the block-sparse reduction method becomes useful for systems with at least :math:`10^4` particles. 

.. GENERATED FROM PYTHON SOURCE LINES 22-27

Set up and benchmarks
---------------------------

First, some standard imports...


.. GENERATED FROM PYTHON SOURCE LINES 27-36

.. code-block:: default


    import copy
    import time 
    import torch
    from matplotlib import pyplot as plt

    use_cuda = torch.cuda.is_available()
    dtype = torch.cuda.FloatTensor if use_cuda else torch.FloatTensor








.. GENERATED FROM PYTHON SOURCE LINES 37-39

Let the :math:`N` particles be uniformly scattered in a box of size :math:`L` with interaction radius  :math:`R` and uniformly sampled velocities. 


.. GENERATED FROM PYTHON SOURCE LINES 39-58

.. code-block:: default


    from sisyphe.models import Vicsek

    N = 100000
    L = 100.  
    R = 1.

    pos = L*torch.rand((N,2)).type(dtype)
    vel = torch.randn(N,2).type(dtype)
    vel = vel/torch.norm(vel,dim=1).reshape((N,1))

    simu=Vicsek(pos=pos,vel=vel,
                v=1.,
                sigma=1.,nu=3.,
                interaction_radius=R,
                box_size=L)

    simu.__next__() #GPU warmup... 





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none


    {'position': tensor([[31.6564,  6.3145],
            [61.3271, 90.9929],
            [88.1392,  1.8941],
            ...,
            [98.8685, 90.9011],
            [28.4884, 29.0810],
            [15.1606, 66.4232]], device='cuda:0'), 'velocity': tensor([[ 0.9581, -0.2866],
            [-0.2292, -0.9734],
            [-0.8279, -0.5608],
            ...,
            [-0.6027, -0.7979],
            [ 0.9594, -0.2821],
            [-0.4047, -0.9145]], device='cuda:0')}



.. GENERATED FROM PYTHON SOURCE LINES 59-60

Without block-sparse reduction, let us compute the simulation time of 100 iterations.

.. GENERATED FROM PYTHON SOURCE LINES 60-71

.. code-block:: default


    simu_copy = copy.deepcopy(simu) # Make a new deepcopy
    s = time.time()
    for k in range(100):
        simu_copy.__next__()
    e = time.time()

    simulation_time = e-s

    print("Average simulation time without block-sparse reduction: " + str(simulation_time) + " seconds.")





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Average simulation time without block-sparse reduction: 2.9741082191467285 seconds.




.. GENERATED FROM PYTHON SOURCE LINES 72-73

Then with block-sparse reduction... First, turn on the attribute :attr:`blocksparse <sispyphe.particles.Particles.blocksparse>`. 

.. GENERATED FROM PYTHON SOURCE LINES 73-76

.. code-block:: default

      
    simu.blocksparse = True
      







.. GENERATED FROM PYTHON SOURCE LINES 77-79

Then, we need to define the maximum number of cells. This can be set by the keyword argument ``number_of_cells`` when an instance of the class :class:`sisyphe.particles.Particles` is created. The number of cells has a strong influence on the efficiency of the method and should be chosen wisely.  When the optimal value is not known a priori, it is recommanded to use the  method :meth:`best_blocksparse_parameters() <sisyphe.particles.Particles.best_blocksparse_parameters>` which will time 100 iterations of the simulation for various numbers of cells and automatically choose the best one. Below, we test all the numbers of cells which are powers of the dimension (here :math:`d=2`) between :math:`10^2` and :math:`70^2`. 


.. GENERATED FROM PYTHON SOURCE LINES 79-85

.. code-block:: default

 
    ncell_min = 10
    ncell_max = 70
    fastest, nb_cells, average_simu_time, simulation_time = simu.best_blocksparse_parameters(ncell_min, ncell_max, step=1, nb_calls=100)






.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none


    Progress:0.0%
    Progress:1.67%
    Progress:3.33%
    Progress:5.0%
    Progress:6.67%
    Progress:8.33%
    Progress:10.0%
    Progress:11.67%
    Progress:13.33%
    Progress:15.0%
    Progress:16.67%
    Progress:18.33%
    Progress:20.0%
    Progress:21.67%
    Progress:23.33%
    Progress:25.0%
    Progress:26.67%
    Progress:28.33%
    Progress:30.0%
    Progress:31.67%
    Progress:33.33%
    Progress:35.0%
    Progress:36.67%
    Progress:38.33%
    Progress:40.0%
    Progress:41.67%
    Progress:43.33%
    Progress:45.0%
    Progress:46.67%
    Progress:48.33%
    Progress:50.0%
    Progress:51.67%
    Progress:53.33%
    Progress:55.0%
    Progress:56.67%
    Progress:58.33%
    Progress:60.0%
    Progress:61.67%
    Progress:63.33%
    Progress:65.0%
    Progress:66.67%
    Progress:68.33%
    Progress:70.0%
    Progress:71.67%
    Progress:73.33%
    Progress:75.0%
    Progress:76.67%
    Progress:78.33%
    Progress:80.0%
    Progress:81.67%
    Progress:83.33%
    Progress:85.0%
    Progress:86.67%
    Progress:88.33%
    Progress:90.0%
    Progress:91.67%
    Progress:93.33%
    Progress:95.0%
    Progress:96.67%
    Progress:98.33%



.. GENERATED FROM PYTHON SOURCE LINES 86-87

We plot the average simulation time as a function of the square root of the number of cells and print the best. 

.. GENERATED FROM PYTHON SOURCE LINES 87-94

.. code-block:: default


    plt.plot(nb_cells,average_simu_time)      
    plt.xlabel("Square root of the number of cells") 
    plt.ylabel("Simulation time") 

    print("Average simulation time with block-sparse reduction: " + str(average_simu_time.min()) + " seconds.")




.. image:: /_auto_tutorials/images/sphx_glr_plot_d_blocksparse_001.png
    :alt: plot d blocksparse
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Average simulation time with block-sparse reduction: 0.4012424945831299 seconds.




.. GENERATED FROM PYTHON SOURCE LINES 95-96

Same experiment with one million particles. 

.. GENERATED FROM PYTHON SOURCE LINES 96-123

.. code-block:: default


    N = 1000000
    L = 100.  
    R = 1.

    pos = L*torch.rand((N,2)).type(dtype)
    vel = torch.randn(N,2).type(dtype)
    vel = vel/torch.norm(vel,dim=1).reshape((N,1))


    simu=Vicsek(pos=pos,vel=vel,
                v=1.,
                sigma=1.,nu=3.,
                interaction_radius=R,
                box_size=L,
                block_sparse_reduction=False)

    simu_copy = copy.deepcopy(simu) # Make a new deepcopy
    s = time.time()
    for k in range(100):
        simu_copy.__next__()
    e = time.time()

    simulation_time = e-s

    print("Average simulation time without block-sparse reduction: " + str(simulation_time) + " seconds.")





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Average simulation time without block-sparse reduction: 274.7271876335144 seconds.




.. GENERATED FROM PYTHON SOURCE LINES 124-125

With block-sparse reduction...

.. GENERATED FROM PYTHON SOURCE LINES 125-131

.. code-block:: default


    simu.blocksparse = True

    fastest, nb_cells, average_simu_time, simulation_time = simu.best_blocksparse_parameters(30, 100, nb_calls=100)






.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none


    Progress:0.0%
    Progress:1.43%
    Progress:2.86%
    Progress:4.29%
    Progress:5.71%
    Progress:7.14%
    Progress:8.57%
    Progress:10.0%
    Progress:11.43%
    Progress:12.86%
    Progress:14.29%
    Progress:15.71%
    Progress:17.14%
    Progress:18.57%
    Progress:20.0%
    Progress:21.43%
    Progress:22.86%
    Progress:24.29%
    Progress:25.71%
    Progress:27.14%
    Progress:28.57%
    Progress:30.0%
    Progress:31.43%
    Progress:32.86%
    Progress:34.29%
    Progress:35.71%
    Progress:37.14%
    Progress:38.57%
    Progress:40.0%
    Progress:41.43%
    Progress:42.86%
    Progress:44.29%
    Progress:45.71%
    Progress:47.14%
    Progress:48.57%
    Progress:50.0%
    Progress:51.43%
    Progress:52.86%
    Progress:54.29%
    Progress:55.71%
    Progress:57.14%
    Progress:58.57%
    Progress:60.0%
    Progress:61.43%
    Progress:62.86%
    Progress:64.29%
    Progress:65.71%
    Progress:67.14%
    Progress:68.57%
    Progress:70.0%
    Progress:71.43%
    Progress:72.86%
    Progress:74.29%
    Progress:75.71%
    Progress:77.14%
    Progress:78.57%
    Progress:80.0%
    Progress:81.43%
    Progress:82.86%
    Progress:84.29%
    Progress:85.71%
    Progress:87.14%
    Progress:88.57%
    Progress:90.0%
    Progress:91.43%
    Progress:92.86%
    Progress:94.29%
    Progress:95.71%
    Progress:97.14%
    Progress:98.57%



.. GENERATED FROM PYTHON SOURCE LINES 132-133

We plot the average simulation time as a function of the square root of the number of cells and print the best. 

.. GENERATED FROM PYTHON SOURCE LINES 133-140

.. code-block:: default


    plt.plot(nb_cells,average_simu_time)      
    plt.xlabel("Square root of the number of cells") 
    plt.ylabel("Simulation time") 

    print("Average simulation time with block-sparse reduction: " + str(average_simu_time.min()) + " seconds.")




.. image:: /_auto_tutorials/images/sphx_glr_plot_d_blocksparse_002.png
    :alt: plot d blocksparse
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Average simulation time with block-sparse reduction: 2.2031195163726807 seconds.




.. GENERATED FROM PYTHON SOURCE LINES 141-143

.. note::
        The optimal parameters chosen initially may not stay optimal in the course of the simulation. This may be the case in particular if there is a strong concentration of particles.

.. GENERATED FROM PYTHON SOURCE LINES 145-178

How does it work 
---------------------------------

Cell size and number of cells
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The cells have a rectangular shape. The length of the cells along each dimension cannot be smaller than the interaction radius :math:`R`. The maximum number of cells is thus equal to: 

.. math::

        n_\mathrm{max} = \prod_{k=1}^d \left\lfloor \frac{L_k}{R} \right\rfloor,

where :math:`L_k` is the length of the (rectangular) domain along dimension :math:`k`. This corresponds to rectangular cells with a length along dimension :math:`k` equal to: 

.. math:: 

        \varepsilon_k = \frac{L_k}{\left\lfloor \frac{L_k}{R} \right\rfloor}.

If the number of cells demanded :math:`n_0` exceeds :math:`n_\mathrm{max}`, this will be the chosen value. Otherwise, we first compute the typical length: 

.. math:: 

        \varepsilon_0 = \left(\frac{\prod_{k=1}^d L_k}{n_0}\right)^{1/d}

Then the length of the cells along dimension :math:`k` is set to

.. math::

        \varepsilon_k = \frac{L_k}{\left\lfloor\frac{L_k}{\varepsilon_0}\right\rfloor}.

In particular, in a square domain :math:`L_k=L` for all :math:`k` and when :math:`n_0` is a power of :math:`d`, then there are exactly :math:`n_0` square cells with length :math:`L/n_0^{1/d}`. 



.. GENERATED FROM PYTHON SOURCE LINES 180-193

The block-sparse parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The initialisation or the method :meth:`best_blocksparse_parameters() <sisyphe.particles.Particles.best_blocksparse_parameters>` define three attributes which are used to speed up the computations. Given a number of cells, they are computed by the method :meth:`compute_blocksparse_parameters() <sisyphe.particles.Particles.compute_blocksparse_parameters>`.

* :attr:`centroids <sisyphe.Particles.particles.centroids>` : the coordinates of the centers of the cells. 
* :attr:`keep <sisyphe.Particles.particles.keep>` : a square BoolTensor which indicates whether two cells are contiguous. 
* :attr:`eps <sisyphe.Particles.particles.keep>` : the length of the cells along each dimension. 

The particles are clustered into the cells using the method :meth:`uniform_grid_separation() <sisyphe.toolbox.uniform_grid_separation>`. 

.. note::
    A drawback of the method is the high memory cost needed to store the boolean mask :attr:`keep <sisyphe.Particles.particles.keep>`. As a consequence, unlike the classical Verlet list method, the optimal number of cells is often **not** the maximum one. In the examples presented in this documentation, the optimal number of cells is always smaller than :math:`10^4`. 


.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 8 minutes  9.001 seconds)


.. _sphx_glr_download__auto_tutorials_plot_d_blocksparse.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download sphx-glr-download-python

     :download:`Download Python source code: plot_d_blocksparse.py <plot_d_blocksparse.py>`



  .. container:: sphx-glr-download sphx-glr-download-jupyter

     :download:`Download Jupyter notebook: plot_d_blocksparse.ipynb <plot_d_blocksparse.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
:orphan:



.. _sphx_glr__auto_tutorials:

=========
Tutorials
=========

Basic features. 


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Tutorial 01: Particles and models">

.. only:: html

 .. figure:: /_auto_tutorials/images/thumb/sphx_glr_plot_a_particlesmodels_thumb.png
     :alt: Tutorial 01: Particles and models

     :ref:`sphx_glr__auto_tutorials_plot_a_particlesmodels.py`

.. raw:: html

    </div>


.. toctree::
   :hidden:

   /_auto_tutorials/plot_a_particlesmodels

.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Many useful **local averages** (see tuto_averages) are pre-defined and may be directly used in ...">

.. only:: html

 .. figure:: /_auto_tutorials/images/thumb/sphx_glr_plot_b_targetsoptions_thumb.png
     :alt: Tutorial 02: Targets and options

     :ref:`sphx_glr__auto_tutorials_plot_b_targetsoptions.py`

.. raw:: html

    </div>


.. toctree::
   :hidden:

   /_auto_tutorials/plot_b_targetsoptions

.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="The particle systems are defined in a rectangular box whose dimensions are specified by the att...">

.. only:: html

 .. figure:: /_auto_tutorials/images/thumb/sphx_glr_plot_c_boundaryconditions_thumb.png
     :alt: Tutorial 03: Boundary conditions

     :ref:`sphx_glr__auto_tutorials_plot_c_boundaryconditions.py`

.. raw:: html

    </div>


.. toctree::
   :hidden:

   /_auto_tutorials/plot_c_boundaryconditions

.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="In many cases, the interaction radius R is much smaller than the size of the domain. Consequent...">

.. only:: html

 .. figure:: /_auto_tutorials/images/thumb/sphx_glr_plot_d_blocksparse_thumb.png
     :alt: Tutorial 04: Block-sparse reduction

     :ref:`sphx_glr__auto_tutorials_plot_d_blocksparse.py`

.. raw:: html

    </div>


.. toctree::
   :hidden:

   /_auto_tutorials/plot_d_blocksparse

.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Simulating swarming models requires expensive mean-field convolution operations of the form: ">

.. only:: html

 .. figure:: /_auto_tutorials/images/thumb/sphx_glr_plot_e_kernels_thumb.png
     :alt: Tutorial 05: Kernels and averages

     :ref:`sphx_glr__auto_tutorials_plot_e_kernels.py`

.. raw:: html

    </div>


.. toctree::
   :hidden:

   /_auto_tutorials/plot_e_kernels
.. raw:: html

    <div class="sphx-glr-clear"></div>



.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-gallery


  .. container:: sphx-glr-download sphx-glr-download-python

    :download:`Download all examples in Python source code: _auto_tutorials_python.zip </_auto_tutorials/_auto_tutorials_python.zip>`



  .. container:: sphx-glr-download sphx-glr-download-jupyter

    :download:`Download all examples in Jupyter notebooks: _auto_tutorials_jupyter.zip </_auto_tutorials/_auto_tutorials_jupyter.zip>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_

.. DO NOT EDIT.
.. THIS FILE WAS AUTOMATICALLY GENERATED BY SPHINX-GALLERY.
.. TO MAKE CHANGES, EDIT THE SOURCE PYTHON FILE:
.. "_auto_tutorials/plot_b_targetsoptions.py"
.. LINE NUMBERS ARE GIVEN BELOW.

.. only:: html

    .. note::
        :class: sphx-glr-download-link-note

        Click :ref:`here <sphx_glr_download__auto_tutorials_plot_b_targetsoptions.py>`
        to download the full example code

.. rst-class:: sphx-glr-example-title

.. _sphx_glr__auto_tutorials_plot_b_targetsoptions.py:


Tutorial 02: Targets and options
=======================================

Many useful **local averages** (see :ref:`tuto_averages`) are pre-defined and may be directly used in the simulation of new or classical models. This tutorial showcases the basic usage of the **target methods** to simulate the variants of the Vicsek model. 

.. GENERATED FROM PYTHON SOURCE LINES 10-24

An example: variants of the Vicsek model
------------------------------------------

In its most abstract form, the Vicsek model reads: 

.. math:: 

        \mathrm{d}X^i_t = c_0 V^i_t \mathrm{d} t

.. math::

        \mathrm{d}V^i_t = \sigma\mathsf{P}(V^i_t)\circ ( J^i_t \mathrm{d}t + \mathrm{d} B^i_t),

where :math:`c_0` is the speed, :math:`\mathsf{P}(v)` is the orthogonal projection on the orthogonal plane to :math:`v` and :math:`\sigma` is the diffusion coefficient. The drift coefficient :math:`J^i_t` is called a **target**. In synchronous and asynchronous Vicsek models, the target is the center of the sampling distribution. A comparison of classical targets is shown below.

.. GENERATED FROM PYTHON SOURCE LINES 27-29

First, some standard imports...


.. GENERATED FROM PYTHON SOURCE LINES 29-40

.. code-block:: default


    import time 
    import torch
    import pprint
    from matplotlib import pyplot as plt
    from sisyphe.models import Vicsek
    from sisyphe.display import display_kinetic_particles

    use_cuda = torch.cuda.is_available()
    dtype = torch.cuda.FloatTensor if use_cuda else torch.FloatTensor








.. GENERATED FROM PYTHON SOURCE LINES 41-52

The classical non-normalised Vicsek model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The target is: 

.. math:: 

        J^i_t = \kappa\frac{\sum_{j=1}^N K(|X^j_t-X^i_t|)V^j_t}{|\sum_{j=1}^N K(|X^j_t-X^i_t|)V^j_t|}.

The parameters of the model... 


.. GENERATED FROM PYTHON SOURCE LINES 52-68

.. code-block:: default



    N = 100000
    L = 100 

    pos = L*torch.rand((N,2)).type(dtype)
    vel = torch.randn(N,2).type(dtype)
    vel = vel/torch.norm(vel,dim=1).reshape((N,1))

    R = 3.
    c = 3.
    nu = 5.
    sigma = 1.

    dt = .01








.. GENERATED FROM PYTHON SOURCE LINES 69-71

The choice of the target is implemented in the keyword argument ``variant``. For the classical normalised target, it is given by the following dictionary.


.. GENERATED FROM PYTHON SOURCE LINES 71-88

.. code-block:: default


    variant = {"name" : "normalised", "parameters" : {}}

    simu = Vicsek(
        pos = pos.detach().clone(),
        vel = vel.detach().clone(), 
        v = c, 
        sigma = sigma, 
        nu = nu, 
        interaction_radius = R,
        box_size = L,
        dt = dt,
        variant = variant,
        block_sparse_reduction = True,
        number_of_cells = 40**2)









.. GENERATED FROM PYTHON SOURCE LINES 89-90

Finally run the simulation over 300 units of time.... 

.. GENERATED FROM PYTHON SOURCE LINES 90-97

.. code-block:: default


    frames = [0., 5., 10., 30., 42., 71., 100, 124, 161, 206, 257, 300]

    s = time.time()
    it, op = display_kinetic_particles(simu, frames, order=True)
    e = time.time()




.. rst-class:: sphx-glr-horizontal


    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_b_targetsoptions_001.png
          :alt: Vicsek(normalised)  Parameters: N=100000 ; R=3.0 ; nu=5.0 ; sigma=1.0 ; v=3.0  Time=0.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_b_targetsoptions_002.png
          :alt: Vicsek(normalised)  Parameters: N=100000 ; R=3.0 ; nu=5.0 ; sigma=1.0 ; v=3.0  Time=5.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_b_targetsoptions_003.png
          :alt: Vicsek(normalised)  Parameters: N=100000 ; R=3.0 ; nu=5.0 ; sigma=1.0 ; v=3.0  Time=10.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_b_targetsoptions_004.png
          :alt: Vicsek(normalised)  Parameters: N=100000 ; R=3.0 ; nu=5.0 ; sigma=1.0 ; v=3.0  Time=30.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_b_targetsoptions_005.png
          :alt: Vicsek(normalised)  Parameters: N=100000 ; R=3.0 ; nu=5.0 ; sigma=1.0 ; v=3.0  Time=42.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_b_targetsoptions_006.png
          :alt: Vicsek(normalised)  Parameters: N=100000 ; R=3.0 ; nu=5.0 ; sigma=1.0 ; v=3.0  Time=71.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_b_targetsoptions_007.png
          :alt: Vicsek(normalised)  Parameters: N=100000 ; R=3.0 ; nu=5.0 ; sigma=1.0 ; v=3.0  Time=100.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_b_targetsoptions_008.png
          :alt: Vicsek(normalised)  Parameters: N=100000 ; R=3.0 ; nu=5.0 ; sigma=1.0 ; v=3.0  Time=124.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_b_targetsoptions_009.png
          :alt: Vicsek(normalised)  Parameters: N=100000 ; R=3.0 ; nu=5.0 ; sigma=1.0 ; v=3.0  Time=161.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_b_targetsoptions_010.png
          :alt: Vicsek(normalised)  Parameters: N=100000 ; R=3.0 ; nu=5.0 ; sigma=1.0 ; v=3.0  Time=206.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_b_targetsoptions_011.png
          :alt: Vicsek(normalised)  Parameters: N=100000 ; R=3.0 ; nu=5.0 ; sigma=1.0 ; v=3.0  Time=257.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_b_targetsoptions_012.png
          :alt: Vicsek(normalised)  Parameters: N=100000 ; R=3.0 ; nu=5.0 ; sigma=1.0 ; v=3.0  Time=300.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_b_targetsoptions_013.png
          :alt: plot b targetsoptions
          :class: sphx-glr-multi-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none


    Progress:0%
    Progress:1%
    Progress:2%
    Progress:3%
    Progress:4%
    Progress:5%
    Progress:6%
    Progress:7%
    Progress:8%
    Progress:9%
    Progress:10%
    Progress:11%
    Progress:12%
    Progress:13%
    Progress:14%
    Progress:15%
    Progress:16%
    Progress:17%
    Progress:18%
    Progress:19%
    Progress:20%
    Progress:21%
    Progress:22%
    Progress:23%
    Progress:24%
    Progress:25%
    Progress:26%
    Progress:27%
    Progress:28%
    Progress:29%
    Progress:30%
    Progress:31%
    Progress:32%
    Progress:33%
    Progress:34%
    Progress:35%
    Progress:36%
    Progress:37%
    Progress:38%
    Progress:39%
    Progress:40%
    Progress:41%
    Progress:42%
    Progress:43%
    Progress:44%
    Progress:45%
    Progress:46%
    Progress:47%
    Progress:48%
    Progress:49%
    Progress:50%
    Progress:51%
    Progress:52%
    Progress:53%
    Progress:54%
    Progress:55%
    Progress:56%
    Progress:57%
    Progress:58%
    Progress:59%
    Progress:60%
    Progress:61%
    Progress:62%
    Progress:63%
    Progress:64%
    Progress:65%
    Progress:66%
    Progress:67%
    Progress:68%
    Progress:69%
    Progress:70%
    Progress:71%
    Progress:72%
    Progress:73%
    Progress:74%
    Progress:75%
    Progress:76%
    Progress:77%
    Progress:78%
    Progress:79%
    Progress:80%
    Progress:81%
    Progress:82%
    Progress:83%
    Progress:84%
    Progress:85%
    Progress:86%
    Progress:87%
    Progress:88%
    Progress:89%
    Progress:90%
    Progress:91%
    Progress:92%
    Progress:93%
    Progress:94%
    Progress:95%
    Progress:96%
    Progress:97%
    Progress:98%
    Progress:99%



.. GENERATED FROM PYTHON SOURCE LINES 98-99

Print the total simulation time and the average time per iteration. 

.. GENERATED FROM PYTHON SOURCE LINES 99-103

.. code-block:: default


    print('Total time: '+str(e-s)+' seconds')
    print('Average time per iteration: '+str((e-s)/simu.iteration)+' seconds')





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Total time: 145.67891550064087 seconds
    Average time per iteration: 0.004855963850021362 seconds




.. GENERATED FROM PYTHON SOURCE LINES 104-105

Plot the histogram of the angles of the directions of motion. 

.. GENERATED FROM PYTHON SOURCE LINES 105-112

.. code-block:: default


    angle = torch.atan2(simu.vel[:,1],simu.vel[:,0])
    angle = angle.cpu().numpy()
    h = plt.hist(angle, bins=1000)
    plt.xlabel("angle")
    plt.show()




.. image:: /_auto_tutorials/images/sphx_glr_plot_b_targetsoptions_014.png
    :alt: plot b targetsoptions
    :class: sphx-glr-single-img





.. GENERATED FROM PYTHON SOURCE LINES 113-114

After an initial clustering phase, the system self-organizes into a uniform flock. 

.. GENERATED FROM PYTHON SOURCE LINES 116-126

Non-normalised Vicsek model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The target is: 

.. math:: 

        J^i_t = \frac{\frac{1}{N}\sum_{j=1}^N K(|X^j_t-X^i_t|)V^j_t}{\frac{1}{\kappa}+\frac{1}{\kappa_0}|\frac{1}{N}\sum_{j=1}^N K(|X^j_t-X^i_t|)V^j_t|}.

Define the corresponding dictionary...

.. GENERATED FROM PYTHON SOURCE LINES 126-145

.. code-block:: default


    kappa_0 = 15.

    variant = {"name" : "max_kappa", "parameters" : {"kappa_max" : kappa_0}}

    simu = Vicsek(
        pos = pos.detach().clone(),
        vel = vel.detach().clone(), 
        v = c, 
        sigma = sigma, 
        nu = nu, 
        interaction_radius = R,
        box_size = L,
        dt = dt,
        variant = variant,
        block_sparse_reduction = True,
        number_of_cells = 40**2)









.. GENERATED FROM PYTHON SOURCE LINES 146-147

Finally run the simulation over 300 units of time.... 

.. GENERATED FROM PYTHON SOURCE LINES 147-154

.. code-block:: default


    frames = [0., 5., 10., 30., 42., 71., 100, 124, 161, 206, 257, 300]

    s = time.time()
    it, op = display_kinetic_particles(simu, frames, order=True)
    e = time.time()




.. rst-class:: sphx-glr-horizontal


    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_b_targetsoptions_015.png
          :alt: Vicsek(max_kappa)  Parameters: N=100000 ; R=3.0 ; nu=5.0 ; sigma=1.0 ; v=3.0  Time=0.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_b_targetsoptions_016.png
          :alt: Vicsek(max_kappa)  Parameters: N=100000 ; R=3.0 ; nu=5.0 ; sigma=1.0 ; v=3.0  Time=5.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_b_targetsoptions_017.png
          :alt: Vicsek(max_kappa)  Parameters: N=100000 ; R=3.0 ; nu=5.0 ; sigma=1.0 ; v=3.0  Time=10.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_b_targetsoptions_018.png
          :alt: Vicsek(max_kappa)  Parameters: N=100000 ; R=3.0 ; nu=5.0 ; sigma=1.0 ; v=3.0  Time=30.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_b_targetsoptions_019.png
          :alt: Vicsek(max_kappa)  Parameters: N=100000 ; R=3.0 ; nu=5.0 ; sigma=1.0 ; v=3.0  Time=42.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_b_targetsoptions_020.png
          :alt: Vicsek(max_kappa)  Parameters: N=100000 ; R=3.0 ; nu=5.0 ; sigma=1.0 ; v=3.0  Time=71.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_b_targetsoptions_021.png
          :alt: Vicsek(max_kappa)  Parameters: N=100000 ; R=3.0 ; nu=5.0 ; sigma=1.0 ; v=3.0  Time=100.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_b_targetsoptions_022.png
          :alt: Vicsek(max_kappa)  Parameters: N=100000 ; R=3.0 ; nu=5.0 ; sigma=1.0 ; v=3.0  Time=124.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_b_targetsoptions_023.png
          :alt: Vicsek(max_kappa)  Parameters: N=100000 ; R=3.0 ; nu=5.0 ; sigma=1.0 ; v=3.0  Time=161.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_b_targetsoptions_024.png
          :alt: Vicsek(max_kappa)  Parameters: N=100000 ; R=3.0 ; nu=5.0 ; sigma=1.0 ; v=3.0  Time=206.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_b_targetsoptions_025.png
          :alt: Vicsek(max_kappa)  Parameters: N=100000 ; R=3.0 ; nu=5.0 ; sigma=1.0 ; v=3.0  Time=257.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_b_targetsoptions_026.png
          :alt: Vicsek(max_kappa)  Parameters: N=100000 ; R=3.0 ; nu=5.0 ; sigma=1.0 ; v=3.0  Time=300.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_b_targetsoptions_027.png
          :alt: plot b targetsoptions
          :class: sphx-glr-multi-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none


    Progress:0%
    Progress:1%
    Progress:2%
    Progress:3%
    Progress:4%
    Progress:5%
    Progress:6%
    Progress:7%
    Progress:8%
    Progress:9%
    Progress:10%
    Progress:11%
    Progress:12%
    Progress:13%
    Progress:14%
    Progress:15%
    Progress:16%
    Progress:17%
    Progress:18%
    Progress:19%
    Progress:20%
    Progress:21%
    Progress:22%
    Progress:23%
    Progress:24%
    Progress:25%
    Progress:26%
    Progress:27%
    Progress:28%
    Progress:29%
    Progress:30%
    Progress:31%
    Progress:32%
    Progress:33%
    Progress:34%
    Progress:35%
    Progress:36%
    Progress:37%
    Progress:38%
    Progress:39%
    Progress:40%
    Progress:41%
    Progress:42%
    Progress:43%
    Progress:44%
    Progress:45%
    Progress:46%
    Progress:47%
    Progress:48%
    Progress:49%
    Progress:50%
    Progress:51%
    Progress:52%
    Progress:53%
    Progress:54%
    Progress:55%
    Progress:56%
    Progress:57%
    Progress:58%
    Progress:59%
    Progress:60%
    Progress:61%
    Progress:62%
    Progress:63%
    Progress:64%
    Progress:65%
    Progress:66%
    Progress:67%
    Progress:68%
    Progress:69%
    Progress:70%
    Progress:71%
    Progress:72%
    Progress:73%
    Progress:74%
    Progress:75%
    Progress:76%
    Progress:77%
    Progress:78%
    Progress:79%
    Progress:80%
    Progress:81%
    Progress:82%
    Progress:83%
    Progress:84%
    Progress:85%
    Progress:86%
    Progress:87%
    Progress:88%
    Progress:89%
    Progress:90%
    Progress:91%
    Progress:92%
    Progress:93%
    Progress:94%
    Progress:95%
    Progress:96%
    Progress:97%
    Progress:98%
    Progress:99%



.. GENERATED FROM PYTHON SOURCE LINES 155-156

Print the total simulation time and the average time per iteration. 

.. GENERATED FROM PYTHON SOURCE LINES 156-160

.. code-block:: default


    print('Total time: '+str(e-s)+' seconds')
    print('Average time per iteration: '+str((e-s)/simu.iteration)+' seconds')





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Total time: 156.8030505180359 seconds
    Average time per iteration: 0.005226768350601196 seconds




.. GENERATED FROM PYTHON SOURCE LINES 161-162

Plot the histogram of the angles of the directions of motion. 

.. GENERATED FROM PYTHON SOURCE LINES 162-169

.. code-block:: default


    angle = torch.atan2(simu.vel[:,1],simu.vel[:,0])
    angle = angle.cpu().numpy()
    h = plt.hist(angle, bins=1000)
    plt.xlabel("angle")
    plt.show()




.. image:: /_auto_tutorials/images/sphx_glr_plot_b_targetsoptions_028.png
    :alt: plot b targetsoptions
    :class: sphx-glr-single-img





.. GENERATED FROM PYTHON SOURCE LINES 170-171

The system self-organizes into a strongly clustered flock with band-like structures. 

.. GENERATED FROM PYTHON SOURCE LINES 174-189

Nematic Vicsek model 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The target is: 

.. math:: 

        J^i_t = \kappa (V^i_t\cdot \overline{\Omega}^i_t)\overline{\Omega}^i_t,

where :math:`\overline{\Omega}^i_t` is any unit eigenvector associated to the maximal eigenvalue of the average Q-tensor:

.. math::

        Q^i_t = \frac{1}{N}\sum_{j=1}^N K(|X^j_t-X^i_t|){\left(V^j_t\otimes V^j_t - \frac{1}{d} I_d\right)}.

Define the corresponding dictionary...

.. GENERATED FROM PYTHON SOURCE LINES 189-205

.. code-block:: default


    variant = {"name" : "nematic", "parameters" : {}}

    simu = Vicsek(
        pos = pos.detach().clone(),
        vel = vel.detach().clone(), 
        v = c, 
        sigma = sigma, 
        nu = nu, 
        interaction_radius = R,
        box_size = L,
        dt = dt,
        variant = variant,
        block_sparse_reduction = True, 
        number_of_cells = 40**2)








.. GENERATED FROM PYTHON SOURCE LINES 206-207

Finally run the simulation over 100 units of time... The color code indicates the angle of the direction of motion between :math:`-\pi` and :math:`\pi`. 

.. GENERATED FROM PYTHON SOURCE LINES 207-214

.. code-block:: default


    frames = [0., 5., 10., 30., 42., 71., 100]

    s = time.time()
    it, op = display_kinetic_particles(simu, frames, order=True, color=True)
    e = time.time()




.. rst-class:: sphx-glr-horizontal


    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_b_targetsoptions_029.png
          :alt: Vicsek(nematic)  Parameters: N=100000 ; R=3.0 ; nu=5.0 ; sigma=1.0 ; v=3.0  Time=0.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_b_targetsoptions_030.png
          :alt: Vicsek(nematic)  Parameters: N=100000 ; R=3.0 ; nu=5.0 ; sigma=1.0 ; v=3.0  Time=5.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_b_targetsoptions_031.png
          :alt: Vicsek(nematic)  Parameters: N=100000 ; R=3.0 ; nu=5.0 ; sigma=1.0 ; v=3.0  Time=10.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_b_targetsoptions_032.png
          :alt: Vicsek(nematic)  Parameters: N=100000 ; R=3.0 ; nu=5.0 ; sigma=1.0 ; v=3.0  Time=30.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_b_targetsoptions_033.png
          :alt: Vicsek(nematic)  Parameters: N=100000 ; R=3.0 ; nu=5.0 ; sigma=1.0 ; v=3.0  Time=42.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_b_targetsoptions_034.png
          :alt: Vicsek(nematic)  Parameters: N=100000 ; R=3.0 ; nu=5.0 ; sigma=1.0 ; v=3.0  Time=71.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_b_targetsoptions_035.png
          :alt: Vicsek(nematic)  Parameters: N=100000 ; R=3.0 ; nu=5.0 ; sigma=1.0 ; v=3.0  Time=100.0
          :class: sphx-glr-multi-img

    *

      .. image:: /_auto_tutorials/images/sphx_glr_plot_b_targetsoptions_036.png
          :alt: plot b targetsoptions
          :class: sphx-glr-multi-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none


    Progress:0%
    Progress:1%
    Progress:2%
    Progress:3%
    Progress:4%
    Progress:5%
    Progress:6%
    Progress:7%
    Progress:8%
    Progress:9%
    Progress:10%
    Progress:11%
    Progress:12%
    Progress:13%
    Progress:14%
    Progress:15%
    Progress:16%
    Progress:17%
    Progress:18%
    Progress:19%
    Progress:20%
    Progress:21%
    Progress:22%
    Progress:23%
    Progress:24%
    Progress:25%
    Progress:26%
    Progress:27%
    Progress:28%
    Progress:29%
    Progress:30%
    Progress:31%
    Progress:32%
    Progress:33%
    Progress:34%
    Progress:35%
    Progress:36%
    Progress:37%
    Progress:38%
    Progress:39%
    Progress:40%
    Progress:41%
    Progress:42%
    Progress:43%
    Progress:44%
    Progress:45%
    Progress:46%
    Progress:47%
    Progress:48%
    Progress:49%
    Progress:50%
    Progress:51%
    Progress:52%
    Progress:53%
    Progress:54%
    Progress:55%
    Progress:56%
    Progress:57%
    Progress:58%
    Progress:59%
    Progress:60%
    Progress:61%
    Progress:62%
    Progress:63%
    Progress:64%
    Progress:65%
    Progress:66%
    Progress:67%
    Progress:68%
    Progress:69%
    Progress:70%
    Progress:71%
    Progress:72%
    Progress:73%
    Progress:74%
    Progress:75%
    Progress:76%
    Progress:77%
    Progress:78%
    Progress:79%
    Progress:80%
    Progress:81%
    Progress:82%
    Progress:83%
    Progress:84%
    Progress:85%
    Progress:86%
    Progress:87%
    Progress:88%
    Progress:89%
    Progress:90%
    Progress:91%
    Progress:92%
    Progress:93%
    Progress:94%
    Progress:95%
    Progress:96%
    Progress:97%
    Progress:98%
    Progress:99%



.. GENERATED FROM PYTHON SOURCE LINES 215-216

Print the total simulation time and the average time per iteration. 

.. GENERATED FROM PYTHON SOURCE LINES 216-220

.. code-block:: default


    print('Total time: '+str(e-s)+' seconds')
    print('Average time per iteration: '+str((e-s)/simu.iteration)+' seconds')





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Total time: 56.330119132995605 seconds
    Average time per iteration: 0.00563301191329956 seconds




.. GENERATED FROM PYTHON SOURCE LINES 221-222

Plot the histogram of the angles of the directions of motion. 

.. GENERATED FROM PYTHON SOURCE LINES 222-231

.. code-block:: default


    # sphinx_gallery_thumbnail_number = -1

    angle = torch.atan2(simu.vel[:,1],simu.vel[:,0])
    angle = angle.cpu().numpy()
    h = plt.hist(angle, bins=1000)
    plt.xlabel("angle")
    plt.show()




.. image:: /_auto_tutorials/images/sphx_glr_plot_b_targetsoptions_037.png
    :alt: plot b targetsoptions
    :class: sphx-glr-single-img





.. GENERATED FROM PYTHON SOURCE LINES 232-233

There are two modes separated by an angle :math:`\pi` which indicates that two groups of equal size are moving in opposite direction. This a *nematic flock*. 

.. GENERATED FROM PYTHON SOURCE LINES 236-240

The target dictionary 
--------------------------------

Several other targets are implemented. The complete list of available targets can be found in the **dictionary of targets** :attr:`target_method <sisyphe.particles.Particles.target_method>`.

.. GENERATED FROM PYTHON SOURCE LINES 240-243

.. code-block:: default


    pprint.pprint(simu.target_method)





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    {'max_kappa': <bound method KineticParticles.max_kappa of <sisyphe.models.Vicsek object at 0x7efb4e14e970>>,
     'mean_field': <bound method KineticParticles.mean_field of <sisyphe.models.Vicsek object at 0x7efb4e14e970>>,
     'morse': <bound method Particles.morse_target of <sisyphe.models.Vicsek object at 0x7efb4e14e970>>,
     'motsch_tadmor': <bound method KineticParticles.motsch_tadmor of <sisyphe.models.Vicsek object at 0x7efb4e14e970>>,
     'nematic': <bound method KineticParticles.nematic of <sisyphe.models.Vicsek object at 0x7efb4e14e970>>,
     'normalised': <bound method KineticParticles.normalised of <sisyphe.models.Vicsek object at 0x7efb4e14e970>>,
     'overlapping_repulsion': <bound method Particles.overlapping_repulsion_target of <sisyphe.models.Vicsek object at 0x7efb4e14e970>>,
     'quadratic_potential': <bound method Particles.quadratic_potential_target of <sisyphe.models.Vicsek object at 0x7efb4e14e970>>}




.. GENERATED FROM PYTHON SOURCE LINES 244-245

Custom targets can be added to the dictionary of targets using the method :meth:`add_target_method() <sisyphe.particles.Particles.add_target_method>`. Then a target can be readily used in a simulation using the method :meth:`compute_target() <sisyphe.particles.Particles.compute_target>`. 

.. GENERATED FROM PYTHON SOURCE LINES 248-252

Options
---------------

Customized targets can also be defined by applying an **option** which modifies an existing target. See :ref:`examplemill`.  


.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 6 minutes  25.941 seconds)


.. _sphx_glr_download__auto_tutorials_plot_b_targetsoptions.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download sphx-glr-download-python

     :download:`Download Python source code: plot_b_targetsoptions.py <plot_b_targetsoptions.py>`



  .. container:: sphx-glr-download sphx-glr-download-jupyter

     :download:`Download Jupyter notebook: plot_b_targetsoptions.ipynb <plot_b_targetsoptions.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_

.. DO NOT EDIT.
.. THIS FILE WAS AUTOMATICALLY GENERATED BY SPHINX-GALLERY.
.. TO MAKE CHANGES, EDIT THE SOURCE PYTHON FILE:
.. "_auto_tutorials/plot_e_kernels.py"
.. LINE NUMBERS ARE GIVEN BELOW.

.. only:: html

    .. note::
        :class: sphx-glr-download-link-note

        Click :ref:`here <sphx_glr_download__auto_tutorials_plot_e_kernels.py>`
        to download the full example code

.. rst-class:: sphx-glr-example-title

.. _sphx_glr__auto_tutorials_plot_e_kernels.py:


.. _tuto_averages:

Tutorial 05: Kernels and averages
============================================

Simulating swarming models requires expensive mean-field convolution operations of the form: 

.. math::
    
    J^i = \frac{1}{N}\sum_{j=1}^N K(|X^j-X^i|) U^j,
    
for :math:`1\leq i\leq N`, where :math:`(X^i)_{1\leq i \leq N}` are the positions of the particles, :math:`(U^j)_{1\leq j\leq N}` are given vectors and :math:`K` is an **observation kernel**. Typically, :math:`K(|X^i-X^j|)` is equal to 1 if :math:`X^i` and :math:`X^j` are at distance smaller than a fixed interaction distance and 0 otherwise. Other kernels are defined in the module :mod:`sisyphe.kernels`. Below, we show a simple application case. 

.. GENERATED FROM PYTHON SOURCE LINES 18-23

Linear local averages
---------------------------

First, some standard imports...


.. GENERATED FROM PYTHON SOURCE LINES 23-32

.. code-block:: default


    import time 
    import math
    import torch
    from matplotlib import pyplot as plt

    use_cuda = torch.cuda.is_available()
    dtype = torch.cuda.FloatTensor if use_cuda else torch.FloatTensor








.. GENERATED FROM PYTHON SOURCE LINES 33-35

Let the :math:`N` particles be uniformly scattered in a box of size :math:`L` with interaction radius  :math:`R`.


.. GENERATED FROM PYTHON SOURCE LINES 35-43

.. code-block:: default


    N = 100000
    L = 1. 
    R = .15

    pos = L*torch.rand((N,2)).type(dtype)









.. GENERATED FROM PYTHON SOURCE LINES 44-45

We can also assume that the particles have a bounded cone of vision around an axis (defined by a unit vector). The default behaviour is a full vision angle equal to :math:`2\pi` in which case the axis is a :data:`None` object. Here we take a cone of vision with angle :math:`\pi/2` around an axis which is sampled uniformly. For the :class:`sisyphe.particles.KineticParticles`, the default axis is the velocity. 

.. GENERATED FROM PYTHON SOURCE LINES 45-50

.. code-block:: default


    angle = math.pi/2
    axis = torch.randn(N,2).type(dtype)
    axis = axis/torch.norm(axis,dim=1).reshape((N,1))








.. GENERATED FROM PYTHON SOURCE LINES 51-52

Let us create an instance of a particle system with these parameters. 

.. GENERATED FROM PYTHON SOURCE LINES 52-62

.. code-block:: default


    from sisyphe.particles import Particles 

    particles = Particles(
        pos = pos,
        interaction_radius = R,
        box_size = L,
        vision_angle = angle,
        axis = axis)








.. GENERATED FROM PYTHON SOURCE LINES 63-65

.. note::
        By default, the system and the operations below are defined with periodic boundary conditions. 

.. GENERATED FROM PYTHON SOURCE LINES 67-72

As a simple application, we can compute the number of neighbours of each particle and print the number of neighbours of the first particle. This operation is already implemented in the method :func:`number_of_neighbours() <sisyphe.particles.Particles.number_of_neighbours>`. It simply corresponds to the average: 

.. math::

        N^i_\mathrm{neigh} = \sum_{j=1}^N K(|X^j-X^i|).

.. GENERATED FROM PYTHON SOURCE LINES 72-79

.. code-block:: default


    Nneigh = particles.number_of_neighbours()

    Nneigh0 = int(Nneigh[0].item())

    print("The first particle sees " + str(Nneigh0) + " other particles.")





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    The first particle sees 1754 other particles.




.. GENERATED FROM PYTHON SOURCE LINES 80-81

For custom objects, the mean-field average can be computed using the method :func:`linear_local_average() <sisyphe.particles.Particles.linear_local_average>`. As an example, let us compute the center of mass of the neighbours of each particle. First we define the quantity :math:`U` that we want to average. Here, since we are working on a torus, there are two: the sine and the cosine of the spatial coordinates. 

.. GENERATED FROM PYTHON SOURCE LINES 81-85

.. code-block:: default


    cos_pos = torch.cos((2*math.pi / L) * particles.pos)
    sin_pos = torch.sin((2*math.pi / L) * particles.pos)








.. GENERATED FROM PYTHON SOURCE LINES 86-87

Then we compute the two mean field averages, i.e. the standard convolution over the :math:`N` particles. The center of mass along each dimension is the argument of the complex number whose coordinates are the average cosine and sine. 

.. GENERATED FROM PYTHON SOURCE LINES 87-97

.. code-block:: default


    average_cos, average_sin = particles.linear_local_average(cos_pos, sin_pos)
    center_x = torch.atan2(average_sin[:,0], average_cos[:,0])
    center_x = (L / (2*math.pi)) * torch.remainder(center_x, 2*math.pi)
    center_y = torch.atan2(average_sin[:,1], average_cos[:,1])
    center_y = (L / (2*math.pi)) * torch.remainder(center_y, 2*math.pi)

    center_of_mass = torch.cat((center_x.reshape((N,1)), center_y.reshape((N,1))),
                                dim=1)





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    [pyKeOps] Compiling libKeOpstorch3dcd0c2195 in /data/and18/.cache/pykeops-1.5-cpython-38:
           formula: Sum_Reduction(((Step((Var(5,1,2) - Sum(Square((((Var(0,2,1) - Var(1,2,0)) + (Step(((Minus(Var(2,2,2)) / Var(3,1,2)) - (Var(0,2,1) - Var(1,2,0)))) * Var(2,2,2))) - (Step(((Var(0,2,1) - Var(1,2,0)) - (Var(2,2,2) / Var(4,1,2)))) * Var(2,2,2))))))) * Step((Var(7,1,2) + (Sum(((((Var(0,2,1) - Var(1,2,0)) + (Step(((Minus(Var(2,2,2)) / Var(3,1,2)) - (Var(0,2,1) - Var(1,2,0)))) * Var(2,2,2))) - (Step(((Var(0,2,1) - Var(1,2,0)) - (Var(2,2,2) / Var(4,1,2)))) * Var(2,2,2))) * Var(6,2,0))) / Sqrt(Sum(Square((((Var(0,2,1) - Var(1,2,0)) + (Step(((Minus(Var(2,2,2)) / Var(3,1,2)) - (Var(0,2,1) - Var(1,2,0)))) * Var(2,2,2))) - (Step(((Var(0,2,1) - Var(1,2,0)) - (Var(2,2,2) / Var(4,1,2)))) * Var(2,2,2)))))))))) * Var(8,4,1)),0)
           aliases: Var(0,2,1); Var(1,2,0); Var(2,2,2); Var(3,1,2); Var(4,1,2); Var(5,1,2); Var(6,2,0); Var(7,1,2); Var(8,4,1); 
           dtype  : float32
    ... 
    Done.




.. GENERATED FROM PYTHON SOURCE LINES 98-99

In the method :func:`linear_local_average() <sisyphe.particles.Particles.linear_local_average>`, the default observation kernel is a :class:`LazyTensor` of size :math:`(N,N)` whose :math:`(i,j)` component is equal to 1 when particle :math:`j` belongs to the cone of vision of particle :math:`i` and 0 otherwise. To retrieve the indexes of the particles which belong to the cone of vision of the first particle, we can use the `K-nearest-neighbours reduction <https://www.kernel-operations.io/keops/_auto_tutorials/knn/plot_knn_mnist.html#sphx-glr-auto-tutorials-knn-plot-knn-mnist-py>`_ provided by the `KeOps <https://www.kernel-operations.io/keops/index.html>`_ library. 

.. GENERATED FROM PYTHON SOURCE LINES 99-118

.. code-block:: default


    from sisyphe.kernels import lazy_interaction_kernel

    interaction_kernel = lazy_interaction_kernel(
        particles.pos, 
        particles.pos, 
        particles.R,
        particles.L,
        boundary_conditions = particles.bc,
        vision_angle = particles.angle,
        axis = particles.axis)

    K_ij = 1. - interaction_kernel 

    neigh0 = K_ij.argKmin(Nneigh0, dim=1)[0]

    print("The indexes of the neighbours of the first particles are: ")
    print(neigh0)





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    [pyKeOps] Compiling libKeOpstorchaf7a666555 in /data/and18/.cache/pykeops-1.5-cpython-38:
           formula: ArgKMin_Reduction((Var(8,1,2) - (Step((Var(5,1,2) - Sum(Square((((Var(0,2,1) - Var(1,2,0)) + (Step(((Minus(Var(2,2,2)) / Var(3,1,2)) - (Var(0,2,1) - Var(1,2,0)))) * Var(2,2,2))) - (Step(((Var(0,2,1) - Var(1,2,0)) - (Var(2,2,2) / Var(4,1,2)))) * Var(2,2,2))))))) * Step((Var(7,1,2) + (Sum(((((Var(0,2,1) - Var(1,2,0)) + (Step(((Minus(Var(2,2,2)) / Var(3,1,2)) - (Var(0,2,1) - Var(1,2,0)))) * Var(2,2,2))) - (Step(((Var(0,2,1) - Var(1,2,0)) - (Var(2,2,2) / Var(4,1,2)))) * Var(2,2,2))) * Var(6,2,0))) / Sqrt(Sum(Square((((Var(0,2,1) - Var(1,2,0)) + (Step(((Minus(Var(2,2,2)) / Var(3,1,2)) - (Var(0,2,1) - Var(1,2,0)))) * Var(2,2,2))) - (Step(((Var(0,2,1) - Var(1,2,0)) - (Var(2,2,2) / Var(4,1,2)))) * Var(2,2,2))))))))))),1754,0)
           aliases: Var(0,2,1); Var(1,2,0); Var(2,2,2); Var(3,1,2); Var(4,1,2); Var(5,1,2); Var(6,2,0); Var(7,1,2); Var(8,1,2); 
           dtype  : float32
    ... 
    Done.
    The indexes of the neighbours of the first particles are: 
    tensor([    0,    29,   130,  ..., 99811, 99868, 99982], device='cuda:0')




.. GENERATED FROM PYTHON SOURCE LINES 119-120

Finally, a fancy display of what we have computed. We plot the full particle system in black, the first particle in orange, its neighbours in blue and the center of mass of the neighbours in red. 

.. GENERATED FROM PYTHON SOURCE LINES 120-142

.. code-block:: default


    xall = particles.pos[:,0].cpu().numpy()
    yall = particles.pos[:,1].cpu().numpy()

    x = particles.pos[neigh0,0].cpu().numpy()
    y = particles.pos[neigh0,1].cpu().numpy()

    x0 = particles.pos[0,0].item()
    y0 = particles.pos[0,1].item()

    xc = center_of_mass[0,0].item()
    yc = center_of_mass[0,1].item()


    fig, ax = plt.subplots(figsize=(6,6))
    ax.scatter(xall, yall, s=.003, c='black')
    ax.scatter(x, y, s=.3)
    ax.scatter(x0, y0, s=24)
    ax.scatter(xc, yc, s=24, c='red')
    ax.axis([0, L, 0, L])
    ax.set_aspect("equal")




.. image:: /_auto_tutorials/images/sphx_glr_plot_e_kernels_001.png
    :alt: plot e kernels
    :class: sphx-glr-single-img





.. GENERATED FROM PYTHON SOURCE LINES 143-162

Nonlinear averages
---------------------------------

In some cases, we need to compute a **nonlinear average** of the form 

.. math::

        J^i = \frac{1}{N}\sum_{j=1}^N K(|X^j-X^i|) b(U^i,V^j),

where :math:`(U^i)_{1\leq i \leq N}` and :math:`(V^j)_{1\leq j \leq N}` are given vectors and :math:`b` is a given function. When the **binary formula** :math:`b` can be written as a :class:`LazyTensor`, this can be computed with the method :func:`nonlinear_local_average() <sisyphe.particles.Particles.nonlinear_local_average>`. 

For instance, let us compute the local mean square distance: 

.. math::

        J^i = \frac{\sum_{j=1}^N K(|X^j-X^i|) |X^j-X^i|^2}{\sum_{j=1}^N K(|X^j-X^i|)}.

In this case, we can use the function :func:`sisyphe.kernels.lazy_xy_matrix` to define a custom binary formula. Given two vectors :math:`X=(X^i)_{1\leq i\leq M}` and :math:`Y = (Y^j)_{1\leq j\leq N}`, respectively of sizes :math:`(M,d)` and :math:`(N,d)`, the :math:`XY` matrix is a :math:`(M,N,d)` LazyTensor whose :math:`(i,j,:)` component is the vector :math:`Y^j-X^i`. 


.. GENERATED FROM PYTHON SOURCE LINES 162-173

.. code-block:: default


    from sisyphe.kernels import lazy_xy_matrix 

    def b(x,y): 
        K_ij = lazy_xy_matrix(x,y,particles.L)
        return (K_ij ** 2).sum(-1)

    x = particles.pos
    y = particles.pos
    mean_square_dist = N/Nneigh.reshape((N,1)) * particles.nonlinear_local_average(b,x,y)








.. GENERATED FROM PYTHON SOURCE LINES 174-180

Since the particles are uniformly scattered in the box, the theoretical value is 

.. math::

        MSD_0 = \frac{\int_0^R \int_0^{\pi/2} r^3 \mathrm{d}r\mathrm{d}\theta}{\int_0^R \int_0^{\pi/2} r \mathrm{d}r\mathrm{d}\theta} = \frac{R^2}{2}


.. GENERATED FROM PYTHON SOURCE LINES 180-193

.. code-block:: default


    print("Theoretical value: " + str(R**2/2))
    print("Experimental value: " + str(mean_square_dist[0].item()))






    
    
    
    
    



.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Theoretical value: 0.01125
    Experimental value: 0.01103969756513834





.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 2 minutes  49.544 seconds)


.. _sphx_glr_download__auto_tutorials_plot_e_kernels.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download sphx-glr-download-python

     :download:`Download Python source code: plot_e_kernels.py <plot_e_kernels.py>`



  .. container:: sphx-glr-download sphx-glr-download-jupyter

     :download:`Download Jupyter notebook: plot_e_kernels.ipynb <plot_e_kernels.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
Kernels
===========================

.. currentmodule:: sisyphe.kernels
.. autosummary::

	lazy_xy_matrix
    lazy_interaction_kernel
    lazy_overlapping_kernel
    lazy_morse
    lazy_quadratic

.. automodule:: sisyphe.kernels
	:members:Display
===========================

.. currentmodule:: sisyphe.display
.. autosummary::

	save
    display_kinetic_particles
    scatter_particles
    
.. automodule:: sisyphe.display
	:members:Particle systems
===========================

.. currentmodule:: sisyphe.particles
.. autosummary::

	Particles
    KineticParticles
    BOParticles

.. automodule:: sisyphe.particles
	:members:Toolbox
===========================

.. currentmodule:: sisyphe.toolbox
.. autosummary::


.. automodule:: sisyphe.toolbox
	:members:Models
===========================

.. currentmodule:: sisyphe.models
.. autosummary::

	AsynchronousVicsek
    Vicsek
    SynchronousVicsek
    BOAsynchronousVicsek
    VolumeExclusion
    AttractionRepulsion
    

.. automodule:: sisyphe.models
	:members:Sampling
===========================

.. currentmodule:: sisyphe.sampling
.. autosummary::

	uniform_sphere_rand
    uniform_angle_rand
    vonmises_rand
    uniformball_unitquat_rand
    vonmises_quat_rand

.. automodule:: sisyphe.sampling
	:members: