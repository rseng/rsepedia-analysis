# QMCTorch

Pytorch Implementation of Real Space Quantum Monte Carlo Simulations of Molecular Systems

[![PyPI version](https://badge.fury.io/py/qmctorch.svg)](https://badge.fury.io/py/qmctorch)
[![Build Status](https://github.com/NLESC-JCER/QMCTorch/workflows/build/badge.svg)](https://github.com/NLESC-JCER/QMCTorch/actions)
[![Coverage Status](https://coveralls.io/repos/github/NLESC-JCER/QMCTorch/badge.svg?branch=master)](https://coveralls.io/github/NLESC-JCER/QMCTorch?branch=master)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/5d99212add2a4f0591adc6248fec258d)](https://www.codacy.com/manual/NicoRenaud/QMCTorch?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=NLESC-JCER/QMCTorch&amp;utm_campaign=Badge_Grade)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3780094.svg)](https://doi.org/10.5281/zenodo.3780094)

## Installation

Clone the repository and install the code from source or use the Python package manager:

`pip install qmctorch`

## Documentation 
https://qmctorch.readthedocs.io/en/latest/intro.html


## Disclaimer
QMCTorch is currently under developmement and most likely won't behave as expected 
---
name: Bug report
about: Something doesn't work like I expected.
title: ''
labels: bug
assignees: ''

---

Use your best judgment to provide a useful level of information. Depending on the nature of the issue, consider including, e.g.

- Which version of the software you're running
- Any logging or error output you see
---
name: Feature request
about: I would like a new feature to be included in the library
title: ''
labels: enhancement
assignees: ''

---

Tell us what would be a nice feature to make the **QMCTorch** library even better!
---
name: Help wanted
about: I need some help with the code
title: ''
labels: help wanted
assignees: ''

---
############################
Contributing guidelines
############################

We welcome any kind of contribution to our software, from simple comment or question to a full fledged `pull request <https://help.github.com/articles/about-pull-requests/>`_. Please read and follow our `Code of Conduct <CODE_OF_CONDUCT.rst>`_.

A contribution can be one of the following cases:

#. you have a question;
#. you think you may have found a bug (including unexpected behavior);
#. you want to make some kind of change to the code base (e.g. to fix a bug, to add a new feature, to update documentation).

The sections below outline the steps in each case.

You have a question
*******************

#. use the search functionality `here <https://github.com/NLESC-JCER/pyCHAMP/issues>`__ to see if someone already filed the same issue;
#. if your issue search did not yield any relevant results, make a new issue;
#. apply the "Question" label; apply other labels when relevant.

You think you may have found a bug
**********************************

#. use the search functionality `here <https://github.com/NLESC-JCER/pyCHAMP/issues>`__ to see if someone already filed the same issue;
#. if your issue search did not yield any relevant results, make a new issue, making sure to provide enough information to the rest of the community to understand the cause and context of the problem. Depending on the issue, you may want to include:
    - the `SHA hashcode <https://help.github.com/articles/autolinked-references-and-urls/#commit-shas>`_ of the commit that is causing your problem;
    - some identifying information (name and version number) for dependencies you're using;
    - information about the operating system;
#. apply relevant labels to the newly created issue.

You want to make some kind of change to the code base
*****************************************************

#. (**important**) announce your plan to the rest of the community *before you start working*. This announcement should be in the form of a (new) issue;
#. (**important**) wait until some kind of consensus is reached about your idea being a good idea;
#. if needed, fork the repository to your own Github profile and create your own feature branch off of the latest master commit. While working on your feature branch, make sure to stay up to date with the master branch by pulling in changes, possibly from the 'upstream' repository (follow the instructions `here <https://help.github.com/articles/configuring-a-remote-for-a-fork/>`__ and `here <https://help.github.com/articles/syncing-a-fork/>`__);
#. make sure the existing tests still work by running ``python setup.py test``;
#. add your own tests (if necessary);
#. update or expand the documentation;
#. `push <http://rogerdudler.github.io/git-guide/>`_ your feature branch to (your fork of) the pyCHAMP repository on GitHub;
#. create the pull request, e.g. following the instructions `here <https://help.github.com/articles/creating-a-pull-request/>`__.

In case you feel like you've made a valuable contribution, but you don't know how to write or run tests for it, or how to generate the documentation: don't let this discourage you from making the pull request; we can help you! Just go ahead and submit the pull request, but keep in mind that you might be asked to append additional commits to your pull request.
###############################################################################
Contributor Covenant Code of Conduct
###############################################################################

Our Pledge
**********

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, gender identity and expression, level of experience,
education, socio-economic status, nationality, personal appearance, race,
religion, or sexual identity and orientation.

Our Standards
*************

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

Our Responsibilities
********************

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

Scope
*****

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

Enforcement
***********

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at n.renaud@esciencecenter.nl. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

Attribution
***********

This Code of Conduct is adapted from the `Contributor Covenant <https://www.contributor-covenant.org>`_, version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html
##########
Change Log
##########

0.1.1 [Unreleased]
******************

Change
------

* Make `Horovod optional <https://github.com/horovod/horovod>`_ (#57)


Add
---
* GitHub Actions CI (#58)
* `Plams dependency <https://github.com/SCM-NV/PLAMS>`_ (#63)
Quantum Monte Carlo
========================

TODO : Quantum Monte Carlo a 10 min introduction.. QMCTorch documentation master file, created by
   sphinx-quickstart on Wed Apr 22 17:16:01 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Quantum Monte Carlo with Pytorch
================================================

.. toctree::
   :maxdepth: 1
   :caption: QMCTorch

   intro
   qmc
   qmctorch/molecule
   qmctorch/wavefunction
   qmctorch/sampler
   qmctorch/optimizer
   qmctorch/solver

.. toctree::
   :maxdepth: 1
   :caption: Tutorial

   tutorial/tutorial_sampling_traj
   tutorial/tutorial_correlation
   tutorial/tutorial_wf_opt
   tutorial/tutorial_jastrow
   tutorial/tutorial_backflow
   # tutorial/tutorial_geo_opt
   tutorial/tutorial_gpus
   tutorial/tutorial_hdf5

.. toctree::
   :maxdepth: 1
   :caption: Examples

   example/example_sp
   example/example_opt
   example/example_gpu
   example/example_horovod

.. toctree::
   :maxdepth: 1
   :caption: API

   source/modules.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
Installation
=========================

Via Python Package
-----------------------------

The latest release of QMCTorch can be installed usint the pypi package manager with :

``pip install qmctorch`` 

If you are planning to only use QMCTorch this is the quickest way to obtain the code


Via GitHub
-------------

For user who would like to contribute, the code is hosted on GitHub_ (https://github.com/NLESC-JCER/QMCTorch)

.. _GitHub: https://github.com/NLESC-JCER/QMCTorch

To install the code

 * clone the repository ``git clone https://github.com/NLESC-JCER/QMCTorch.git``
 * go there ``cd QMCTorch``
 * install the module ``pip insall -e ./``

You can then test the installation :

 * ``cd test``
 * ``pytest``


Single point calculation
==============================================

.. literalinclude:: ../../example/single_point/h2o_sampling.pyWave function optimization
==============================================

.. literalinclude:: ../../example/optimization/h2.pyGPU support
==============================================

.. literalinclude:: ../../example/gpu/h2.pyMulti-GPU support
==============================================

.. literalinclude:: ../../example/horovod/h2.pyGeometry Optimization
====================================

We present here a complete example on how to use QMCTorch on a H2O molecule.
As previously the firs task is to import all the modules needed

>>> from torch import optim
>>> from torch.optim import Adam
>>>  from qmctorch.wavefunction import SlaterJastrow
>>> from qmctorch.solver import SolverSlaterJastrow
>>> from qmctorch.samplerimport Metropolis
>>> from qmctorch.scf import Molecule
>>> from qmctorch.utils import plot_energy

We then define the molecule. We read here an xyz file of a water molecule
where the three atoms are on the same line.

>>> # define the molecule
>>> mol = Molecule(atom='water_line.xyz', unit='angs',
>>>                calculator='pyscf', basis='sto-3g')


The QMCNet wave function is defined from the molecule object. We only consider here the
ground state of the molecule in the CI expansion.

>>> # define the wave function
>>> wf = SlaterJastrow(mol, kinetic='jacobi',
>>>          configs='ground_state',
>>>           use_jastrow=True)

We use a Metropolis Hasting sampler with 100 walkers each performing 2000 steps.

>>> # sampler
>>> sampler = Metropolis(nwalkers=1000, nstep=2000, step_size=0.5,
>>>                      nelec=wf.nelec, ndim=wf.ndim,
>>>                      init=mol.domain('normal'),
>>>                      move={'type': 'one-elec', 'proba': 'normal'})

As an opimizer we use ADAM and define a simple linear scheduler.

>>> # optimizer
>>> opt = Adam(wf.parameters(), lr=0.005)
>>>
>>> # scheduler
>>> scheduler = optim.lr_scheduler.StepLR(opt, step_size=20, gamma=0.75)

We can now assemble the solver

>>> # solver
>>> solver = SolverSlaterJastrow(wf=wf,
>>>                        sampler=sampler,
>>>                        optimizer=opt,
>>>                        scheduler=scheduler)

To optimize the geometry of the molecule we must specify `task=geo_opt`. This will
freeze all the parameters of the wave function at the exception of th atomic coordinate.
We can then run the optimization here using 50 epochs, the energy as loss function and the low variance expression
of the gradients.

>>> # optimize the geometry
>>> solver.configure(task='geo_opt')
>>> obs = solver.run(50, loss='energy', grad='manual')
>>> solver.save_traj('h2o_traj.xyz')

We can then plot the energy

>>> # plot the data
>>> plot_energy(obs.local_energy)

.. image:: ../../pics/h2_go_opt.png


GPUs and multi-GPUs Support
==============================================

.. warning::
    The use of GPU and mutli-GPU is under developpement and hasn't been
    thoroughly tested yet. Proceed with caution !


Using pytorch as a backend, QMCTorch can leverage GPU cards available on your hardware.
You of course must have the CUDA version of pytorch installed (https://pytorch.org/)


Running on a single GPU
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The use of GPU acceleration has been streamlined in QMCTorch, the only modification
you need to do on your code is to specify `cuda=True` in the declaration of the wave function :


>>> # define the wave function
>>> wf = SlaterJastrow(mol, kinetic='jacobi',
>>>             configs='cas(2,2)',
>>>             cuda=True)

This will automatically port all the necesaary tensors to the GPU and offload all the corresponding operation
there.

Multi-GPU support
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The use of multiple GPUs is made possible through the `Horovod` library : https://github.com/horovod/horovod
A dedicated QMCTorch Solver has been developped to handle multiple GPU. To use this solver simply import it
and use is as the normal solver and only a few modifications are required to use horovod :


>>> import horovod.torch as hvd
>>> from deepqmc.solver.solver_slater_jastrow_horovod import SolverSlaterJastrowHorovod
>>>
>>> hvd.init()
>>> if torch.cuda.is_available():
>>>    torch.cuda.set_device(hvd.local_rank())
>>>
>>> # define the molecule
>>> mol = Molecule(atom='H 0 0 -0.69; H 0 0 0.69',
>>>                calculator='pyscf', basis='sto-3g',
>>>               unit='bohr', rank=hvd.local_rank())
>>>
>>> ...
>>>
>>> solver = SolverSlaterJastrowHorovod(wf=wf, sampler=sampler,
>>>                              optimizer=opt, scheduler=scheduler,
>>>                              rank=hvd.rank())
>>> ....

As you can see some classes need the rank of the process when they are defined. This is simply
to insure that only the master process generates the HDF5 files containing the information relative to the calculation.

The code can then be launched using the `horovodrun` executalbe :

::

    horovodrun -np 2 python <example>.py

See the horovod documentation for more details : https://github.com/horovod/horovod


This solver distribute the `Nw` walkers over the `Np` process . For example specifying 2000 walkers
and using 4 process will lead to each process using only 500 walkers. During the optimizaiton of the wavefunction
each process will compute the gradients of the variational parameter using their local 500 walkers.
The gradients are then averaged over all the processes before the optimization step takes place. This data parallel
model has been greatly succesfull in machine learning applications (http://jmlr.org/papers/volume20/18-789/18-789.pdf)Creating your own Backflow transformation
==============================================

We present here how to create your own backflow transformation. During the import you must import the base class of the backflow kernel

>>> from qmctorch.scf import Molecule
>>> from qmctorch.wavefunction import SlaterJastrowBackFlow
>>> from qmctorch.wavefunction.orbitals.backflow.kernels import BackFlowKernelBase


We can then use this base class to create a new backflow transformation kernel.
This is done in the same way one would create a new neural network layer in pytorch

>>> class MyBackflow(BackFlowKernelBase):
>>>
>>>     def __init__(self, mol, cuda, size=16):
>>>         super().__init__(mol, cuda)
>>>         self.fc1 = nn.Linear(1, size, bias=False)
>>>         self.fc2 = nn.Linear(size, 1, bias=False)
>>>
>>>     def forward(self, x):
>>>         original_shape = x.shape
>>>         x = x.reshape(-1,1)
>>>         x = self.fc2(self.fc1(x))
>>>         return x.reshape(*original_shape)

This backflow transformation consists of two fully connected layers. The calculation of the first and second derivative are then done via automatic differentiation
as implemented in the `BackFlowKernelBase` class. To use this new kernel in the `SlaterJastrowBackFlow` wave function ansatz we simply pass the class name as argument of the `backflow_kernel` keyword argument :

>>> # define the wave function
>>> wf = SlaterJastrowBackFlow(mol, kinetic='jacobi',
>>>                    backflow_kernel=MyBackflow,
>>>                    backflow_kernel_kwargs={'size' : 64},
>>>                    configs='single_double(2,2)')
Correlation and Blocking
========================================

One important part of the sampling is to estimate the correlation between the different sampling points.
To this end let's import the following modules


>>> from qmctorch.scf import Molecule
>>> from qmctorch.wavefunction import SlaterJastrow
>>> from qmctorch.sampler import Metropolis
>>> from qmctorch.utils import set_torch_double_precision
>>> from qmctorch.utils import blocking, plot_blocking_energy
>>> from qmctorch.utils import plot_correlation_coefficient, plot_integrated_autocorrelation_time


To obtain a bettter accuracy on our results we can switch to a double precision default
tensor type for pytorch :

>>> set_torch_double_precision()


Correlation coefficient
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The simplest way to estimate the decorelation time is to compute the autocorrelation coefficient of the local energy.
After having defined the molecule and the wave function as in the previous tutorial, let's define the sampler as follow:

>>> sampler = Metropolis(nwalkers=100, nstep=500, step_size=0.25,
>>>                     nelec=wf.nelec, ndim=wf.ndim,
>>>                     init=mol.domain('normal'),
>>>                     ntherm=0, ndecor=1,
>>>                     move={'type': 'all-elec', 'proba': 'normal'})

Compared to before we here record all the walker position along the trajectory. We can then define the solver as before
and run the following commands:

>>> pos = solver.sampler(solver.wf.pdf)
>>> obs = solver.sampling_traj(pos)
>>> rho, tau = plot_correlation_coefficient(obs.local_energy)

Which gives the following plot:

.. image:: ../../pics/correlation_coeff.png


On this picture is represented the autocorrelation coefficient of the local energy of all the walkers (transparent colorful line)
and the average of the autocorrelation coefficient (thick black line). This mean value is fitted with an exponential decay
to obtain the autocorrelation time that is here equal to 4.13 MCs.

Integrated autocorrelation time
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Another way to estimate the correlation time is to compute the integrated autocorrelation time. This can be done with

>>> plot_integrated_autocorrelation_time(obs.local_energy)

That leads to the following plot

.. image:: ../../pics/iac.png

A conservative estimate of the correlation time can be obtain when the iac cross the dashed line, leading here to a value of about 7 steps.


Energy blocking
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is also common practice to use blocking of the local energy values to reduce the variance. This can easily be done with

>>> eb = plot_blocking_energy(obs.local_energy, block_size=100,
>>>                           walkers='mean')

leading to the plot

.. image:: ../../pics/blocking.png

That shows the raw and blocked values of the mean local energy values.



Exploring the results with h5x
=====================================

The results and input of any the calculation performed by QMCTorch is stored in a dedicated HDF5 File
that can be explored using `h5x`. `h5x` is hosted at https://github.com/DeepRank/h5xplorer.
To install the browser, you clone and install the repository or use the PyPi package manager :

::

    pip install h5xplorer


To launch h5x simply execute the python file in the h5x folder

::

    cd QMCTroch/h5x
    python h5x.py

This will start the browser where you can load and explore data files generated byt QMCTorch.

SCF calculation
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The result of any SCF calculation generates an hdf5 file named by default

::

    <molecule_name>_<calculator>_<basis>.hdf5

Another name can be specified via the Molecule argument. This file contains all the data of the Molecule instance, in particular the basis set information, used in the calculation.
By browsing the file using h5x you can retreive any information needed. This file is also reused if possible to avoid computing the SCF of a previously studied molecular system.

.. image:: ../../pics/h5xmol_qmctorch.png



QMCTorch calculation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The results of QMCTorch calculation are also stored in a hdf5 file named

::

    <molecule_name>_<calculator>_<basis>_QMCTorch.hdf5

Another name can be specified via the Solver argument. This file contains all the information about the
wavefunction, the sampler the solver and the results of all the calculations performed in the script.
Context menus have been incorporate in order to plot the results directly from the hdf5 file.

.. image:: ../pics/h5x_solver.PNGCreating your own Jastrow Factor
====================================

We present here how to create your own electron-electron Jastrow factor and use it in QMCTorch

During the import you must import the base class of the electron-electron Jastrow.

>>> from qmctorch.scf import Molecule
>>> from qmctorch.wavefunction import SlaterJastrow,
>>> from qmctorch.wavefunction.jastrows.elec_elec.kernels import JastrowKernelElectronElectronBase


We can then use this base class to create a new Jastrow Factor. This is done in the same way one would create
a new neural network layer in pytorch.

>>> class MyJastrow(JastrowKernelElectronElectronBase)
>>>
>>>     def __init__(self, nup, ndown, cuda, size=16):
>>>         super().__init__(nup, ndown, cuda)
>>>         self.fc1 = nn.Linear(1, size, bias=False)
>>>         self.fc2 = nn.Linear(size, 1, bias=False)
>>>
>>>     def forward(self, x):
>>>         nbatch, npair = x.shape
>>>         x = x.reshape(-1,1)
>>>         x = self.fc2(self.fc1(x))
>>>         return x.reshape(nbatch, npair)

As seen above the prototype of the class constructor must be:

>>> def __init__(self, nup, ndown, cuda, **kwargs)

The list of keyword argument can contain any pairs such as ``size=16``.

This Jastrow use two fully connected layers. The size of the hidden layer is here controlled by a keyword argument ``size`` whose defauilt value is 16
It is important to note that the calculation of the first and second derivative of the jastrow kernel wrt the electronic positions are then done via automatic differentiation
as implemented in the `JastrowKernelElectronElectronBase` class. Hence there is no need to derive and implement these derivatives. However it
is necessary that the ``forward`` function, which takes as input a ``torch.tensor`` of
dimension ``[Nbatch, Npair]`` first reshape this tensor to ``[Nbatch*Npair,1]``, then applies the transformation on this tensor and finally reshape
the output tensor to ``[Nbatch, Npair]``.

To use this new Jastrow in the `SlaterJastrow` wave function ansatz we simply pass the class name as argument of the `jastrow_kernel` keyword argument. It is also
possible to specify the values of the keyword argument ``size`` with the ``jastrow_kernel_kwargs``. As seen below the pair of keyword argument and its value is passed as
a python dictionary :

>>> # define the wave function
>>> wf = SlaterJastrow(mol, kinetic='jacobi',
>>>                    jastrow_kernel=MyJastrow
>>>                    jastrow_kernel_kwargs={'size' : 64}
>>>                    configs='single_double(2,2)')

As previously, `mol` is an instance of the `qmctorch.scf.Molecule` class.Single point and Sampling trajectory
========================================

In this first tutorial we are going to define a simple molecule and sample its wave function.
Let's start by importing all the modules we need :

>>> from qmctorch.scf import Molecule
>>> from qmctorch.wavefunction import SlaterJastrow
>>> from qmctorch.sampler import Metropolis
>>> from qmctorch.utils import set_torch_double_precision
>>> from qmctorch.utils import plot_walkers_traj

To obtain a bettter accuracy on our results we can switch to a double precision default
tensor type for pytorch :

>>> set_torch_double_precision()


Gaussian orbitals : STO-3G with pyscf
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We now  need to define the molecule. We are here going to read the position of a water molecule from an xyz file.
The coordinates of an XYZ files are in Angstrom and we therefore specify this unit when defining the molecule.
We choose here to use pyscf with a samall sto-3g basis set to perform the SCF calculation on the molecule.

>>> # define the molecule
>>> mol = Molecule(atom='water.xyz', unit='angs',
>>>               calculator='pyscf', basis='sto-3g')

We can now define the QMCNet wavefunction. We choose here to include only the ground state of the of the system in the CI expansion.
The kinetic energy of the system will be computed using the Jacobi formula instead of the hessian of the wave function  wrt the electron positions through automatic differentiation.
We also include a Jastrow factor in the calculation of the wave function.

>>> # define the wave function
>>> wf = SlaterJastrow(mol, kinetic='jacobi',
>>>             configs='ground_state')

We then define the sampler as Metropolis Hasting samper, using only 100 walkers. Each walker is initially in a shpere around the molecule.
During the MC sampling each walker performs 500 steps of 0.25 atomic units. We here move only 1 electron per MC step per walkers using a normal distribution centred around their current position.

>>> # sampler
>>> sampler = Metropolis(nwalkers=100, nstep=500, step_size=0.25,
>>>                     nelec=wf.nelec, ndim=wf.ndim,
>>>                     init=mol.domain('normal'),
>>>                     move={'type': 'one-elec', 'proba': 'normal'})

We can now assemble the solver.

>>> # solver
>>> solver = SolverSlaterJastrow(wf=wf, sampler=sampler)

We can first perform  a single point calculation.
In this calculation the solver will sample the wave function and compute the energy of the system.

>>> obs = solver.single_point()

::

    Energy   :  -71.28987884521484  +/-  1.5550774335861206
    Variance :  241.82656860351562

This leads to an energy of about E = -71 hartrees witha  very large error bar of several hartrees. This bad result is due to the basis set
we are using here (sto-3g) that is a poor approximation of the slater orbitals.

To understand how the sampling works we can compute and visualize the sampling trajectory. To this end we are going to change the parameters of th sampler so that it keeps
the position of the walkers during the MC sampling. We will therefore keep the walker positions each 5 MC step starting from their initial positions.

>>> solver.sampler.ntherm = 0
>>> solver.sampler.ndecor = 5

We then have to sample the wave function and computes the energy for each sampling interval.

>>> pos = solver.sampler()
>>> obs = solver.sampling_traj(pos)

We can now visualize the results.

>>> plot_walkers_traj(obs.local_energy, walkers='mean')

.. image:: ../../pics/samplig_sto.png


Slater orbitals : DZP with ADF
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is possible to improve the results bu using ADF as a calculator and a DZP basis set. We can hence change the molecule to :

>>> # define the molecule
>>> mol = Molecule(atom='water.xyz', unit='angs',
>>>               calculator='adf', basis='dzp')

Without changing anything else in the calculation we now obtain a smaller error bar on the energy (few tenth of hartrees).

::

    Energy   :  -72.27135467529297  +/-  0.36148950457572937
    Variance :  13.06746768951416

and the following graph:

.. image:: ../../pics/samplig_dzp.png




Wave function Optimization
====================================

We present here a complete example on how to use QMCTorch on a H2 molecule.
We first need to import all the relevant modules :

>>> from torch import optim
>>> from qmctorch.scf import Molecule
>>> from qmctorch.wavefunction import SlaterJastrow,
>>> from qmctorch.solver import SolverSlaterJastrow
>>> from qmctorch.sampler import Metropolis
>>> from qmctorch.utils import set_torch_double_precision
>>> from qmctorch.utils import (plot_energy, plot_data)

To obtain a bettter accuracy on our results we can switch to a double precision default
tensor type for pytorch :

>>> set_torch_double_precision()

The first step is to define a molecule. We here use a H2 molecule with both hydrgen atoms
on the z-axis and separated by 1.38 atomic unit. We also choose to use ADF as SCF calculator using
a double zeta + polarization (dzp) basis set.

>>> # define the molecule
>>> mol = Molecule(atom='H 0 0 -0.69; H 0 0 0.69',
>>>                calculator='adf',
>>>                basis='dzp',
>>>                unit='bohr')

We then define the QMCNet wave function relative to this molecule. We also specify here
the determinants we want to use in the CI expansion. We use here a to include all the single
and double excitation with 2 electrons and 2 orbitals

>>> # define the wave function
>>> wf = SlaterJastrow(mol, kinetic='jacobi',
>>>              configs='single_double(2,2)')

As a sampler we use a simple Metropolis Hasting with 1000 walkers. The walkers are initially localized around the atoms.
Each walker will perform 2000 steps of size 0.2 atomic unit and will only keep the last position of each walker (`ntherm=-1`).
During each move all the the electrons are moved simultaneously within a normal distribution centered around their current location.

>>> # define the sampler
>>> sampler = Metropolis(nwalkers=1000,
>>>                      nstep=2000, step_size=0.2,
>>>                      ntherm=-1, ndecor=100,
>>>                      nelec=wf.nelec, init=mol.domain('atomic'),
>>>                      move={'type': 'all-elec', 'proba': 'normal'})


We will use the ADAM optimizer implemented in pytorch with custom learning rate for each layer.
We also define a linear scheduler that will decrease the learning rate after 100 steps

>>> # optimizer
>>> lr_dict = [{'params': wf.jastrow.parameters(), 'lr': 1E-3},
>>>            {'params': wf.ao.parameters(), 'lr': 1E-6},
>>>            {'params': wf.mo.parameters(), 'lr': 1E-3},
>>>            {'params': wf.fc.parameters(), 'lr': 1E-3}]
>>> opt = optim.Adam(lr_dict, lr=1E-3)
>>>
>>> # scheduler
>>> scheduler = optim.lr_scheduler.StepLR(opt, step_size=100, gamma=0.90)

We can now assemble the solver :

>>> # QMC solver
>>> solver = SolverSlaterJastrow(wf=wf, sampler=sampler, optimizer=opt, scheduler=None)

The solver needs to be configured. Many parameters of the optimization can be controlled as illustrated below

We can specify which observale to track during the optimization. Here only the local energies will be recorded
but one can also record the variational parameters
>>> solver.configure(track=['local_energy'])

Some variational parameters can be frozen and therefore not optimized. We here freeze the MO coefficients and the AO parameters
and therefore only the jastrow parametres and the CI coefficients will be optmized
>>> solver.configure(freeze=['ao', 'mo'])

Either the mean or the variance of local energies can be used as a loss function. We here use the mean
>> solver.configure(loss='energy')

The gradients can be computed directly via automatic differntiation or via a reduced noise formula
>>> solver.configure(grad='auto')

We also configure the resampling so that the positions of the walkers are updated by performing
25 MC steps from their previous positions after each optimization step.

>>> solver.configure(resampling={'mode': 'update',
>>>                             'resample_every': 1,
>>>                             'nstep_update': 25})




We can now run the optimization. We use here 250 optimization steps (epoch), using all the points
in a single mini-batch.

>>> data = solver.run(250)

Once the optimization is done we can plot the results

>>> plot_energy(solver.observable.local_energy, e0=-1.1645, show_variance=True)
>>> plot_data(solver.observable, obsname='jastrow.weight')

.. image:: ../../pics/h2_dzp.png
Creating a molecule
========================================

In this tutorial we present how to create a molecule and run the SCF calculation. First, the `Molecule` class must be imported :

>>> from qmctorch.scf import Molecule

This class can interface with `pyscf` and `ADF` to perform SCF calculations. Of course both software use different types of
atomic orbitals, respectively Gaussian type orbitals for `pyscf` and Slater type orbitals for `ADF`.


Geometry of the molecule
------------------------------------------

The geometry of the molecule can be specified through the `atom` keyword of the `Molecule` class. The units of the positions, `bohr` or `angs` (default is 'bohr')
can also be specified via the `unit` keyword argument.  The geometry can be passed as a single string

>>> Molecule(atom = 'H 0. 0. 0; H 0. 0. 1.', unit='bohr')

or via an XYZ file containing the geomtry of the molecular structure

>>> Molecule(atom='h2.xyz', unit='angs')

SCF calculations
--------------------------------------------

As mentionned above `QMCTorch` can use `pyscf` or `ADF` to perform SCF calculation on the molecular structure. At the moment only Hartree-Fock calculations
are supported but DFT calculations will be implemented later. We present here how to perform these SCF calculations


Gaussian orbitals with pyscf
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To use `pyscf` to compute the molecular orbitals of the system the `calculator` simply need to be set to `pyscf` as shown below.

>>> # define the molecule
>>> mol = Molecule(atom='H 0. 0. 0; H 0. 0. 1.', unit='angs',
>>>                calculator='pyscf', basis='sto-3g')

The `basis` keyword specify which basis set to use in the calculation. We use here a small `STO-3G` basis set. The exhaustive list of supported basis
set can be found here : https://pyscf.org/user/gto.html


Slater orbitals with ADF
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If a valid SCM license is found  QMCTorch can use `ADF` as a backend by simply switching the `calculator` to 'adf' as below :

>>> # define the molecule
>>> mol = Molecule(atom='H 0. 0. 0; H 0. 0. 1.', unit='angs',
>>>               calculator='adf', basis='dzp')

Here as well the ``basis`` keyword argument specifies the basis set used in the scf calculation.
The list of supported basis set can be found here : https://www.scm.com/doc/ADF/Input/Basis_sets_and_atomic_fragments.html

Additional basis sets, namely VB1, VB2, VB3, CVB1, CVB2 and CVB3, are available. These are STO valence and core-valence basis set presented by Ema et. al
in  "Polarized basis sets for Slater-type orbitals: H to Ne atoms", https://doi.org/10.1002/jcc.10227. Changing the ``basis``
 keyword argument to : ``basis=VB1``` will for examle use the small VB1 basis set during the SCF calculation.

Loading a SCF calculation
----------------------------------

By default QMCTorch will create a HDF5 file containing all the required information about the molecule and SCF calculation. The name of
this file is given by the name of the molecule, the calculator name and the basis set, e.g. `LiH_adf_dzp.hdf5` or 'water_pyscf_sto3g.xyz'. This files
can be loaded to instanciate the molecule object through the `load` keyword argument:

>>> mol = Molecule(load='LiH_adf_dzp.hdf5')

FermiNet Wave Function
----------------------------------------Samplers
--------------------------------

`QMCTorch` offers different sampler to propagate the walkers. The default sampler is a Metropolis-Hasting
that can be defined as follows :

>>> from qmctorch.sampler import Metropolis
>>> sampler = Metropolis(nwalkers=500,
>>>                  nstep=2000, step_size=0.2,
>>>                  ntherm=-1, ndecor=100,
>>>                  nelec=wf.nelec, init=mol.domain('atomic'),
>>>                  move={'type': 'all-elec', 'proba': 'normal'})

Optimizers
===========================

`QMCTorch` allows to use all the optimizers included in `pytorch` to opmtimize the QMCNet wave function.
The list of optimizers can be found here :  https://pytorch.org/docs/stable/optim.html 

For example to use the ADAM optimizer with different learning rate for each layer of the QMCNet wave funciton, one can simply define :

>>> from torch import optim
>>> lr_dict = [{'params': wf.jastrow.parameters(), 'lr': 1E-3},
>>>        {'params': wf.ao.parameters(), 'lr': 1E-6},
>>>        {'params': wf.mo.parameters(), 'lr': 1E-3},
>>>        {'params': wf.fc.parameters(), 'lr': 1E-3}]
>>> opt = optim.Adam(lr_dict, lr=1E-3)

Scheduler
==============================

Similarly QMCTorch allows to use scheduler to gradually decrease the learning rate during the optimization.
There as well all the scheduler of pytorch can be used : https://pytorch.org/docs/stable/optim.html
For example a simple scheudler that decrease the learning rate every number of epoch is simply defined as :

>>> from torch import optim
>>> scheduler = optim.lr_scheduler.StepLR(opt, step_size=100, gamma=0.90)



Neural Network Wave Function
================================

QMCTorch expresses the wave function of a molecular structure as a neural network. In most cases this network takes the electronic positions as input
and returns the corresponding value of the wave function. Several neural network wave function ansatz are implemented in QMCTorch :


.. toctree::
   :maxdepth: 1
   :caption: Wavefunction Ansatz:

   slaterjastrow
   slaterjastrow_backflow


Solvers
=========================

Solvers are responsibe to orchestrate the calculations by combining the different elements, Molecule, QMCNet wave function, samplers and optimizers/schedulers
The main solver is caled `SolverSlaterJastrow` and is defined as

>>> from qmctorch.solver import SolverSlaterJastrow
>>> solver = SolverSlaterJastrow(wf=wf, sampler=sampler,
>>>                       optimizer=opt, scheduler=scheduler, output='output.hdf5'

As soon as the solver its content is defined the HDF5 file specficied byt `out`. This file will contain all the parameter
of the solver and can be explored using the dedicated `h5x` browser. This solver allows to perform different calculations as detailled below

Single point calculation
----------------------------

A single point calculation sample the current wave function and computes the energy & variance of the system.
It can simply be done by:

>>> obs = solver.single_point()
>>> plt.hist(obs.local_energy)

`obs` is a `SimpleNamespace` instance with the following attributes:
  * `obs.pos` : position of the walkers
  * `obs.local_energy` : values of the local energies for each sampling point
  * `obs.energy` : energy of the systems (i.e. the mean of the local energy)
  * `obs.variance` : variance of the local energy values
  * `obs.error` : error on the energy

The result of the calculation will also be stored in the hdf5 output file. The energy distribution can be vizualised
with the matplotlib histogram function.

Sampling trajectory
----------------------------

It is possible to compute the local energy during the propagation of the wlakers to assess the convergence of the sampling
To this end the sampler must be configured to output the walker position after each `N` steps.
For example to start recording the walkers positions after 1000 MC steps and then record their position each 100 MC steps one can use :

>>> from qmctorch.utils import plot_walkers_traj
>>> solver.sampler.ntherm = 1000
>>> solver.sampler.ndecor = 100
>>> pos = solver.sampler(solver.wf.pdf)
>>> obs = solver.sampling_traj(pos)
>>> plot_walkers_traj(obs.local_energy)

There as well the results are returned in the `obs` SimpleNamespace and are stored in the hdf5 file.
The trajectory can be visualized with the `plot_wakers_traj` routine of QMCTorch

Wave function optimization
-------------------------------

Optimizing the wave function is the main task of the solver. Before otpimization starts the solver needs to be
configured properly.

>>> solver.configure(task='wf_opt', freeze=['ao', 'mo'])

To main task are available wave function optimization (`wf_opt`) and geometry optimization (`geo_opt`).
If a wave function optimization is selected the atom coordinate will be frozen while all the other parameters of the QMCNet will be optimized.
If a geometry optimization is selected only the atom coordinates will be optimized. One cal also freeze (i.e. not optimize) certain parameter groups.
In the example above the parameters of the atomic orbitals and molecular orbitals will be frozen,

One can specified the observale that needs to be recorded during the optimization.

>>> solver.track_observable(['local_energy'])

By default the local energy and all the variational parameters will be recorded.

As the system is optimized, one can resample the wave function by changing the positions of the walkers.
Several strategies are available to resample the wave function. The preferred one is to update the walkers by performing a small number of MC steps after each optimization step.
This can be specified with :

>>> solver.configure_resampling(mode='update', resample_every=1, nstep_update=25)

Finally we can now optimize the wave function using the `.run()` method of the solver.
This methods takes a few arguments, the number of optimization step, the batchsize, and some parameters to compute the gradients.

>>> data = solver.run(5, batchsize=None,
>>>                  loss='energy',
>>>                  grad='manual',
>>>                  clip_loss=False)

The results are returned in a SimpleNamespace and can be visualized with dedicated routines :

>>> plot_energy(solver.observable.local_energy, e0=-
>>>              1.1645, show_variance=True)

>>> plot_data(solver.observable, obsname='jastrow.weight')Slater Jastrow Backflow Wave Function
----------------------------------------

The Slater Jastrow Backflow wave function builds on the the Slater Jastrow wavefunction but adds a backflow transformation to
the electronic positions. Following this transformation, each electron becomes a quasi-particle whose position depends on all
electronic positions. The backflow transformation is given by :

.. math::

    q(x_i) = x_i + \sum_{j\neq i} \text{Kernel}(r_{ij}) (x_i-x_j)

The kernel of the transformation can be any function that depends on the distance between two electrons. A popular kernel
is simply the inverse function :

.. math::
    \text{Kernel}(r_{ij}) = \frac{\omega}{r_{ij}}

and is the default value in QMCTorch. However any other kernel function can be implemented and used in the code.

The wave function is then constructed as :

.. math::

    \Psi(R) = J(R) \sum_n c_n D_n^{\uparrow}(Q) D_n^{\downarrow}(Q)

The Jastrow factor is still computed using the original positions of the electrons while the determinant part uses the
backflow transformed positions.

Orbital Dependent Backflow Transformation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The backflow transformation can be different for each atomic orbitals.

.. math::

    q^\alpha(x_i) = x_i + \sum_{j\neq i} \text{Kernel}^\alpha(r_{ij}) (x_i-x_j)

where each orbital has its dedicated backflow kernel. This provides much more flexibility when optimizing the wave function.

Usage
^^^^^^^^^^^^^^^^^^^^^

This wave function can be used with

>>> from qmctorch.wavefunction import SlaterJastrowBackFlow
>>> from qmctorch.wavefunction.orbitals.backflow.kernels import BackFlowKernelInverse
>>> from qmctorch.wavefunction.jastrows.elec_elec.kernels.pade_jastrow_kernel import PadeJastrowKernel
>>>
>>> wf = SlaterJastrowBackFlow(self.mol, kinetic='jacobi',
>>>                            configs='single_double(2,2)',
>>>                            jastrow_kernel=PadeJastrowKernel,
>>>                            orbital_dependent_backflow=False,
>>>                            backflow_kernel=BackFlowKernelInverse)

Compared to the ``SlaterJastrow`` wave function, the kernel of the backflow transformation must be specified.
By default the inverse kernel will be used. Orbital dependent backflow orbitals can be easily achieved by using ``orbital_dependent_backflow=True``Slater Jastrow Wave Function
-----------------------------------


The ``SlaterJastrow`` neural-network wavefunction ansatz matches closely the wave function ansatz commonly used in QMC simulations. The wave function
is here expressed as

.. math::

    \Psi(R) = J(R) \sum_n c_n D_n^{\uparrow} D_n^{\downarrow}

The term `J(R)` is the so called Jastrow factor that captures the electronic correlation. The Jastrow factor is given by :

.. math::

    J(R) = \exp\left(  \sum_{i<j} \text{Kernel}(r_{ij}) \right)

where the sum runs over all the electron pairs and where the kernel function defines the action of the Jastrow factor. A common expression for the
kernel (and the default option for QMCTorch) is the Pade-Jastrow form given by:

.. math::

    \text{Kernel}(r_{ij}) = \frac{\omega_0 r_{ij}}{1+\omega r_{ij}}

where :math: `\omega_0` is a fixed coefficient equals to 0.25(0.5) for antiparallel(parallel) electron spins and :math: `\omega` a variational parameter.

The determinantal parts in the expression of :math: `\Psi` are given by the spin-up and spin-down slater determinants e.g. :

.. math::

    D_n^{\uparrow} = \frac{1}{\sqrt{N}} \begin{vmatrix} & & \\ & \phi_j(r_i) & \\ & & \end{vmatrix}


Implementation
^^^^^^^^^^^^^^^^^^^^^^^^

The Slater Jastrow function is implemented in QMCTorch as represented in the figure below :

.. image:: ../../pics/mol_nn.png

As seen on this figure, the input of the network are the positions of the walkers that should be sampled from the density of the wavefunction.
The first layer of the network computes the values of all the atomic orbitals at the position of each electron contained in a given walker configuration.
This new layer has been implemented in Pytorch and allows the optimizaiton of the atomic position (geometry optimization), and of the atomic basis parameters (basis exponents and coefficients).
The second layer transforms the atomic orbitals in molecular orbitals using a linear map whose weights are the molecular orbitals coefficients. The third layer, dubbed Slater Pooling,
computes the values of all the desired Slater determinants from the values of the molecular orbilals. The fina layer sums up the slater determinant using the CI coefficients as weight

In parallel a Jastrow layer has been implemented in Pytorch and allows the calculation of the Jastrow factor directly from the walkers positions.
The Jastrow factor is multiplied with the sum of the slater determinant to obtain the value of the wave function.

Usage
^^^^^^^^^^^^^^^^^^^^^^^
The ``SlaterJastrow`` wave function can instantiated following :

>>> wf = SlaterJastrow(mol, configs='single_double(2,2)', jastrow_kernel=PadeJastrowKernel)

The ``SlaterJastrow`` takes as first mandiatory argument a ``Molecule`` instance. The Slater determinants required in the calculation
are specified with the ``configs`` arguments which can take the following values :

  * ``configs='ground_state'`` : only the ground state SD
  * ``configs='cas(n,m)'`` : complete active space using n electron and m orbitals
  * ``configs='single(n,m)'`` : only single excitation using n electron and m orbitals
  * ``configs='single_double(n,m)'`` : only single/double excitation using n electron and m orbitals

Finally the kernel function of the Jastrow factor can be specifed using the ``jastrow_kernel``
The ``SlaterJastrow`` class accepts other initialisation arguments to fine tune some advanced settings. The default values
of these arguments are adequeate for most cases.

Orbital dependent Jastrow factor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Jastrow factor can be made orbital dependent with the ``SlaterOrbitalDependentJastrow``

>>> from qmctorch.wavefunction import SlaterOrbitalDependentJastrow
>>> from qmctorch.wavefunction.jastrows.elec_elec.kernels import PadeJastrowKernel
>>> wf = SlaterOrbitalDependentJastrow(mol, configs='single_double(2,4)'
>>>                                    jastrow_kernel=PadeJastrowKernel)
Pauli Net Wave Function
----------------------------------------qmctorch.utils.interpolate module
=================================

.. automodule:: qmctorch.utils.interpolate
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.jastrows.elec\_elec.kernels.pade\_jastrow\_kernel module
==============================================================================

.. automodule:: qmctorch.wavefunction.jastrows.elec_elec.kernels.pade_jastrow_kernel
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.utils package
======================

Submodules
----------

.. toctree::
   :maxdepth: 4

   qmctorch.utils.algebra_utils
   qmctorch.utils.hdf5_utils
   qmctorch.utils.interpolate
   qmctorch.utils.plot_data
   qmctorch.utils.stat_utils
   qmctorch.utils.torch_utils

Module contents
---------------

.. automodule:: qmctorch.utils
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.orbitals package
======================================

Subpackages
-----------

.. toctree::
   :maxdepth: 4

   qmctorch.wavefunction.orbitals.backflow

Submodules
----------

.. toctree::
   :maxdepth: 4

   qmctorch.wavefunction.orbitals.atomic_orbitals
   qmctorch.wavefunction.orbitals.atomic_orbitals_backflow
   qmctorch.wavefunction.orbitals.atomic_orbitals_orbital_dependent_backflow
   qmctorch.wavefunction.orbitals.norm_orbital
   qmctorch.wavefunction.orbitals.radial_functions
   qmctorch.wavefunction.orbitals.spherical_harmonics

Module contents
---------------

.. automodule:: qmctorch.wavefunction.orbitals
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction package
=============================

Subpackages
-----------

.. toctree::
   :maxdepth: 4

   qmctorch.wavefunction.jastrows
   qmctorch.wavefunction.orbitals
   qmctorch.wavefunction.pooling

Submodules
----------

.. toctree::
   :maxdepth: 4

   qmctorch.wavefunction.slater_combined_jastrow
   qmctorch.wavefunction.slater_jastrow
   qmctorch.wavefunction.slater_jastrow_backflow
   qmctorch.wavefunction.slater_jastrow_base
   qmctorch.wavefunction.slater_orbital_dependent_jastrow
   qmctorch.wavefunction.wf_base

Module contents
---------------

.. automodule:: qmctorch.wavefunction
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.jastrows.distance.electron\_nuclei\_distance module
=========================================================================

.. automodule:: qmctorch.wavefunction.jastrows.distance.electron_nuclei_distance
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.jastrows.elec\_elec\_nuclei.jastrow\_factor\_electron\_electron\_nuclei module
====================================================================================================

.. automodule:: qmctorch.wavefunction.jastrows.elec_elec_nuclei.jastrow_factor_electron_electron_nuclei
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.jastrows.elec\_nuclei.jastrow\_factor\_electron\_nuclei module
====================================================================================

.. automodule:: qmctorch.wavefunction.jastrows.elec_nuclei.jastrow_factor_electron_nuclei
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.orbitals.backflow.kernels.backflow\_kernel\_fully\_connected module
=========================================================================================

.. automodule:: qmctorch.wavefunction.orbitals.backflow.kernels.backflow_kernel_fully_connected
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.jastrows.elec\_elec.orbital\_dependent\_jastrow\_kernel module
====================================================================================

.. automodule:: qmctorch.wavefunction.jastrows.elec_elec.orbital_dependent_jastrow_kernel
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.utils.torch\_utils module
==================================

.. automodule:: qmctorch.utils.torch_utils
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.pooling.slater\_pooling module
====================================================

.. automodule:: qmctorch.wavefunction.pooling.slater_pooling
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.jastrows.elec\_elec.kernels package
=========================================================

Submodules
----------

.. toctree::
   :maxdepth: 4

   qmctorch.wavefunction.jastrows.elec_elec.kernels.fully_connected_jastrow_kernel
   qmctorch.wavefunction.jastrows.elec_elec.kernels.jastrow_kernel_electron_electron_base
   qmctorch.wavefunction.jastrows.elec_elec.kernels.pade_jastrow_kernel
   qmctorch.wavefunction.jastrows.elec_elec.kernels.pade_jastrow_polynomial_kernel

Module contents
---------------

.. automodule:: qmctorch.wavefunction.jastrows.elec_elec.kernels
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.jastrows.distance.scaling module
======================================================

.. automodule:: qmctorch.wavefunction.jastrows.distance.scaling
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.jastrows.distance.electron\_electron\_distance module
===========================================================================

.. automodule:: qmctorch.wavefunction.jastrows.distance.electron_electron_distance
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.utils.stat\_utils module
=================================

.. automodule:: qmctorch.utils.stat_utils
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.jastrows.elec\_elec.jastrow\_factor\_electron\_electron module
====================================================================================

.. automodule:: qmctorch.wavefunction.jastrows.elec_elec.jastrow_factor_electron_electron
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.orbitals.backflow.orbital\_dependent\_backflow\_kernel module
===================================================================================

.. automodule:: qmctorch.wavefunction.orbitals.backflow.orbital_dependent_backflow_kernel
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.scf.calculator.pyscf module
====================================

.. automodule:: qmctorch.scf.calculator.pyscf
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.jastrows.elec\_elec.kernels.fully\_connected\_jastrow\_kernel module
==========================================================================================

.. automodule:: qmctorch.wavefunction.jastrows.elec_elec.kernels.fully_connected_jastrow_kernel
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.slater\_jastrow\_backflow module
======================================================

.. automodule:: qmctorch.wavefunction.slater_jastrow_backflow
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.utils.plot\_data module
================================

.. automodule:: qmctorch.utils.plot_data
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.scf.calculator package
===============================

Submodules
----------

.. toctree::
   :maxdepth: 4

   qmctorch.scf.calculator.adf
   qmctorch.scf.calculator.calculator_base
   qmctorch.scf.calculator.pyscf

Module contents
---------------

.. automodule:: qmctorch.scf.calculator
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.jastrows package
======================================

Subpackages
-----------

.. toctree::
   :maxdepth: 4

   qmctorch.wavefunction.jastrows.distance
   qmctorch.wavefunction.jastrows.elec_elec
   qmctorch.wavefunction.jastrows.elec_elec_nuclei
   qmctorch.wavefunction.jastrows.elec_nuclei

Submodules
----------

.. toctree::
   :maxdepth: 4

   qmctorch.wavefunction.jastrows.jastrow_factor_combined_terms

Module contents
---------------

.. automodule:: qmctorch.wavefunction.jastrows
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.slater\_orbital\_dependent\_jastrow module
================================================================

.. automodule:: qmctorch.wavefunction.slater_orbital_dependent_jastrow
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.orbitals.spherical\_harmonics module
==========================================================

.. automodule:: qmctorch.wavefunction.orbitals.spherical_harmonics
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.jastrows.elec\_nuclei.kernels.pade\_jastrow\_kernel module
================================================================================

.. automodule:: qmctorch.wavefunction.jastrows.elec_nuclei.kernels.pade_jastrow_kernel
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.utils.hdf5\_utils module
=================================

.. automodule:: qmctorch.utils.hdf5_utils
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.jastrows.jastrow\_factor\_combined\_terms module
======================================================================

.. automodule:: qmctorch.wavefunction.jastrows.jastrow_factor_combined_terms
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.scf.molecule module
============================

.. automodule:: qmctorch.scf.molecule
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.orbitals.backflow.kernels.backflow\_kernel\_power\_sum module
===================================================================================

.. automodule:: qmctorch.wavefunction.orbitals.backflow.kernels.backflow_kernel_power_sum
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.orbitals.norm\_orbital module
===================================================

.. automodule:: qmctorch.wavefunction.orbitals.norm_orbital
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.pooling.orbital\_projector module
=======================================================

.. automodule:: qmctorch.wavefunction.pooling.orbital_projector
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.orbitals.radial\_functions module
=======================================================

.. automodule:: qmctorch.wavefunction.orbitals.radial_functions
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.slater\_jastrow module
============================================

.. automodule:: qmctorch.wavefunction.slater_jastrow
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.sampler.hamiltonian module
===================================

.. automodule:: qmctorch.sampler.hamiltonian
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.orbitals.backflow.kernels.backflow\_kernel\_inverse module
================================================================================

.. automodule:: qmctorch.wavefunction.orbitals.backflow.kernels.backflow_kernel_inverse
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.jastrows.elec\_elec\_nuclei.kernels package
=================================================================

Submodules
----------

.. toctree::
   :maxdepth: 4

   qmctorch.wavefunction.jastrows.elec_elec_nuclei.kernels.fully_connected_jastrow_kernel
   qmctorch.wavefunction.jastrows.elec_elec_nuclei.kernels.jastrow_kernel_electron_electron_nuclei_base

Module contents
---------------

.. automodule:: qmctorch.wavefunction.jastrows.elec_elec_nuclei.kernels
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.sampler package
========================

Submodules
----------

.. toctree::
   :maxdepth: 4

   qmctorch.sampler.generalized_metropolis
   qmctorch.sampler.hamiltonian
   qmctorch.sampler.metropolis
   qmctorch.sampler.sampler_base
   qmctorch.sampler.walkers

Module contents
---------------

.. automodule:: qmctorch.sampler
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.slater\_combined\_jastrow module
======================================================

.. automodule:: qmctorch.wavefunction.slater_combined_jastrow
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.scf.calculator.adf module
==================================

.. automodule:: qmctorch.scf.calculator.adf
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.jastrows.elec\_elec\_nuclei.kernels.fully\_connected\_jastrow\_kernel module
==================================================================================================

.. automodule:: qmctorch.wavefunction.jastrows.elec_elec_nuclei.kernels.fully_connected_jastrow_kernel
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.orbitals.atomic\_orbitals module
======================================================

.. automodule:: qmctorch.wavefunction.orbitals.atomic_orbitals
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.orbitals.backflow.kernels package
=======================================================

Submodules
----------

.. toctree::
   :maxdepth: 4

   qmctorch.wavefunction.orbitals.backflow.kernels.backflow_kernel_autodiff_inverse
   qmctorch.wavefunction.orbitals.backflow.kernels.backflow_kernel_base
   qmctorch.wavefunction.orbitals.backflow.kernels.backflow_kernel_fully_connected
   qmctorch.wavefunction.orbitals.backflow.kernels.backflow_kernel_inverse
   qmctorch.wavefunction.orbitals.backflow.kernels.backflow_kernel_power_sum
   qmctorch.wavefunction.orbitals.backflow.kernels.backflow_kernel_square

Module contents
---------------

.. automodule:: qmctorch.wavefunction.orbitals.backflow.kernels
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.orbitals.backflow.kernels.backflow\_kernel\_autodiff\_inverse module
==========================================================================================

.. automodule:: qmctorch.wavefunction.orbitals.backflow.kernels.backflow_kernel_autodiff_inverse
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.solver.solver\_base module
===================================

.. automodule:: qmctorch.solver.solver_base
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.jastrows.elec\_elec.kernels.jastrow\_kernel\_electron\_electron\_base module
==================================================================================================

.. automodule:: qmctorch.wavefunction.jastrows.elec_elec.kernels.jastrow_kernel_electron_electron_base
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.slater\_jastrow\_base module
==================================================

.. automodule:: qmctorch.wavefunction.slater_jastrow_base
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.jastrows.elec\_elec package
=================================================

Subpackages
-----------

.. toctree::
   :maxdepth: 4

   qmctorch.wavefunction.jastrows.elec_elec.kernels

Submodules
----------

.. toctree::
   :maxdepth: 4

   qmctorch.wavefunction.jastrows.elec_elec.jastrow_factor_electron_electron
   qmctorch.wavefunction.jastrows.elec_elec.orbital_dependent_jastrow_kernel

Module contents
---------------

.. automodule:: qmctorch.wavefunction.jastrows.elec_elec
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.pooling package
=====================================

Submodules
----------

.. toctree::
   :maxdepth: 4

   qmctorch.wavefunction.pooling.orbital_configurations
   qmctorch.wavefunction.pooling.orbital_projector
   qmctorch.wavefunction.pooling.slater_pooling

Module contents
---------------

.. automodule:: qmctorch.wavefunction.pooling
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.wf\_base module
=====================================

.. automodule:: qmctorch.wavefunction.wf_base
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.sampler.generalized\_metropolis module
===============================================

.. automodule:: qmctorch.sampler.generalized_metropolis
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.solver.solver\_slater\_jastrow module
==============================================

.. automodule:: qmctorch.solver.solver_slater_jastrow
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.jastrows.elec\_nuclei package
===================================================

Subpackages
-----------

.. toctree::
   :maxdepth: 4

   qmctorch.wavefunction.jastrows.elec_nuclei.kernels

Submodules
----------

.. toctree::
   :maxdepth: 4

   qmctorch.wavefunction.jastrows.elec_nuclei.jastrow_factor_electron_nuclei

Module contents
---------------

.. automodule:: qmctorch.wavefunction.jastrows.elec_nuclei
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.orbitals.backflow.kernels.backflow\_kernel\_square module
===============================================================================

.. automodule:: qmctorch.wavefunction.orbitals.backflow.kernels.backflow_kernel_square
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.jastrows.elec\_nuclei.kernels package
===========================================================

Submodules
----------

.. toctree::
   :maxdepth: 4

   qmctorch.wavefunction.jastrows.elec_nuclei.kernels.fully_connected_jastrow_kernel
   qmctorch.wavefunction.jastrows.elec_nuclei.kernels.jastrow_kernel_electron_nuclei_base
   qmctorch.wavefunction.jastrows.elec_nuclei.kernels.pade_jastrow_kernel

Module contents
---------------

.. automodule:: qmctorch.wavefunction.jastrows.elec_nuclei.kernels
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.jastrows.elec\_nuclei.kernels.fully\_connected\_jastrow\_kernel module
============================================================================================

.. automodule:: qmctorch.wavefunction.jastrows.elec_nuclei.kernels.fully_connected_jastrow_kernel
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.orbitals.backflow.kernels.backflow\_kernel\_base module
=============================================================================

.. automodule:: qmctorch.wavefunction.orbitals.backflow.kernels.backflow_kernel_base
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.orbitals.backflow.backflow\_transformation module
=======================================================================

.. automodule:: qmctorch.wavefunction.orbitals.backflow.backflow_transformation
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.orbitals.atomic\_orbitals\_orbital\_dependent\_backflow module
====================================================================================

.. automodule:: qmctorch.wavefunction.orbitals.atomic_orbitals_orbital_dependent_backflow
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.jastrows.elec\_nuclei.kernels.jastrow\_kernel\_electron\_nuclei\_base module
==================================================================================================

.. automodule:: qmctorch.wavefunction.jastrows.elec_nuclei.kernels.jastrow_kernel_electron_nuclei_base
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.jastrows.elec\_elec\_nuclei.kernels.jastrow\_kernel\_electron\_electron\_nuclei\_base module
==================================================================================================================

.. automodule:: qmctorch.wavefunction.jastrows.elec_elec_nuclei.kernels.jastrow_kernel_electron_electron_nuclei_base
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.pooling.orbital\_configurations module
============================================================

.. automodule:: qmctorch.wavefunction.pooling.orbital_configurations
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.scf.calculator.calculator\_base module
===============================================

.. automodule:: qmctorch.scf.calculator.calculator_base
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.sampler.walkers module
===============================

.. automodule:: qmctorch.sampler.walkers
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.solver package
=======================

Submodules
----------

.. toctree::
   :maxdepth: 4

   qmctorch.solver.solver_base
   qmctorch.solver.solver_slater_jastrow
   qmctorch.solver.solver_slater_jastrow_horovod

Module contents
---------------

.. automodule:: qmctorch.solver
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.orbitals.backflow.orbital\_dependent\_backflow\_transformation module
===========================================================================================

.. automodule:: qmctorch.wavefunction.orbitals.backflow.orbital_dependent_backflow_transformation
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.jastrows.elec\_elec\_nuclei package
=========================================================

Subpackages
-----------

.. toctree::
   :maxdepth: 4

   qmctorch.wavefunction.jastrows.elec_elec_nuclei.kernels

Submodules
----------

.. toctree::
   :maxdepth: 4

   qmctorch.wavefunction.jastrows.elec_elec_nuclei.jastrow_factor_electron_electron_nuclei

Module contents
---------------

.. automodule:: qmctorch.wavefunction.jastrows.elec_elec_nuclei
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.orbitals.backflow package
===============================================

Subpackages
-----------

.. toctree::
   :maxdepth: 4

   qmctorch.wavefunction.orbitals.backflow.kernels

Submodules
----------

.. toctree::
   :maxdepth: 4

   qmctorch.wavefunction.orbitals.backflow.backflow_transformation
   qmctorch.wavefunction.orbitals.backflow.orbital_dependent_backflow_kernel
   qmctorch.wavefunction.orbitals.backflow.orbital_dependent_backflow_transformation

Module contents
---------------

.. automodule:: qmctorch.wavefunction.orbitals.backflow
   :members:
   :undoc-members:
   :show-inheritance:

Subpackages
-----------

.. toctree::
   :maxdepth: 4

   qmctorch.sampler
   qmctorch.scf
   qmctorch.solver
   qmctorch.utils
   qmctorch.wavefunction

Module contents
---------------

.. automodule:: qmctorch
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.sampler.sampler\_base module
=====================================

.. automodule:: qmctorch.sampler.sampler_base
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.orbitals.atomic\_orbitals\_backflow module
================================================================

.. automodule:: qmctorch.wavefunction.orbitals.atomic_orbitals_backflow
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.utils.algebra\_utils module
====================================

.. automodule:: qmctorch.utils.algebra_utils
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.scf package
====================

Subpackages
-----------

.. toctree::
   :maxdepth: 4

   qmctorch.scf.calculator

Submodules
----------

.. toctree::
   :maxdepth: 4

   qmctorch.scf.molecule

Module contents
---------------

.. automodule:: qmctorch.scf
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.jastrows.elec\_elec.kernels.pade\_jastrow\_polynomial\_kernel module
==========================================================================================

.. automodule:: qmctorch.wavefunction.jastrows.elec_elec.kernels.pade_jastrow_polynomial_kernel
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.sampler.metropolis module
==================================

.. automodule:: qmctorch.sampler.metropolis
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.solver.solver\_slater\_jastrow\_horovod module
=======================================================

.. automodule:: qmctorch.solver.solver_slater_jastrow_horovod
   :members:
   :undoc-members:
   :show-inheritance:
qmctorch.wavefunction.jastrows.distance package
===============================================

Submodules
----------

.. toctree::
   :maxdepth: 4

   qmctorch.wavefunction.jastrows.distance.electron_electron_distance
   qmctorch.wavefunction.jastrows.distance.electron_nuclei_distance
   qmctorch.wavefunction.jastrows.distance.scaling

Module contents
---------------

.. automodule:: qmctorch.wavefunction.jastrows.distance
   :members:
   :undoc-members:
   :show-inheritance:
