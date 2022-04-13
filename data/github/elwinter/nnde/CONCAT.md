# Contributing to nnde

If you have modified the nnde package for your own work, and wish to contribute to the project, we welcome your contributions of code and documentation. To contribute, fork the repository on GitHub, and submit a pull request at https://github.com/elwinter/nnde.

# Contact

Submit bug reports and feature requests on the [repository issue tracker](https://github.com/elwinter/nnde/issues).

Eric Winter <eric.winter62@gmail.com>
# nnde

`nnde` is a package of Python modules which implement a set of neural networks which solve ordinary and partial differential equations using technique of [Lagaris et al. (1998)](https://dx.doi.org/10.1109/72.712178).

The `nnde` package provides a pure-Python implementation of one of the earliest approaches to using neural networks to solve differential equations - the trial function method [Lagaris et al. (1998)](https://dx.doi.org/10.1109/72.712178). It was initially developed primarily as a vehicle for understanding the internal workings of feedforward neural networks, without the restrictions imposed by an existing neural network framework. The target audiences are machine learning researchers who want to develop an understanding of the basic steps needed to solve differential equations with neural networks and researchers from general engineering fields who want to solve a differential equation without using a complex neural network library.

The package is described in [Winter and Weigel, 2021](https://github.com/elwinter/blob/master/paper.pdf). Additional details on the method are given in [Winter 2020](https://github.com/elwinter/proposal/blob/master/proposal.pdf).

Demos are available in the repository [`nnde_demos`](https://github.com/elwinter/nnde_demos).

# Install and Use

```bash
pip install nnde
git clone https://github.com/elwinter/nnde_demos
pip install matplotlib # optional
cd nnde_demos; python lagaris01_demo.py
```

# Developer

```bash
git clone https://github.com/elwinter/nnde
cd nnde
pip install -e .
python setup.py test  # Deprecated
pytest  # Preferred
```

# An overview of nnde

The `nnde` package was developed as a tool to help understand how neural networks can be used to solve differential equations. The effort was originally inspired by work done in the 1990s showing how feedforward neural networks could be used as universal approximators. That work led to efforts to show how neural networks could be applied to differential equations. A good example of this technique can be found in Lagaris et al (1993).

The basic technique is straightforward:

* Arrange the differential equation in a standard form:

  G(x, y, dy/dx, ...) = 0

* Define a trial solution `y_trial` that can be substituted into the differential equation. This trial solution includes a component computed by a neural network, as well as a term which incorporates all boundary conditions.

* Use the analytical definition of `y_trial` to determine the analytical forms of the various derivatives of `y` used in the differential equation.

* Define the structure of the neural network so that it has one input per independent variable, and a single output.

* Using a set of training points defined on the domain of the differential equation, run the network and use the output to compute the value of the trial solution and its derivatives.

* Compute the value of the standardized differential equation `G()` at each training point.

* Compute the loss function as the sum of the squared values of `G()`.

* Train a neural network to solve the equation for the trial solution by minimizing the loss function at each training point: `L = SUM(G_i^2)`. This is done by adjusting the network parameters (weights and biases) until a satifactory solution is obtained.

The `nnde` package is divided into 3 major packages (`differentialequation`, `neuralnetwork`, and `trialfunction`), and two auxiliary packages (`exceptions` and `math`).

The `differentialequation` package is divided into two sub-packages: `ode` (for ordinary differential equations) and `pde` (for partial differential equations). These sub-packages are very lightweight - they are composed of abstract classes used to define the methods that the user-defined differential equation  must implement. These are primarily methods to evaluate the differential equation itself, and the various derivatives required for its evaluation. Each of these functions depends on the training points, and computes the values and derivatives of the trial solution during evaluation The ode and pde packages provide classes for 1st-order ODE initial-value problems, 1st-order PDE initial value problems, and diffusion problems in 1, 2, and 3 spatial dimensions. These classes can be used as-is to solve the corresponding equation types by defining the required methods.

The `neuralnetwork` package provides the core of the `nnde` functionality. Each module provides a customizable neural network tailored to the needs of a specific problem type. Currently support types are the same as those supported in the `differentialequation` package. A different class was used for each equation type because the details of the computation of the network output differ slightly for each equation type.

The `trialfunction` package provides a previously unavailable capability: determine the structure of the boundary condition component of the trial solution automatically. The algorithm for determining the form of this component of the trial solution is recursive, and the number of terms grows rapidly as the dimensionality of the problem is increased. However, if the user provides the boundary conditions (as another set of functions), these modules can automatically construct the boundary condition term, greatly easing the problem definition burden on the user. The modules also provide the option to short-circuit this process by allowing the user to define an optimized form of the boundary condition function that can greatly reduce the required amount of computation. Trial function classes are provided for diffusion problems in 1, 2, and 3 spatial dimensions. The code which solves 1st-order ODE IVP uses a trivial form of the boundary condition component (a constant), and therefore its use is coded directly, without a separate class.

A Jupyter notebook ([`tutorial.ipynb`](https://github.com/elwinter/nnde/tree/master/examples)) providing a structured walkthrough of the process of defining and solving a problem using the `nnde` package is available in the `examples` directory.

# Contribute

If you discover bugs in the `nnde` package, please create an issue at the project repository on GitHub at https://github.com/elwinter/nnde.

If you find the nnde package useful, we welcome your contributions of code and documentation. To contribute, fork the repository on GitHub, and submit a pull request at https://github.com/elwinter/nnde.

# Contact

Submit bug reports and feature requests on the [repository issue tracker](https://github.com/elwinter/nnde/issues).

Eric Winter <eric.winter62@gmail.com>
---
title: '`nnde`: A Python package for solving differential equations using neural networks'
tags:
  - Python
  - neural networks
  - differential equations
authors:
  - name: Eric Winter
    orcid: 0000-0001-5226-2107
    affiliation: 1
  - name: R.S. Weigel
    orcid: 0000-0002-9521-5228
    affiliation: 1
affiliations:
  - name: Department of Physics and Astronomy, George Mason University, Fairfax, Virginia, USA
    index: 1
date: 14 February 2022
bibliography: paper.bib
---

# Summary

Neural networks have been shown to have the ability to solve differential equations [@Yadav:2015; @Chakraverty:2017]. `nnde` is a pure-Python package for the solution of ordinary and partial differential equations of up to second order. We present the results of sample runs showing the effectiveness of the software in solving the two-dimensional diffusion problem.

# Statement of need

The `nnde` package provides a pure-Python implementation of one of the earliest approaches to using neural networks to solve differential equations - the trial function method [@Lagaris:1998]. The `nnde` package was initially developed as a vehicle for understanding the internal workings of feedforward neural networks, without the constraints imposed by an existing neural network framework. It has since been enhanced to provide the capability to solve differential equations of scientific interest, such as the diffusion equation described here. The ultimate goal of the package is to provide the capability to solve systems of coupled partial differential equations, such as the equations of magnetohydrodynamics.

Development of the `nnde` package began before the widespread adoption of modern neural network software frameworks. In the Python ecosystem, the most popular packages are TensorFlow (`https://tensorflow.org`) and PyTorch (`https://pytorch.org`). These frameworks are designed to be application-neutral - they can be used to develop neural networks with arbitrary architectures for arbitrary learning objectives. The primary advantages of these frameworks are autodifferentiation and distributed computing. By recording the sequence of mathematical operations performed in the forward pass through the network, autodifferentiation can automatically compute the gradients of the loss function with respect to each of the network parameters, as well as the network inputs. The latter capability is central to solving differential equations. Autodifferentiation also greatly reduces the volume of code that must be developed to solve a given problem. The distributed computing capability allows a network to take advantage of GPU-enabled hardware, and multiple compute nodes, to speed the calculation, with little or no code changes required. The `nnde` package uses a more direct method - precomputed derivative functions for the components of the differential equations of interest and the trial solution. This code is typically faster than TensorFlow or PyTorch, but requires more hand-crafted code to solve a given problem.

The most commonly used methods for solving differential equations are the Finite Element Method (FEM) and Finite Difference Method (FDM). However, these methods can be difficult to parallelize due to the need for communication between computational elements at the boundaries of the allocated subgrids. These models can also  have large storage requirements for model outputs. The neural network method is straightforward to parallelize due to the independent characteristics of the computational nodes in each network layer. Additionally, the trained network solution is more compact than an FDM or FEM solution because storage of only the network weights and biases is required. The neural network solution is mesh-free and does not require interpolation to retrieve the solution at a non-grid point, as is the case with FDM or FEM. Once the network is trained, computing a solution at any spatial or temporal scale requires only a series of matrix multiplications, one per network layer. The trained solution is a sum of arbitrary differentiable basis functions, and therefore the trained solution is also differentiable, which is particularly useful when computing derived quantities such as gradients and fluxes. This approach has led to several different classes of methods for solving ODEs and PDEs with neural networks. The recent surge in interest in "physics-informed neural networks" [@Raissi:2019] is an indication of the dynamic nature of the field.

# Description

`nnde` implements a version of the trial function algorithm described by @Lagaris:1998. This software also incorporates a modification of the trial function algorithm to automatically incorporate arbitrary Dirichlet boundary conditions of the problem directly into the neural network solution.

As a concrete example of the sort of problem that can be solved using `nnde`, consider the diffusion equation in two dimensions:

$$\frac {\partial \psi} {\partial t} - D \left( \frac {\partial^2 \psi} {\partial x^2} + \frac {\partial^2 \psi} {\partial y^2} \right) = 0$$

With all boundaries fixed at $0$ and with an initial condition of

$$\psi(x,y,0) = \sin(\pi x) \sin(\pi y)$$

the analytical solution is

$$\psi_a(x,y,t) = e^{-2\pi^2 D t} \sin(\pi x) \sin(\pi y)$$

The `nnde` package was used to create a neural network with a single hidden layer and 10 hidden nodes and trained to solve this problem. The error in the trained solution for the case of $D=0.1$ is shown as a function of time in \autoref{fig:diff2d_error}.

![Difference between the trained neural network solution $\psi_t(x,y,t)$ and the analytical solution $\psi_a(x,y,t)$ of the diffusion problem in 2 spatial dimensions using `nnde` with 10 nodes.\label{fig:diff2d_error}](figures/Y_e.png)

# Software repository

The `nnde` software is available at https://github.com/elwinter/nnde.

A collection of example python scripts using `nnde`  is available at https://github.com/elwinter/nnde_demos.

A collection of example Jupyter notebooks using `nnde` is available at https://github.com/elwinter/nnde_notebooks.

# References
