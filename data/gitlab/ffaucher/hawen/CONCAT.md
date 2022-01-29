# **hawen: time-HArmonic Wave modEling and INversion using Hybridizable Discontinuous Galerkin Discretization.**

___

## version 1.0.1

___

> ### Dedicated website with more infos and tutorials: [https://ffaucher.gitlab.io/hawen-website/](https://ffaucher.gitlab.io/hawen-website/).
>
> ### The software solves time-harmonic wave problems for acoustic and elastic media using the Hybridizable Discontinuous Galerkin (HDG) method for discretization. It combines mpi and OpenMP parallelism to solve large-scale applications such as Earth's imaging and helioseismology.
> ### It can handle the forward problem (propagation of waves) as well as the inverse problem (parameter identification).
> ### For comments or additional information, contact <florian.faucher@univie.ac.at>.



___

## License
> ### The code is distributed with no warranty, under the GNU General Public License v3.0, see the LICENSE file.

## Folder structure
> Once the repository is cloned/downloaded, the following items must be in the folder:
> - **CHANGELOG** 
> - **code**
> - **doc**
> - **examples**
> - **LICENSE**
> - **README.md**

## Installation
### Mandatory Dependencies
> The following dependencies are necessary before the compilation of hawen:
>   - **make**, **MPI** libraries and **Fortran90** compiler with mpi.
>   - **BLAS**, **LAPACK** and **scaLAPACK** libraries (to link with MUMPS). It is also possible to use **Openblas** instead of the two first.
>
>   - [**MUMPS**](http://mumps.enseeiht.fr/) to solve the linear system.
>   - [**Metis**](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) is used for mesh partition.

> ___

>> **NOTE FOR LINUX USERS**: All of the packages can be retrieved from the official repositories, on `Ubuntu 20.04` these are obtained with
>>
>> 	sudo apt-get install  build-essential
>> 	sudo apt-get install  libopenmpi-dev
>> 	sudo apt-get install  libblas-dev libscalapack-openmpi-dev liblapack-dev
>> 	sudo apt-get install  libmetis-dev
>> 	sudo apt-get install  libmumps-dev
>>
>> **REMARK**: as of July 2020, the latest version of MUMPS is 5.3.3 while the one deployed on the linux repository is the 5.2.1. We recommend to install the latest version to benefit from the latest features.

> ___

### Optional Dependencies
> We refer to the [User's documentation](https://ffaucher.gitlab.io/hawen-website) for the information on the optional dependencies, which include
>  - [**arb**](http://arblib.org/) toolbox for the efficient computation of special functions.
>  - [**ARPACK/pARPACK**](https://www.caam.rice.edu/software/ARPACK/)  for the computation of eigenvalues and eigenvectors.


### Edit file **config/make.version_config**
> Once the dependencies are installed, you can move the the `code` repository: 
>
>		cd code
>
> Before to run `make`, one needs to create the file **make.version_config** in the config/ folder. Templates are given in the folder **config/config-files_make/**, for different architectures. This file contains the path towards the dependencies, that is, where to find the LAPACK, scaLAPACK, MUMPS and metis libraries. Note that any dependencies that has been added to MUMPS must be added in hawen (e.g., scotch partitioner, etc.).
> The variable `ARCH` is used to indicate the architecture: it must be one of the following: `GNU`,`INTEL` or `PGI`, respectively for GNU, Intel or PGI based compilation. All have been successfully tested on supercomputers.
>
>> - For instance, the template file **config/config-files_make/make.version_config__GNU.default** can be used directly if one has installed the dependencies from the linux repository as given above.
>>
>>			cp config/config-files_make/make.version_config__GNU.default config/make.version_config
>>
> The main keywords are described below for the `config/config-files_make/make.version_config` file.
> Other keywords are for the non-mandatory libraries and can remain blank if these are not linked.
>		
>		# F90 indicates the compiler, e.g., mpif90 or mpiifort.
>		F90  := mpif90
>		# ARCH indicates the architecture, it can either be
>		# GNU, INTEL or PGI
>		ARCH := GNU
>			
>		# metis is a mandatory dependency, both the library
>		# and include folder must be linked.
>		METISDIR  = /usr/lib
>		LMETIS    = -L$(METISDIR) -lmetis
>		IMETIS    = -I/usr/include/
>		
>		# MUMPS is a mandatory dependency, both the library
>		# -lcmumps -lzmumps -lmumps_common and -lpord must
>		# be linked, as well as the include folder
>		MUMPSDIR= /usr/lib
>		LMUMPS  = -L$(MUMPSDIR) -lcmumps -lzmumps
>		LMUMPS += -L$(MUMPSDIR) -lmumps_common -lpord
>		IMUMPS  = -I/usr/include/
>		
>		# additional dependencies for Mumps for linear algebra
>		# are required with LAPACK and SCALAPACK
>		LMUMPS += -L/usr/lib -lblas -lscalapack-openmpi -llapack
>


### (Optional) Edit file **config/make.version**
>
> The file  **config/make.version** contains information on the compilation, such as the choice of optimization or debug flags, as well as the choice of method to be compiled. The choice of problem is made prior to the compilation, such that one executable is linked with one problem (forward or inverse, acoustic or elastic medium, etc). By default, the forward problem associated with acoustic media is compiled.
>
> 		vi config/make.version
> 
> You can choose the propagator that is compiled and the problem, by default, we have the acoustic isotropic forward problem:
>		# PROBLEM indicates the type of problem: this can be `forward` or `inversion`.
>		PROBLEM        := forward
>
>		# The different choices are propagator are currently:
>		#   - `helmholtz_acoustic-iso`   for the acoustic isotropic case.
>		#   - `helmholtz_elastic-iso`    for the elastic isotropic case.
>		#   - `helio_ODE-radial`         for the 1D helioseismology ODE under spherical symmetry.
>		#   - `helio_scalar-conjugated`  for the helioseismology scalar wave problem. 
>		PROPAGATOR     := helmholtz_acoustic-iso
>
>		# MUMPS_LATEST can be set to 1 if you have compiled MUMPS version 5.3.3 or higher.
>		MUMPS_LATEST=0
>		# DEPENDENCY_ARB can be set to 1 if you have linked arb library.
>		DEPENDENCY_ARB=0
>		# DEPENDENCY_ARPACK can be set to 1 if you have linked ARPACK library.
>		DEPENDENCY_ARPACK=0
>		# DEPENDENCY_PARPACK can be set to 1 if you have linked PARPACK library.
>		DEPENDENCY_PARPACK=0
>
> We refer to the [User's documentation](https://ffaucher.gitlab.io/hawen-website) for more information on the corresponding equations and options.
> Note that the dependencies are **not** activated by default, and one must change the value to **1** is  it is linked with the compilation, it also requires that the proper libraries are linked via the `config/make.version_config` file.
>

### (Optional) Edit file **src/m_define_precision.f90**
> The file  **src/m_define_precision.f90** contains the choice of arithmetric precision, which has to be selected prior to the compilation (for instance, simple or double precision matrix operations). It is detailed in the [User's documentation](https://ffaucher.gitlab.io/hawen-website)

### Compilation

> ___
>> **WARNING**: in order to compile several methods (`helmholtz_acoustic-iso` or `helmholtz_elastic-iso`, `forward` or `inversion`), you always have to run
>>
>>		make clean
>> 
>> **in between sucessive compilation.**

> ___

> Once you have
> 1. Create the file `config/make.version_config` to link the mandatory libraries,
> 2. Possibly adjusted the compilation in the file `config/make.version`
> 3. Possibly adjusted the arithmetic precision in the file `src/m_define_precision.f90`
>
> you can compile with (still in the `code` repository):
>
>		make
>
> Upon successful operation, the generated executable is created in the folder `bin` at root directory:
>
>		ls ../bin
>
> To compile the different methods you have to:
> 1. compile using `make`
> 2. run `make clean`
> 3. choose the method in the `config/make.version` file (e.g. elastic, acoustic, forward, inversion)
> 4. go to step 1 and repeat steps 1--4 for all methods you need.


## Utilization
> Once the executable(s) have been successfully compiled, you can try it in the folder `examples`. Assuming you are in the `code` folder:
>
>		cd ../examples
>
> Note that the code works with an **input parameter files** that provides all information. These are detailed in the dedicated section of the [User's documentation](https://ffaucher.gitlab.io/hawen-website).
>
> A benchmark to run the inverse problem can also be downloaded from the [hawen website](https://ffaucher.gitlab.io/hawen-website), in the section **Tutorial**.
---
title: "`hawen`: time-harmonic wave modeling and inversion using hybridizable discontinuous Galerkin discretization"
tags:
  - Fortran
  - wave equations
  - inverse problems
  - geophysics
  - helioseismology
  - Hybridizable Discontinuous Galerkin method
  - viscoelasticity
  - MPI and OpenMP parallelism 
authors:
  - name: Florian Faucher
    orcid: 0000-0003-4958-7511
    affiliation: "1"
affiliations:
 - name: Faculty of Mathematics, University of Vienna, Oskar-Morgenstern-Platz 1, A-1090 Vienna, Austria.
   index: 1
date: September 2020
bibliography: paper.bib

---

# Summary

Many applications such as seismic and medical imaging, material sciences, 
or helioseismology and planetary science, aim to reconstruct properties 
of a non directly accessible or non-visible interior.
For this purpose, they rely on waves whose propagation through a medium 
interrelates with the physical properties (density, sound speed, etc.) of 
this medium.
The methodology for imaging with waves comprises of two main stages 
illustrated in \autoref{fig:setup}.
In the data acquisition stage (\autoref{fig:setup}a), the 
medium response to probing waves is recorded (e.g., seismic waves 
from Earthquakes recorded by ground network).
In the second stage, we rely on a reconstruction procedure which 
iteratively updates an initial model of physical parameters, 
so that numerical simulations approach the measurements (\autoref{fig:setup}b).
This procedure is employed, for instance, for seismic (reconstruction of 
subsurface layers) and medical (disease diagnostic) imaging,
see the introduction of @Faucher2020adjoint and the references therein.


![Illustration of inverse wave problems. 
\textbf{a)} Acquisition stage: probing waves are sent though the medium, 
and the medium's response is recorded by devices positioned on a portion 
of the domain. 
\textbf{b)} The reconstruction algorithm starts from initial properties 
and compares simulations of wave propagation with the measurements, then
iteratively updates those properties. The green block corresponds to
the forward problem (modeling the propagation of waves) and the orange 
ones to the iterative correction. `hawen` solves both the forward and 
inverse problems associated with time-harmonic waves.
  \label{fig:setup}](figures/setup_global.pdf){width=13cm} 


# Statement of need 

`hawen` is designed to address large-scale inverse wave problems, and includes
\begin{enumerate}
  \item simulating time-harmonic waves in heterogeneous 
        media (e.g., visco-acoustic and visco-elastic 
        propagation, modal equations in helioseismology); 
  \item performing iterative minimization to solve the inverse 
        problem, cf. \autoref{fig:setup}.
\end{enumerate}
It combines MPI and OpenMP parallelism, it is deployed on
supercomputers and has been used in studies on seismic imaging 
[@Faucher2020FRgWIacoustic; @Faucher2020adjoint; @Faucher2020DAS], 
as well as in helioseismology [@Pham2020Siam; @Pham2019Esaim; @Pham2019RR; @Pham2020RRscalar].
The software works with an input \emph{parameter file}: 
a text file listing the problem configuration (dimension, choice of viscous 
model, etc.), such that the computational details are transparent for the users.
`hawen` has its dedicated 
\href{https://ffaucher.gitlab.io/hawen-website/download/hawen_users-guide_documentation.pdf}{User's guide} [@HawenUserGuide] that details
its utilization and all the available functionalities, it is available 
on the software dedicated \href{https://ffaucher.gitlab.io/hawen-website/}{website}, 
which contains illustrations, installation guide and tutorials.


One specificity of `hawen` is to implement the hybridizable 
discontinuous Galerkin (HDG) method [@Arnold2002; @Cockburn2009] 
for the discretization of the wave equations.
It helps reduce the cost of the computations by providing smaller 
linear systems compared to, e.g., finite elements, depending on the 
order of approximation [@Kirby2012; @Faucher2020adjoint].
For seismic applications, current software mostly relies on the
spectral element method [@Komatitsch1998; @Komatitsch1999], such
as \href{https://geodynamics.org/cig/software/specfem3d/}{specfem},
or on finite differences discretization.

The software handles unstructured domains (allowing input meshes 
of different formats) and thus accounts for complex geometry. 
In addition, `hawen` relies on the massively parallel sparse direct 
solver MUMPS [@Amestoy2001; @Amestoy2006] for the matrix factorization,
hence handles problems with many right-hand sides.
Despite the technicality of the underlying methods, the purpose of `hawen` 
is to provide a unified and evolutive framework to address large-scale 
visco-elastic and visco-acoustic inverse problems, with a user-friendly 
interface using the input parameter files.


# Modeling the propagation of waves

The first feature of `hawen` is to simulate the propagation of 
time-harmonic waves in different types of medium: this is the
forward problem in \autoref{fig:setup}. 
The propagation of waves is characterized by a Partial Differential 
Equation which depends on the type of the waves (e.g., mechanical or 
electromagnetic), and on the type of medium considered (e.g., fluid or solid)
[@Carcione2007; @Slawinski2010; @Faucher2017].
`hawen` handles visco-acoustic and visco-elastic propagation,
as well as waves propagation for helioseismology, these are 
described in the \href{https://ffaucher.gitlab.io/hawen-website/download/hawen_users-guide_documentation.pdf}{User's documentation} [@HawenUserGuide].
Our implementation is based upon the HDG method and the computational 
steps are illustrated in \autoref{fig:forward}.

![Illustration of the computational steps for the 
  numerical resolution of the forward problem. 
  For simplicity, we illustrate with an homogeneous 
  medium.
\label{fig:forward}](figures/forward.pdf){width=13cm} 


# Inverse problem via iterative minimization

In the forward problem, the objective is to 
simulate the propagation of waves, given the 
medium physical properties.
Conversely, the objective of the inverse 
problem is to recover the physical properties, 
given some observations of waves, see \autoref{fig:setup}.
To address the inverse problem, `hawen` solves the 
minimization problem:
\begin{equation}\label{eq:misfit_generic}
  \min_{\boldsymbol{m}} \mathcal{J}(\boldsymbol{m}) \, \qquad \text{with} \quad
  \mathcal{J}(\boldsymbol{m}) \,=\, \mathrm{dist}\big(\mathcal{F}(\boldsymbol{m}), \, \boldsymbol{d}\big) \, .
\end{equation}
The misfit function $\mathcal{J}$ is defined to evaluate 
a distance function comparing the observed data $\boldsymbol{d}$ 
with simulations of wave propagation restricted to the position 
of the measurements: $\mathcal{F}(\boldsymbol{m})$, where $\boldsymbol{m}$ represents the 
physical properties. 
The iterative minimization scheme successively updates the 
physical properties used for the simulations (\autoref{fig:setup}b),
and follows a Newton-type method [@Nocedal2006]. In the context of 
seismic imaging, it is the Full Waveform Inversion method, cf. @Virieux2009.

`hawen` offers several options to conduct the 
iterative minimization, such as the choice of misfit function
and method to conduct the minimization.
These are further listed and detailed in 
the \href{https://ffaucher.gitlab.io/hawen-website/download/hawen_users-guide_documentation.pdf}
{software documentation} [@HawenUserGuide].

# Acknowledgements

The research of FF is supported by the Austrian Science 
Fund (FWF) under the Lise Meitner fellowship M 2791-N.
We acknowledge the use of the 
\href{https://vsc.ac.at/home/}{Vienna Scientific Clusters vsc3 and vsc4}
for the numerical applications, as well as with the TGCC cluster via 
the GENCI resource allocation project AP010411013.

# References
# **hawen: time-HArmonic Wave modEling and INversion using Hybridizable Discontinuous Galerkin Discretization.**

___

> ### Dedicated website with more infos and tutorials: [https://ffaucher.gitlab.io/hawen-website/](https://ffaucher.gitlab.io/hawen-website/).
>
> 
> ### This code is used to solve time-harmonic wave problem using Hybridizable Discontinuous Galerkin (HDG) discretization. It can be used for the forward problem (modeling) and the inverse one (quantitative parameter identification).
> ### For comments or additional information, contact <florian.faucher@univie.ac.at>.

___

## Examples folder
>
> This folder contains an acoustic and elastic modeling benchmarks to illustrate the code. It assumes that the executables for forward problem in acoustic and elastic are compiled, respectively `helmholtz_acoustic-iso` and `helmholtz_elastic-iso`. As indicated in the README of the root directory, do not forget to use `make clean` between the compilations.
> Start by untaring the benchmark:
>
>		tar -Jxf benchmark.tar.xz
>
> Navigate to the folder
>
>		cd benchmark
>
> There are two parameter files in the folder: 
> - `par.modeling_acoustic` for the acoustic case.
> - `par.modeling_elastic` for the elastic case.
>
> You can open the files to see the parameters (e.g., the values of the wave speeds, the frequency solved, etc). For more information on the parameter files, we refer to the [User's documentation](https://ffaucher.gitlab.io/hawen-website), which is also in `doc` folder.
>

## Launching the code
>
> You can launch the code with or without parallelism, we refer to ${HAWEN_BIN_PATH} for the path where the executables have been compiled (which should be `../../bin`)
> The acoustic case launches with
>
>		mpirun -np 4 -x OMP_NUM_THREADS=2 ${HAWEN_BIN_PATH}/forward_helmholtz_acoustic-iso_hdg.out parameter=par.modeling_acoustic
>
> using 4 mpi and 2 threads per mpi. The user can select the combination depending on its architecture.
>
> The elastic test-case runs with
>
>		mpirun -np 4 -x OMP_NUM_THREADS=1 ${HAWEN_BIN_PATH}/forward_helmholtz_elastic-iso_hdg.out parameter=par.modeling_elastic
>
> using 4 mpi and 1 thread per mpi. 
>
> Upon running, several information are printed on screen, which follows the computational operations.
>

## Visualization of the results
>
> After the computations are done, the `results' folder contains the time-harmonic wave solutions. One can visualize these solutions using [ParaView](https://www.paraview.org/). For instance, for the elastic test-case, one can see the solution at 0.4Hz with
>
>		paraview results_elastic/wavefield/wavefield_frequency_0.00000E+00_4.00000E-01Hz_src000001.pvtu
>
> It is recommended to adjust the scale (colorbar) to better see the wave. Note that all solutions (velocities and stress components, real and imaginary parts) are encompassed in one file.

___

## Inversion benchmark
>
> A benchmark to run the inverse problem can also be downloaded from the [hawen website](https://ffaucher.gitlab.io/hawen-website), in the section **Tutorial**. 

___

## this folder contains the documentation for hawen:

### user's guide pdf

### configuration file to generate doxygen documentation:
### it requires doxygen version at least 1.8.18
doxygen doxygen-config-file.txt

creates the doxygen folder which contains the html documentation and 
the latex to gnerate a pdf file as well


