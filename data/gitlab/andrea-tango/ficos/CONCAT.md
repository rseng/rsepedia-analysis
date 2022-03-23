# FiCoS: a fine-grained and coarse-grained GPU-powered deterministic simulator for biochemical networks

FiCoS (**Fi**ne- and **Co**arse-grained **S**imulator) is a novel efficient **deterministic simulator** of biological systems based on two different integration methods belonging to the Runge-Kutta family: the Dormand–Prince (DOPRI) method [[1](#1), [2](#2)], used in the absence of stiffness, and the Radau IIA method [[3](#3), [4](#4)] exploited when the system is stiff.\
Specifically, FiCoS exploits the **DOPRI5** and **RADAU5** methods, which are explicit and implicit **adaptive Runge-Kutta methods of order 5**, capable of varying the integration step-size during the resolution of the **Ordinary Differential Equation (ODE)** systems.

FiCoS has been designed and developed to simulate **Reaction-Based Models (RBMs)**, which are widely used to propose a mechanistic description of the biological process, allowing for achieving a detailed comprehension of the underlying behaviour of the cellular system under analysis.

RBMs are generally employed within different computational tasks, such as **Parameter Estimation**, **Sensitivity Analysis**, and **Parameter Sweep Analysis**. Unfortunately, these tasks require a huge amount of simulations of the same RBM with distinct parameterizations. Moreover,  RBMs can be easily characterised by hundreds or thousands of molecular species. In both cases, the required simulation(s) of the dynamics can be extremely demanding, rapidly overtaking the capability of modern **Central Processing Units (CPUs)**.

To deal with these problems, FiCoS fully exploits the parallelism provided by modern **Graphics Processing Units (GPUs)**, by using a **combined fine- and coarse-grained parallelization strategy**. In particular, the fine-grained parallelization is applied to distribute the calculations required to solve each ODE over multiple GPU cores, while the coarse-grained parallelization is exploited to perform several independent simulations in a parallel fashion.

  1. [How to cite FiCoS](#cite) 
  2. [Compilation and input parameters](#compInp)
  3. [Models](#models)
  4. [Jupyter Notebooks](#notebooks)
  5. [License](#lic)
  6. [Contacts](#cont)
  7. [References](#ref)

## <a name="cite"></a>How to cite FiCoS ##

A detailed description of FiCoS, as well as a complete experimental comparison against the state-of-the-art simulators (i.e., LSODA [[5](#5)], VODE [[6](#6)], cupSODA [[7](#7)], and LASSIE [[8](#8)]), by using the models described in [Models](#Models), can be found in:

- Tangherloni A., Nobile M.S., Cazzaniga P., Capitoli G., Spolaor S., Rundo R., Mauri G., and Besozzi D.: _FiCoS: a fine-grained and coarse-grained GPU-powered deterministic simulator for biochemical networks_, <a href="https://www.biorxiv.org/content/10.1101/2021.01.15.426855v3">biorxiv</a>, 2021.01.15.426855, 2021. doi: 10.1101/2021.01.15.426855.

## <a name="compInp"></a>Compilation and input parameters ##

We provide the source code of a **single-core CPU version** of FiCoS as well as the source code of the **GPU-accelared version** of FiCoS.
The CPU version has been written in vanilla `C++`, while the GPU version rely on the `NVIDIA CUDA library` (**version <=9.2**). 
In order to compile both the CPU and GPU versions of FiCoS, we provide a bash script (compile.sh) in the folders <a href="https://gitlab.com/andrea-tango/ficos/-/tree/master/CPU">CPU</a> and <a href="https://gitlab.com/andrea-tango/ficos/-/tree/master/GPU">GPU</a>, respectively.

The resulting executable files (named FiCoS_gpu and FiCoS_cpu) are standalone executable programs, which can be easily executed by providing the following parameters:

- path to the input folder;
- output folder.

Optional parameters can be provided:
- `-v1`: to enable the verbose modality (basic information, model summary, and integrator parameters);
- `-v2`: to enable a more detailed verbose modality (the previous information and the ODE system);
- `-v3`: to enable the most detailed verbose modality (the previous information, the Jacobian of ODE system, the integrator steps);
- `-f`: to enable the calculation of the fitness function, described in [[9](#9), [10](#10)], for the Parameter Estimation.

The followiong two parameters can also be provided to the GPU version of FiCoS:
- `-h`: to enable a heuristic to select the number of threads and blocks (default: 32 threads);
- the GPU number (default: 0).

The input folder must contain the input files describing the biochemical systems under investigation.
The specification of the input files can be read <a href="https://gitlab.com/andrea-tango/ficos/-/blob/master/InputFiles.pdf">here</a>

The CPU version of FiCoS can also be compiled on any supported architecture (e.g., GNU/Linux, Microsoft Windows, macOS) using the following compilation command (note that FiCoS **requires** the **LAPACK** (Linear Algebra PACKage) library to be compiled and executed):
```
g++ -std=c++11 src/ficos.cpp src/simulation.cpp src/dopri.cpp src/radau.cpp -O3 -o  FiCoS_cpu -llapack
```

The GPU version of FiCoS can also be compiled on any supported architecture (e.g., GNU/Linux, Microsoft Windows, macOS) using the following compilation command:
```
nvcc src/ficos.cu src/simulation.cu -O3 -Xptxas -O3 -std=c++11 -use_fast_math -rdc=true -o FiCoS_gpu -lcublas -lcublas_device -lcudadevrt -gencode arch=compute_61,code=compute_61
```
The command above creates a binary executable file runnable on GPUs with at least a compute capability equal to 6.1.

## <a name="models"></a>Models ##

The folder <a href="https://gitlab.com/andrea-tango/ficos/-/tree/master/Models">Models</a> contains all the RBMs used to collect the results showed in the paper. In particular, we uploaded the synthetic models (both symmetric and asymmetric), and the real models:
- Autophagy/Translation switch based on the mutual inhibition of MTORC1 and ULK1 (<a href="https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0116550">paper</a>);
- Ras/cAMP/PKA pathway in the yeast _Saccharomyces cerevisiae_ (<a href="https://www.sciencedirect.com/science/article/abs/pii/S0168165607016380">paper</a>);
- A human intracellular core metabolic pathway in a red blood cell (<a href="https://link.springer.com/chapter/10.1007/978-3-030-34585-3_17">paper</a>).

## <a name="notebooks"></a>Jupyter Notebooks ##

The folder <a href="https://gitlab.com/andrea-tango/ficos/-/tree/master/Notebooks">Notebooks</a> contains a Jupyter Notebook showing how FiCoS can be used to perform a 2D Parameter Sweep Analysis, a Jupyter Notebook showing how FiCoS can be used to perform a Sensitivity Analysis, and a Jupyter Notebook showing how FiCoS can be combined with a meta-heuristic to perform a Parameter Estimation.

## <a name="lic"></a>License ##

FiCoS is licensed under the terms of the GPL-3 license.

## <a name="cont"></a>Contacts ##

For questions or support, please contact Andrea Tangherloni (<andrea.tangherloni@unibg.it>).

## <a name="ref"></a>References ##
[<a name="1"></a>1] Dormand J.R. and Prince P.J.: _A family of embedded Runge-Kutta formulae_, J. Comput. Appl. Math., 1980.\
[<a name="2"></a>2] Hairer E., Nørsett S.P., and Wanner G.: _Solving ordinary differential equations I_, Springer-Verlag, 2008.\
[<a name="3"></a>3] Hairer E. and Wanner G.: _Stiff differential equations solved by Radau methods_, J. Comput. Appl. Math. 1999.\
[<a name="4"></a>4] Hairer E. and Wanner G.: _Solving ordinary differential equations II_, Springer-Verlag, 2002.\
[<a name="5"></a>5] Petzold L.: _Automatic selection of methods for solving stiff and nonstiff systemsof ordinary differential equations_, SIAM J. Sci. Stat. Comp., 1983.\
[<a name="6"></a>6] Brown P.N., Byrne G.D., and Hindmarsh A.C.: _VODE: A variable-coefficient ODE solver_, SIAM J. Sci. Stat. Comp., 1989.\
[<a name="7"></a>7] Nobile M.S., Cazzaniga P., Besozzi D., and Mauri G.: _GPU-accelerated simulations of mass-action kinetics models with cupSODA_, J. Supercomput., 2014.\
[<a name="8"></a>8] Tangherloni A., Nobile M.S., Besozzi D., Mauri G., and Cazzaniga P.: _LASSIE: simulating large-scale models of biochemical systems on GPUs_, BMC Bioinform., 2017.\
[<a name="9"></a>9] Nobile M.S., Tangherloni A., Rundo L., Spolaor S., Besozzi D., Mauri G., and Cazzaniga P.: _Computational intelligence for parameter estimation of biochemical systems_, Proc. of the IEEE Congress on Evolutionary Computation (CEC), 2018.\
[<a name="10"></a>10] Tangherloni A., Spolaor S., Cazzaniga P., Besozzi D., Rundo L., Mauri G., and Nobile M.S.: _Biochemical parameter estimation vs. benchmark functions: a comparative study of optimization performance and representation design_, Appl. Soft. Comput., 2019.
