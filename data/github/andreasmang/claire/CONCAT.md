# CLAIRE

[![status](https://img.shields.io/badge/arXiv-1808.04487-red)](https://arxiv.org/abs/1808.04487)
[![status](https://joss.theoj.org/papers/d4af0a15946fd2161935be018820243b/status.svg)](https://joss.theoj.org/papers/d4af0a15946fd2161935be018820243b)
[![status](https://img.shields.io/badge/licencse-GPL-blue)](https://github.com/andreasmang/claire/tree/master/LICENSE)

* Are you looking for **examples**? Check the [doc/examples](https://github.com/andreasmang/claire/tree/master/doc/examples) folder.
* Are in interested in **how CLAIRE works**? Check the [documentation](#clairedoc).
* Are you interested in **what CLAIRE is**? Read the [about](#claireabout) section.


## Documentation <a name="clairedoc"></a>
* [News](doc/README-NEWS.md): Recent updates and news are listed in [doc/README-NEWS.md](doc/README-NEWS.md)
* [Installation](doc/README-INSTALL.md): An installation guide can be found in [doc/README-INSTALL.md](doc/README-INSTALL.md)
* [Examples](doc/README-RUNME.md): A description of how to execute and use CLAIRE on your system can be found in [doc/README-RUNME.md](doc/README-RUNME.md).
* [Contributing](doc/CONTRIBUTING.md): If you have a bug report, a feature request, or you are interested in contributing to CLAIRE visit [doc/CONTRIBUTING.md](doc/CONTRIBUTING.md).
* [References](doc/README-REFERENCES.md): If you use CLAIRE as part of your research, please refer to [doc/README-REFERENCES.md](doc/README-REFERENCES.md) for information on citing our work.

The links above point to individual `markdown` files. These files can be found in the [doc](https://github.com/andreasmang/claire/tree/master/doc) subfolder. Basic examples for how to execute CLAIRE can be found in the [doc/examples](https://github.com/andreasmang/claire/tree/master/doc/examples) folder. The NIREP dataset used to test CLAIRE can be downloaded [here](https://github.com/andreasmang/nirep).


## About <a name="claireabout"></a>
**CLAIRE** stands for *Constrained Large Deformation Diffeomorphic Image Registration*. It is a C/C++ software package for velocity-based diffeomorphic image registration in three dimensions. Its performance is optimized for multi-core CPU systems (`cpu` branch) and multi-node, multi-GPU architectures (`gpu` branch; default). The CPU version uses MPI for data parallelism, and has been demonstrated to scale on several supercomputing platforms. CLAIRE can be executed on large-scale state-of-the-art computing systems as well as on local compute systems with limited resources.

Notice that the CPU version is accurate and running but new features are currently only being added to the GPU version. The GPU code is a major revision and therefore considered the default and recommended for use. 


<p align="center">
<img src="doc/figs/claire4brains.jpg" alt="CLAIRE4Brains"  width="800"/>
</p>

If there are any issues, you have questions, you would like to give us feedback or you have feature requests, do not hesitate to send an email to <andreas@math.uh.edu>.

If you plan on using CLAIRE in your research please cite the following manuscript:
A. Mang, A. Gholami, C. Davatzikos & G. Biros. *CLAIRE: A distributed-memory solver for constrained large deformation diffeomorphic image registration*. SIAM Journal on Scientific Computing 41(5):C548--C584, 2019 [[arxiv](https://arxiv.org/abs/1808.04487), [sisc](https://epubs.siam.org/doi/abs/10.1137/18M1207818)]. Additional references are listed [here](doc/README-REFERENCES.md).


## Contributors
George Biros, Malte Brunn, Amir Gholami, James Herring, Naveen Himthani, [Andreas Mang](mailto:andreas@math.uh.edu), and Miriam Mehl.

If you want to contribute to CLAIRE, read the guidlines (see [doc/CONTRIBUTING.md](doc/CONTRIBUTING.md)). 


## Code of Conduct
See [doc/CODE_OF_CONDUCT.md](doc/CODE_OF_CONDUCT.md).


## License
Read the [LICENSE](https://github.com/andreasmang/claire/tree/master/LICENSE) file for more details.

# Libmorton
[![Build Status](https://travis-ci.org/Forceflow/libmorton.svg?branch=master)](https://travis-ci.org/Forceflow/libmorton) [![license](https://img.shields.io/github/license/mashape/apistatus.svg)](https://opensource.org/licenses/MIT)

 * Libmorton is a **C++ header-only library** with methods to efficiently encode/decode between 64, 32 and 16-bit morton codes and coordinates, in 2D and 3D. *Morton order* is also known as *Z-order* or *[the Z-order curve](https://en.wikipedia.org/wiki/Z-order_curve)*.
 * Libmorton is a **lightweight and portable** library - in its most basic form it only depends on standard C++ headers. Architecture-specific optimizations are implemented incrementally.
 * This library is under active development. SHIT WILL BREAK.
 * More info and some benchmarks in these blogposts: [*Morton encoding*](http://www.forceflow.be/2013/10/07/morton-encodingdecoding-through-bit-interleaving-implementations/), [*Libmorton*](http://www.forceflow.be/2016/01/18/libmorton-a-library-for-morton-order-encoding-decoding/) and [*BMI2 instruction set*](http://www.forceflow.be/2016/11/25/using-the-bmi2-instruction-set-to-encode-decode-morton-codes/)

## Usage
Just include *libmorton/morton.h*. This will always have functions that point to the most efficient way to encode/decode Morton codes. If you want to test out alternative (and possibly slower) methods, you can find them in *libmorton/morton2D.h* and *libmorton/morton3D.h*.

<pre>
// ENCODING 2D / 3D morton codes, of length 32 and 64 bits
inline uint_fast32_t morton2D_32_encode(const uint_fast16_t x, const uint_fast16_t y);
inline uint_fast64_t morton2D_64_encode(const uint_fast32_t x, const uint_fast32_t y);
inline uint_fast32_t morton3D_32_encode(const uint_fast16_t x, const uint_fast16_t y, const uint_fast16_t z);
inline uint_fast64_t morton3D_64_encode(const uint_fast32_t x, const uint_fast32_t y, const uint_fast32_t z);
// DECODING 2D / 3D morton codes, of length 32 and 64 bits
inline void morton2D_32_decode(const uint_fast32_t morton, uint_fast16_t& x, uint_fast16_t& y);
inline void morton2D_64_decode(const uint_fast64_t morton, uint_fast32_t& x, uint_fast32_t& y);
inline void morton3D_32_decode(const uint_fast32_t morton, uint_fast16_t& x, uint_fast16_t& y, uint_fast16_t& z);
inline void morton3D_64_decode(const uint_fast64_t morton, uint_fast32_t& x, uint_fast32_t& y, uint_fast32_t& z);
</pre>

If you want to take advantage of the BMI2 instruction set (only available on Intel Haswell processors and newer), make sure `__BMI2__` is defined before you include `morton.h`.

## Testing
The *test* folder contains tools I use to test correctness and performance of the libmorton implementation. This section is under heavy re-writing, but might contain some useful code for advanced usage.

## Thanks
 * To [@gnzlbg](https://github.com/gnzlbg) and his Rust implementation [bitwise](https://github.com/gnzlbg) for finding bugs in the Magicbits code 
 * Everyone making comments and suggestions on the [original blogpost](http://www.forceflow.be/2013/10/07/morton-encodingdecoding-through-bit-interleaving-implementations/)

## TODO
 * Write better test suite (with L1/L2 trashing, better tests, ...)
 * BMI2 implementation for 2D methods
 * A better naming system for the functions, because m3D_e_sLUT_shifted? That escalated quickly.
# CLAIRE: References

Go back to [README.md](../README.md).

If you plan on using CLAIRE in your research please cite the following two manuscripts:

* A. Mang, A. Gholami, C. Davatzikos & G. Biros. *CLAIRE: A distributed-memory solver for constrained large deformation diffeomorphic image registration*. SIAM Journal on Scientific Computing 41(5):C548--C584, 2019 [[arxiv](https://arxiv.org/abs/1808.04487), [sisc](https://epubs.siam.org/doi/abs/10.1137/18M1207818)].
* M. Brunn, N. Himthani, G. Biros, M. Mehl & A. Mang. *CLAIRE: Constrained Large Deformation Diffeomorphic Image Registration on Parallel Computing Architectures*. Journal of Open Source Software, 6(61), 3038, 2021. [[joss](https://doi.org/10.21105/joss.03038)]

A BibTeX entry for LaTeX users are
```TeX
@article{Mang:2019a,
author = "A. Mang and A. Gholami and C. Davatzikos and G. Biros",
title = "{CLAIRE}: A distributed-memory solver for constrained large deformation diffeomorphic image registration",
journal = "SIAM Journal on Scientific Computing",
volume = "41",
number = "5",
pages = "C548--C584",
year = "2019"}
```

```TeX
@article{Brunn:2021a,
author = "M. Brunn and N. Himthani and G. Biros and M. Mehl and A. Mang",
title = "{CLAIRE}: Constrained large deformation diffeomorphic image registration on parallel computing architectures",
journal = "Journal of Open Source Software",
volume = "6",
number = "61",
pages = "3038",
year = "2021"}
```

If you want to refer to particular implementation or development aspects consider citing one of the works below.

## Parallel CPU Implementation
* A. Mang, A. Gholami, C. Davatzikos & G. Biros. *CLAIRE: A distributed-memory solver for constrained large deformation diffeomorphic image registration*. SIAM Journal on Scientific Computing 41(5):C548--C584, 2019 [[arxiv](https://arxiv.org/abs/1808.04487), [sisc](https://epubs.siam.org/doi/abs/10.1137/18M1207818)].
* A. Mang, A. Gholami & G. Biros. *Distributed-memory large-deformation diffeomorphic 3D image registration*. Proc ACM/IEEE Conference on SuperComputing, #72, 2016 [[arxiv](https://arxiv.org/abs/1608.03630), [ieee](http://dx.doi.org/10.1109/SC.2016.71)].


## Parallel GPU Implementation
* M. Brunn, N. Himthani, G. Biros, M. Mehl & A. Mang. *Multi-node multi-GPU diffeomorphic image registration for large-scale imaging problems*. Proc ACM/IEEE Conference on Supercomputing 2020. [[arxiv](https://arxiv.org/abs/2008.12820), [ieee](https://doi.ieeecomputersociety.org/10.1109/SC41405.2020.00042)].
* M. Brunn, N. Himthani, G. Biros, M. Mehl & A. Mang. *Fast GPU 3D diffeomorphic image registration*. Journal of Parallel and Distributed Computing, 149:149-162, 2021 [[arxiv](https://arxiv.org/abs/2004.08893), [jpdc](https://doi.org/10.1016/j.jpdc.2020.11.006)].



## Algorithmic Developments
* A. Mang & G. Biros. *An inexact Newton-Krylov algorithm for constrained diffeomorphic image registration*. SIAM Journal on Imaging Sciences, 8(2):1030--1069, 2015. [[arxiv](https://arxiv.org/abs/1408.6299v3), [siims](http://epubs.siam.org/doi/10.1137/140984002)].
* A. Mang & G. Biros. *Constrained H1 regularization schemes for diffeomorphic image registration*. SIAM Journal on Imaging Sciences, 9(3):1154--1194, 2016 [[arxiv](https://arxiv.org/abs/1503.00757), [siims](http://epubs.siam.org/doi/10.1137/15M1010919)].
* A. Mang & G. Biros. *A semi-Lagrangian two-level preconditioned Newton-Krylov solver for constrained diffeomorphic image registration*. SIAM Journal on Scientific Computing, 39(6):B1064--B1101, 2017. [[arxiv](https://arxiv.org/abs/1604.02153), [sisc](http://epubs.siam.org/doi/abs/10.1137/16M1070475)].
* M. Brunn, N. Himthani, G. Biros, M. Mehl & A. Mang. *CLAIRE: Constrained Large Deformation Diffeomorphic Image Registration on Parallel Computing Architectures*. Journal of Open Source Software, 6(61), 3038, 2021. [[joss](https://doi.org/10.21105/joss.03038)]
# CLAIRE: The Binaries

Go back to [README.md](../README.md).

## Content
* [Overview](#clairebins)
* [Get Help](#clairehelp)
* [Input Data](#clairedata)
* [Simple Examples: `claire`](#claireexmp)
	* [Synthetic Problem](#clairexmp1)
	* [Synthetic Problem (Parallel Execution)](#clairexmp2)
	* [Real Data (Parallel Execution)](#clairexmp3)
	* [Regularization Parameter Estimation](#clairexmp4)
	* [Parameter Continuation](#clairexmp5)
	* [Output Velocities](#clairexmp6)
* [Simple Examples: `clairetools`](#toolsxmp)
	* [Transporting Images](#toolsxmp1)
	* [Computing Jacobians](#toolsxmp2)
	* [Converting Data](#toolsxmp3)
* [Job Submission on Dedicated Systems](#hpc)
* [Testing and Benchmarks](#testing)


## Overview <a name="clairebins"></a>

CLAIRE is a software for 3D diffeomorphic image registration. CLAIRE has two main binaries: `claire` and `clairetools`.

  * `claire`: perform registrations
  * `clairetools`: post and pre-processing

We provide **several examples** for executing these binaries in the [doc/examples](https://github.com/andreasmang/claire/tree/gpu/examples) subfolder. We briefly explain these examples below.

In addition to that we have added two binaries named `test` and `benchmark` for developers to test the code. They are described in greater detail in the [Testing and Benchmarks](#testing) section below.

These binaries can be found in the `bin` folder after CLAIRE has been built successfully. To learn more about building claire take a look at our [installation guide](README-INSTALL.md) found in [doc/README-INSTALL.md](README-INSTALL.md).

If you are trying to execute CLAIRE on a dedicated HPC system (multi-GPU / multi-CPU envoriment) take a look at the examples provided in the [Job Submission on Dedicated Systems](#hpc) section.

## Get Help <a name="clairehelp"></a>

To learn more about the options available in `claire` and `clairetools` add the `-help` flag:

```bash
$BINDIR/claire -help
$BINDIR/clairetools -help
```

`$BINDIR` is the directory in which the binary is located.

## Input Data <a name="clairedata"></a>

CLAIRE is a software for 3D diffeomorphic image registration. Supported input formats for CLAIRE are images stored as '*.nii', '*.nii.gz' or '*.hdr' and '*.img(.gz)'. This is the default data format for inputs and outputs generated by CLAIRE.


## Simple Examples: `claire` <a name="claireexmp"></a>

### Example 01: Synthetic Problem <a name="clairexmp1"></a>

In [runclaire01.sh](https://github.com/andreasmang/claire/tree/gpu/doc/examples/runclaire01.sh) we execute CLAIRE for a synthetic test problem of size 32x32x32. We use default settings for our solver:
```bash
$BINDIR/claire -synthetic 0
```

`$BINDIR` is the directory in which the binary is located. The flag `-synthetic` allows one to select several smooth test problems. To change the problem size simply add the `-nx <n1xn2xn3>` flag, where `<n$i$>` represents the problem size in each spatial direction (i.e., `-nx 128x128x128` executes CLAIRE with a problem size of `128x128x128`.) We recommend executing CLAIRE in parallel for larger problem sizes (see [example 2](#clairexmp2))

### Example 02: Synthetic Problem (Parallel Execution) <a name="clairexmp2"></a>

In [runclaire02.sh](examples/runclaire02.sh) we execute CLAIRE for a synthetic test problem of size 128x128x128 in parallel. As an example, we use 20 MPI tasks. We use default settings for our solver:

```bash
mpirun -n 20 $BINDIR/claire -synthetic 0 -nx 128
```

The options used with `claire` are explained in [example 1](#clairexmp1). The key difference is the instruction `mpirun -np 20` infront of the executable. This instructs your compute node to use 20 MPI tasks. CLAIRE will determine the processor layout for you. We recommend executing CLAIRE in parallel. If you use `mpiexec` replace `mpirun -np 20` with `mpiexec -n 20`.


### Example 03: Real Data <a name="clairexmp3"></a>

In [runclaire03.sh](https://github.com/andreasmang/claire/tree/gpu/doc/examples/runclaire03.sh) we execute CLAIRE for real medical images (in NIfTI format) of size 128x150x128. We use 20 MPI tasks. The data can be found in the [docs/data](data) subdirectory. We use default settings for our solver:

```bash
$BINDIR/claire -mr $datdir/brain01.nii.gz \
               -mt $datdir/brain02.nii.gz
```

**<span style="color:red">Important:</span>** The images have to be **affinely preregistered** (i.e., have the same grid size).

`$datdir` points to the location where the data is located. The `-mr` flag identifies the image to be used as a **reference image** (alas, *fixed* or *target* image) and the `-mt` flag identifies the image to be used as a **template image** (i.e., the image to be registered; alas *floating* or *moving* image). The line break (backslash `\`) is only added for readability.


### Example 04: Regularization Parameter Estimation <a name="clairexmp4"></a>

In [runclaire04.sh](https://github.com/andreasmang/claire/tree/gpu/doc/examples/runclaire04.sh) we execute CLAIRE to automatically identify an adequate regularization parameter for a given set of images. We use default settings for our solver:

```bash
$BINDIR/claire -mr $datdir/brain01.nii.gz \
               -mt $datdir/brain02.nii.gz -train binary
```

Running `claire` on real image data is explained in [example 3](#clairexmp3). We use a method based on parameter continuation to identify an adequate regularization parameter. The search uses the determinant of the deformation gradient as a "metric". The user can define a lower bound for the Jacobian via the `-jbound <dbl>` option (the upper bound is `1/<dbl>`). To perform the search specify the `-train <type>` option. There are two strategies implemented: A simple reduction of the regularization parameter until the bound is hit (use `-train reduce`) and a binary search (use `-train binary`). The line break (backslash `\`) is only added for readability.

Notice that if you compile CLAIRE in single precision only H1-type regularization operators function properly (there are numerical accuracy issues for H2- and H3-type regularization operators in single precision). To use higher order regularization operators, CLAIRE needs to be compiled in double precision (slower).

### Example 05: Parameter Continuation <a name="clairexmp5"></a>

In [runclaire05.sh](https://github.com/andreasmang/claire/tree/gpu/doc/examples/runclaire05.sh) we show how to execute CLAIRE using a parameter continuation scheme with a target regularization parameter for the velocity. We use default settings for our solver:

```bash
$BINDIR/claire -mr $datdir/brain01.nii.gz \
               -mt $datdir/brain02.nii.gz \
               -betacont 7.75e-04
```

We have observed that a parameter continuation scheme speeds up the rate of convergence of our solver. We recommend using it in practical settings. We show how to estimate an adequate regularization parameter in  [example 4](#clairexmp4). The line breaks (backslashes `\`) are only added for readability.


### Example 06: Output Velocities <a name="clairexmp6"></a>

In [runclaire06.sh](https://github.com/andreasmang/claire/tree/gpu/doc/examples/runclaire06.sh) we show how to store the computed velocity field on file. We use default settings for our solver:

```bash
$BINDIR/claire -mr $datdir/brain01.nii.gz \
               -mt $datdir/brain02.nii.gz \
               -betacont 7.75e-04  -x ./ -velocity
```

This example is a direct extension of [example 4](#clairexmp4). The only difference is that we added an output. We need to provide an output folder. This is done with the `-x <folder>` option. We write the output to the current directory (`./`). **The outputs in CLAIRE will have default names**. If you prefer to store all files for multiple runs in a single folder, we recommend to use a prefix:
```bash
-x /my/output/folder/name/PREFIX_
```

The `-velocity` option tells CLAIRE to write out the velocity field. There are multiple other outputs available, most of which can be computed from the velocity field. This can be done using `clairetools`. To learn more about how to use `clairetools` continue reading. The line breaks (backslashes `\`) are only added for readability.


## Simple Examples: `clairetools` <a name="toolsxmp"></a>

### Example 01: Transporting Images <a name="toolsxmp1"></a>

In [runtools01.sh](https://github.com/andreasmang/claire/tree/gpu/doc/examples/runtools01.sh) we show how to transport an image (i.e., e.g., compute the deformed template image after a velocity has been computed using `claire`.)

```bash
$BINDIR/clairetools -v1 velocity-field-x1.nii.gz       \
                    -v2 velocity-field-x2.nii.gz       \
                    -v3 velocity-field-x3.nii.gz       \
                    -ifile $datdir/brain01.nii.gz      \
                    -xfile brain01-transported.nii.gz -deformimage
```

The input are the three components of the computed velocity (`-v$i$ velocity-field-x$i$.nii.gz `) and the image to be transported (`-ifile $datdir/brain01.nii.gz`; `$datdir` points to the folder the data is located in, i.e., [doc/data](https://github.com/andreasmang/claire/tree/gpu/doc/data)). The output is the transported brain image (`-xfile brain01-transported.nii.gz`). The user can add a path as prefix if desired. The command to tell `clairetools` that we are interested in solving the forward problem (i.e., transporting/deforming an image) is `-deformimage`. The line breaks (backslashes `\`) are only added for readability.


### Example 02: Computing Jacobians <a name="toolsxmp2"></a>

In [runtools02.sh](https://github.com/andreasmang/claire/tree/gpu/doc/examples/runtools02.sh) we show how to compute the determinant of the deformation gradient (alas Jacobian) from a velocity field that has been computed using `claire`.

```bash
$BINDIR/clairetools -v1 velocity-field-x1.nii.gz       \
                    -v2 velocity-field-x2.nii.gz       \
                    -v3 velocity-field-x3.nii.gz       \
                    -x ./ -detdefgrad
```


### Example 03: Converting Data <a name="toolsxmp3"></a>

In [runtools03.sh](https://github.com/andreasmang/claire/tree/gpu/doc/examples/runtools03.sh) we show how to convert data from `*.nii.gz` to `*.nc`.

```bash
$BINDIR/clairetools -ifile $datdir/brain01.nii.gz      \
                    -convert 2nc
```

The output file will be located in the same folder as the input file, with the same filename (for example, `$datdir/brain01.nii.gz` becomes `$datdir/brain01.nc`). Notice that you need to compile CLAIRE with nifti and pnetcdf enabled.


## Job Submission on Dedicated Systems <a name="hpc"></a>


### TACC's Longhorn System

More information about TACC's Longhorn system can be found at [https://www.tacc.utexas.edu/systems/longhorn](https://www.tacc.utexas.edu/systems/longhorn).

To execute CLAIRE on a single GPU node, interactively, do the following:
```bash
idev -N1
ibrun $BINDIR/claire -synthetic 0 -nx 128
```

An example for a script for executing CLAIRE on Longhorn can be found here:
* multi-GPU (2 GPUs one node): [doc/examples/longhorn_mgpu.slurm](examples/longhorn_mgpu.slurm)


## Testing and Benchmarks <a name="testing"></a>

We have implemented numerical tests for the main computational kernels available in CLAIRE to study the performance and accuracy of the mathematical operators that appear in the optimality system. For reproducability, we also posted the NIREP data at [https://github.com/andreasmang/nirep](https://github.com/andreasmang/nirep). We have used this data extensively in our [prior work](README-REFERENCES.md).

To build **binaries for testing** set `BUILD_TEST=yes` in the `makefile` to compile CLAIRE (see [makefile](../makefile)). This will build two binaries: `benchmark` and `test`. To inspect available options use the `-help` flag
```bash
$BINDIR/test -help
$BINDIR/benchmarks -help
```

`$BINDIR` is the directory in which the binary is located.


The `test` application allows users to check the main computational kernels:
* the interpolation kernels (to check the interpolation accuracy)
```bash
$BINDIR/test [other args] -interp
```
* the regularization operators (e.g., biharmonic or laplacian regularization operators)
```bash
$BINDIR/test [other args] -reg
```
* numerical differentiation (to check the accuracy of the differentiation operators)
```bash
$BINDIR/test [other args] -diff
```
* correctness of gradient and hessian (see below for a more detailed description; can also be used within `claire`) 
```bash
$BINDIR/test [other args] -gradient
$BINDIR/test [other args] -hessian
```
* accuracy of trajectory computation (semi-lagrangian time integration scheme) 
```bash
$BINDIR/test [other args] -trajectory
```
These tests are implemented in the `*.cpp` files available in the [UnitTests](https://github.com/andreasmang/claire/tree/gpu/src/UnitTests) subfolder.

We have also implemented several **high-level numerical checks** to assess the performance of our methodology and ensure that the mathematical operators are correct.

The `benchmark` binary allows users to (i) check the accuracy of the forward operator and (ii) report the runtime for evaluating several key mathematical operators. Use the `-help` flag to see all options. The main tests are the following:
* check error of forward operator (solve the forward problem for $v$ and $-v$ and check error with respect to initial condition for $t=0$)
```bash
$BINDIR/benchmark [other args] -terror
```
* report runtimes for evaluating the forward operator, the gradient operator and the Hessian matvec.

```bash
$BINDIR/benchmark [other args] -forward
$BINDIR/benchmark [other args] -gradient
$BINDIR/benchmark [other args] -hessmatvec
```

The **tests/debug options** directly available within the `claire` binary are the following (Use the `-help` flag to see all options. Notice that the help provides a `other parameters/debugging` section that describes the options mentioned below):
* The default test in CLAIRE is to consider synthetic test problems. The user can select between several test problems of varying complexity by setting the flag `-synthetic i`, where `i` selects the particular test case (valid values for `i` are `0`, `1`, ..., `5`).
```bash
$BINDIR/claire [other args] -synthetic 0 
```
* The user can control the verbosity level of `claire` by using the `-verbose 2` flag (debug mode verbosity). This will, e.g., enable command window outputs such as the residual in each iteration of the Krylov subspace method used to compute the search direction (and much more).
```bash
$BINDIR/claire [other args] -verbose 2
```
* The Newton--Krylov solver monitors several critical values during the course of the iterations. The user can see outputs such as
	* number of (Gauss--)Newton iterations.
	* trend of the objective value (has to decrease monotonically)
	* trend of the gradient norm (should decrease but not necessarily monotonically)
	* number of line search steps (should be 1 for a Newton method subject to accuracy requirements)
	* and much more...
* The accuracy of the symmetry of the discretized Hessian operator can be monitored by enabling the `-checksymmetry` flag in `claire`: Notice that we consider an optimize-then-discretize approach and numerical schemes that do not preserve symmetry; consequently, in our current implementation, the Hessian is only symmetric up to the discretization error of the adjoint operators (~1e-2).
```bash
$BINDIR/claire [other args] -checksymmetry
```
* The approximation accuracy of the gradient and Hessian can be monitored by enabling the `-derivativecheck` flag in `claire`. We report the assymptotic behavior of the Taylor expansion. The approximation error should decrease with decreasing perburbation (i.e., the error should converge). However, we do, in general not expect to observe quadratic or cubic convergence (as we would if we considered a discretize-then-optimize approach).
```bash
./bin/claire [other args] -derivativecheck
```
# CLAIRE: Contributor Covenant Code of Conduct

Go back to [README.md](../README.md).

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, sex characteristics, gender identity and expression,
level of experience, education, socio-economic status, nationality, personal
appearance, race, religion, or sexual identity and orientation.

## Our Standards

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

## Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

## Scope

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting [Andreas Mang](mailto:andreas@math.uh.edu). All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see
https://www.contributor-covenant.org/faq
# CLAIRE: Contributing Guidlines

Go back to [README.md](../README.md).

## Content
* [Overview](#overview)
* [Reporting Bugs and Issues](#bugs)
* [Feature Requests](#features)
* [Contributing to CLAIRE](#contribute)
* [Testing and Benchmarks](#testing)
* [Coding Conventions](#conventions)

## Overview <a name="overview"></a>

We welcome all contributions to CLAIRE. Contributions can be in the form of new code functions, improvements to the documentation, or by pointing out a bug or potential improvement. For general questions or requests, you can contact [Andreas Mang](http://math.uh.edu/~andreas) by email: [andreas [at] math [dot] uh.edu](mailto:andreas@math.uh.edu).

This project and everyone participating in it is governed by the code of conduct found in [doc/CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md). By participating, you are expected to uphold this code. Please report unacceptable behavior to [Andreas Mang](mailto:andreas@math.uh.edu).



## Reporting Bugs and Issues <a name="bugs"></a>

Bugs and issues can be reported at [https://github.com/andreasmang/claire/issues](https://github.com/andreasmang/claire/issues). Before submitting a ticket, make sure the bug has not already been reported by searching existing issues on GitHub under [issues](https://github.com/andreasmang/claire/issues). We also list some known issues in [doc/README-INSTALL.md](README-INSTALL.md).

If you are reporting a bug or an issue, please include detailed information to help maintainers reproduce the problem.


## Feature Requests <a name="features"></a>

Additional features can be requested at [https://github.com/andreasmang/claire/issues](https://github.com/andreasmang/claire/issues). Before submitting a feature request, make sure the feature has not already been requested by searching existing issues on GitHub under [issues](https://github.com/andreasmang/claire/issues).


## Contributing to CLAIRE <a name="contribute"></a>

To contribute to the software:

1. [Fork](https://docs.github.com/en/free-pro-team@latest/github/getting-started-with-github/fork-a-repo) the repository.
2. Clone the forked repository, add your contributions and push the changes to your fork.
3. Create a [pull request](https://github.com/andreasmang/claire/pulls).


## Testing and Benchmarks <a name="testing"></a>

We have implemented several tests to check the accuracy of our numerical implementation. These are described in more detail in [doc/README-RUNME.md](https://github.com/andreasmang/claire/blob/gpu/doc/README-RUNME.md#testing-and-benchmarks-).


## Coding Conventions <a name="conventions"></a>

Our source code follows the [Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html). Please adhere to these coding conventions if you would like to contribute to CLAIRE. Notice that all our routines use the error mechanism implemented in PETSc. We strongly encourage all contributers to include this mechanism to make debugging easier. We have added `tags` to identify the major releases of our software (related mostly to publications for reproducability).
# CLAIRE: News Feed

Go back to [README.md](../README.md).
* 05/2021 Our JOSS contribution has been accepted for publication. The reference can be found [here](doc/README-REFERENCES.md).
* 03/2021 We have added contributing guidlines: [doc/CONTRIBUTING.md](CONTRIBUTING.md)
* 03/2021 The `cpu` branch has become the default branch of the repository.
* 03/2021 We have renamed the `master` branch to `cpu` (the legacy CPU implementation of CLAIRE).
* 28/2020 We have released a GPU version of CLAIRE. If you are interested in using our new (multi-node multi-)GPU version, switch to the **GPU branch**. If you are interested in learning more about the GPU version of CLAIRE, check out our [publications](README-REFERENCES.md).
* 11/2018 Support for Normalized-Cross-Correlation (NCC) has been added to CLAIRE.
# CLAIRE: Installation and Requirements

Go back to [README.md](../README.md).

## Content

* [Installation Overview](#installation)
	* [One Shot](#oneshot)
	* [Step by Step](#stepbystep)
* [Detailed Installation Guide](#verboseinstall)
	* [Requirements](#requirements)
	* [Dependencies](#dependencies)
		* Required Dependencies and Compatibility
		* Step 1: Downloading and Installing Dependencies
		* Step 2: Setting Environment Variables
	* [Building CLAIRE](#buildclaire)
	* [Executing CLAIRE](#execclaire)
* [Additional Info for Dependencies](#depsinf)
* [CLAIRE on Specific Systems](#clairesys)
* [Troubleshooting and Known Issues](#faq)


## Installation Overview (Quick Guide)<a name="installation"></a>

In this section, we provide a minimal installation guide. We provide a make environment to download and install the dependencies using generic settings that have worked on most of our systems. If this brief installation guide does not work for you, please consult the [detailed installation guide](#verboseinstall) below.

### One Shot <a name="oneshot"></a>

To use the default settings to build dependencies and CLAIRE itself do the following:

```bash
cd deps
make
source env_source.sh
cd ..
make -j
./bin/claire -synthetic 0
```

The enviroment variables need to be sourced every time you log out of your computer or start a new bash (`source env_source.sh`). As an alternative, you can add the content of `env_source.sh` (for example) to your `.bashrc` or `.bash_profile`.

### Step by Step  <a name="stepbystep"></a>

Next, we go over the steps outlined above step by step. Again, more details are provided [below](#verboseinstall).

#### Step 1) Installing Dependencies

To install the dependencies (the PETSc and NIFTI libraries) go to the top level directory of CLAIRE in your command window and execute the following commands within your command window:

```bash
cd deps
make
```

This makefile downloads and compiles the dependencies for CLAIRE. To add these dependencies to your environment type the following into your command line and press return:

```bash
source env_source.sh
```

Notice that the enviroment variables need to be sourced every time you log out of your computer or start a new bash. As an alternative, you can add the content of `env_source.sh` to your `.bashrc` or `bash_profile`.


#### Step 2) Compiling CLAIRE

Assuming that you are in the top level directory of CLAIRE, all you need to do is to type

```bash
make -j
```

#### Step 3) Executing CLAIRE

If you would like to verify if CLAIRE has been installed correctly run the following command in your command window:

```bash
./bin/claire -synthetic 0
```

Additional examples for executing CLAIRE are described in [doc/README-RUNME.md](README-RUNME.md).


## Detailed Installation Guide <a name="verboseinstall"></a>

In this section we provide a more detailed description of the installation process to help users with troubleshooting. The following table lists system configurations on which we have successfully installed CLAIRE.

|Test   | Compiler  | MPI            | CUDA | PETSc  | CPU    | GPU   | System       |
|---    |---------- |-----           |------|------- |---     |---    |---           |
|b5213fa| GCC 9.3   | OpenMPI 4.0.3  | 11.0 | 3.14.2 | x86_64 | GA102 | Ubuntu 20.04 |
|6f40316| GCC 9.3   | OpenMPI 4.0.3  | 11.1 | 3.14.2 | x86_64 | GK110 | Ubuntu 20.04 |
|4967052| GCC 8.4   | OpenMPI 1.10.2 | 10.1 | 3.12.4 | x86_64 | GK110 | Ubuntu 16.04 |
|4967052| GCC 5.4.0 | OpenMPI 1.10.2 | 10.0 | 3.12.4 | x86_64 | GM200 | Ubuntu 16.04 |
|4967052| GCC 7.4   | OpenMPI 4.0.1  | 10.1 | 3.12.4 | x86_64 | GP100 | Ubuntu 16.04 |
|4967052| GCC 4.8.5 | OpenMPI 3.1.6  | 10.2 | 3.12.4 | Power9 | GV100 | CentOS 7.8   |
|4967052| XLC 16.1  | Spectrum 10.3  | 10.2 | 3.12.4 | Power9 | GV100 | RHEL 7.8     |


### Requirements <a name="requirements"></a>

The minimal requirements for compiling CLAIRE on your system are:
* MPI (Open MPI; MVAPICH; Intel MPI; ...; required by [PETSc](https://www.mcs.anl.gov/petsc), and CLAIRE)
* cmake (see [https://cmake.org](https://cmake.org); required by niftilib)
* python (see [https://www.python.org](https://www.python.org); required by [PETSc](https://www.mcs.anl.gov/petsc) and the optional pyclaire bindings)
* zlib (see [https://www.zlib.net](https://www.zlib.net); required by niftilib)
* CUDA-API

Make sure that the standard *wrappers* for `mpicc`, `mpicxx`, and `nvcc` are available on your system (either by loading the appropriate modules and/or by setting up the appropriate `PATH` and `LD_LIBRARY_PATH` definitions). The compilation has been tested with Open MPI, MVAPICH, and Intel MPI.



### Dependencies <a name="dependencies"></a>


#### Required Dependencies and Compatibility

The compiler needs `C++11` support. The GPU version of CLAIRE requires the following libraries to be installed on your system:

* MPI (with GPU support (CUDA-aware MPI) for multi-GPU multi-node)
* PETSc with CUDA support (see [https://www.mcs.anl.gov/petsc](https://www.mcs.anl.gov/petsc))
* niftilib (see [https://sourceforge.net/projects/niftilib/files/nifticlib](https://sourceforge.net/projects/niftilib/files/nifticlib))
* zlib (see [http://zlib.net](http://zlib.net))

We provide functionality to build PETSc, niftilib, and zlip on your system (see next section).


#### Step 1: Downloading and Installing Dependencies

To download and compile the libraries we provide a `makefile` (see [deps/makefile](../deps/makefile)). Simply run `make` with this script in your command window to download *tarball* files of the libraries identified above.

```bash
cd deps
make
```

The *compressed* tarball files (i.e, `LIBRARY-NAME.tar.gz`) should remain located in or be added to the [deps](../deps) folder. Make sure that all libraries are downloaded (the progress bar of `wget` should be full). To view the urls for the libraries you can take a look at the [deps/makefile](../deps/makefile). We provide additional information about these libraries [below](#depsinf). This also includes links to versions for these libraries that we have used to compile the GPU version of CLAIRE before.

The [makefile](../deps/makefile) has some optional parameters to configure the build. The parameters can be set by `make PARAMETER=value`. Multiple parameter-value pairs can be passed to the make command. The [makefile](../deps/makefile) to compile the dependencies has the following parameters.

| PARAMETER       | Description                                           | Default | Valid Values  |
| --------------- | ----------------------------------------------------- | ------- | ------        |
| BUILD_PETSC     | PETSc version to download and compile; empty for none | 3.12.4  | PETSc Version |
| BUILD_NIFTI     | Download and build `niftilib`                         | yes     | yes, no       |
| WITH_BATCH      | Option to build petsc on a batch system, e.g. slurm   | no      | yes, no       |
| WITH_CUDA_MPI   | MPI is CUDA-aware                                     | yes     | yes, no       |
| CC              | Path to C compiler                                    | mpicc   | file path     |
| CXX             | Path to CXX compiler                                  | mpicxx  | file path     |
| NVCC            | Path to CUDA compiler                                 | nvcc    | file path     |
| WITH_PETSC_OPTS | additional PETSC compile options                      |         |               |

The libraries will be extracted and build in the `deps/lib` subfolder.


#### Step 2: Setting Environment Variables

Before you are able to compile and run CLAIRE you need to add *environment variables* to your system. When building the libraries a file called `env_source.sh` is created. This file should be located in the [debs](../deps) subfolder. To add the environment variables temporarily (for the current session) to your system, do

```bash
source env_source.sh
```

To add them permanently, copy the content of `env_source.sh` to your `~/.bashrc`. Notice that `env_source.sh` defines *absolute paths*.


## Building CLAIRE <a name="buildclaire"></a>

Before you can build CLAIRE you need to

* Make sure that you have installed all *dependencies* (see prior sections).
* Make sure all paths and compilers needed in the `makefile` are available on your system, i.e. `mpicxx`, `nvcc`, and the dependencies.

To inspect all options used in the `makefile` for CLAIRE (see [makefile](../makefile)) do (in the top level directory):

```bash
make VERBOSE=1 VVERBOSE=1 config
```

To build the code using the `make` system do (in the top level directory):

```bash
make -j
```

If you build in parallel using `make -j`, on certain systems to many threads will be used. This will result in compilation errors. To fix this, run `make -j 12` instead (for quick access, you may want to define an alias in your `~/.bashrc`).

The [makefile](../makefile) also contains optional parameters to configure the build. The parameters can be set by `make PARAMETER=value`. Multiple parameter-value pairs can be passed to the `make` command. The `makefile` to compile the dependencies has following parameters.

| PARAMETER      | Description                                           | Default | Valid Values  |
| -------------- | ----------------------------------------------------- | ------- | ------        |
| BUILD_TEST     | build the unit test applications                      | no      | yes; no       |
| BUILD_PYTHON   | build `pyclaire` python bindings                      | no      | yes; no       |
| BUILD_SHARED   | build CLAIRE as shared library                        | no      | yes; no       |
| WITH_NIFTI     | build with `niftilib`                                 | yes     | yes; no       |
| WITH_DEBUG     | build with additional debug informations              | no      | yes; no       |
| WITH_DEVELOP   | build CLAIRE additional development informations      | no      | yes; no       |
| WITH_CUDA_MPI  | MPI is CUDA-aware                                     | yes     | yes, no       |
| BUILD_TARGET   | target CPU architecture                               | X86     | X86; POWER9   |
| GPU_VERSION    | GPU CUDA version to compile, e.g. 35, 60, 70, 75      |         | Compute Capability |
| CPP_VERSION    | C++ Standard to use                                   | c++11   | c++11; c++14  |
| LD_FLAGS       | additional flags for the linker                       |         |               |
| CXX_FLAGS      | additional flags for the C++ compiler                 |         |               |
| NVCC_FLAGS     | additional flags for the CUDA compiler                |         |               |
| MPI_DIR        | main path to the MPI include and lib directory        |         |               |
| CUDA_DIR       | main path to the CUDA include and lib directory       |         |               |
| PETSC_DIR      | main path to the PETSc include and lib directory      |         |               |
| NIFTI_DIR      | main path to the libnifti include and lib directory   |         |               |
| ZLIB_DIR       | main path to the zlib include and lib directory       |         |               |
| PYTHON_DIR     | main path to the Python3 include and lib directory    |         |               |
| VERBOSE        | if set to any value the make command is verbose       |         |               |
| VVERBOSE       | if set to any value the make command is very verbose  |         |               |

Not that the `makefile` generates a cache (`make.cache`) to detect if a complete rebuild of CLAIRE is needed. If this file is removed or does not exsit the next build will first do `make clean` automatically.


## Executing CLAIRE <a name="execclaire"></a>

If you would like to verify if CLAIRE has been installed correctly, run the following command in your command window:

```bash
./bin/claire -synthetic 0
```

Additional examples for executing CLAIRE are described in [doc/README-RUNME.md](README-RUNME.md)


## Additional Info for Dependencies <a name="depsinf"></a>

### PETSc
* PETSc webpage: [https://www.mcs.anl.gov/petsc](https://www.mcs.anl.gov/petsc)
* description: library for numerics, linear algebra, and optimization
* source code also available on bitbucket: [https://bitbucket.org/petsc/petsc](https://bitbucket.org/petsc/petsc)
* versions that have been succesfully used by our group:
	* [petsc-lite-3.11.4.tar.gz](http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.11.4.tar.gz)

### nifticlib
* NIFTICLIB webpage: [http://niftilib.sourceforge.net](http://niftilib.sourceforge.net)
* description: library to read and write NIFTI images
* see [https://sourceforge.net/projects/niftilib/files/nifticlib](https://sourceforge.net/projects/niftilib/files/nifticlib)
* versions that have been succesfully used by our group:
	* [nifticlib-2.0.0.tar.gz](https://sourceforge.net/projects/niftilib/files/nifticlib/nifticlib_2_0_0)



## CLAIRE on Specific Systems <a name="clairesys"></a>

### TACC's Longhorn System (03/17/21)

More information about TACC's Longhorn system can be found at [https://www.tacc.utexas.edu/systems/longhorn](https://www.tacc.utexas.edu/systems/longhorn).

Modules loaded:
```bash
1) xl/16.1.1             4) autotools/1.2   7) TACC
2) spectrum_mpi/10.3.0   5) cmake/3.16.1    8) cuda/10.2 (g)
3) git/2.24.1            6) xalt/2.10.2
```

Compilation of CLAIRE and its dependencies:

```bash
cd deps
make WITH_BATCH=yes
source deps/env_source.sh
make BUILD_TARGET=POWER9 GPU_VERSION=70
```

To test if the compilation worked check if the binaries are available in the `bin` folder. To execute CLAIRE using an interactive job do

```bash
cd bin
idev -N1 # launch an interactive session with one node
ibrun ./claire -help
ibrun ./claire -synthetic 0 -nx 128
```

A job submission file for TACC's Longhorn system (for multi-GPU exection) can be found in [doc/examples/longhorn_mgpu.slurm](examples/longhorn_mgpu.slurm).



## Troubleshooting / Known Issues <a name="faq"></a>

* if MPI is not compiled with CUDA-aware options, add the file `.petscrc` to the working directory and add the option `-use_gpu_aware_mpi 0`
* CUDA >= 11.0 is only supported with PETSc >= 3.14.
* Kepler GPUs work with PETSc 3.12.4  (others not tested)
* Compiling PETSc with CUDA support on cluster login nodes without GPUs might fail
* PNETCDF is currently not tested for GPUs
* The GPU version of CLAIRE can currently only be compiled in single precision. This limits the selection of regularization operators to H1-type regularization only. There are issues with the numerical accuracy of H2- and H3-type regularization operators for single precision. Applying these operators requires a compilation in double precision (available on the GPU branch) 
---
title: 'CLAIRE: Constrained Large Deformation Diffeomorphic Image Registration on Parallel Computing Architectures'
tags:
  - C++
  - GPUs
  - parallel computing
  - high performance computing
  - large deformation diffeomorphic registration
  - optimal control
  - medical imaging
authors:
  - name: Malte Brunn
    affiliation: "1"
  - name: Naveen Himthani
    affiliation: "2"
  - name: George Biros
    affiliation: "2"
  - name: Miriam Mehl
    affiliation: "1"
  - name: Andreas Mang^[Corresponding Author]
    orcid: 0000-0003-4605-3290
    affiliation: "3"
affiliations:
 - name: Institute for Parallel and Distributed Systems, University Stuttgart
   index: 1
 - name: Oden Institute for Computational Engineering and Sciences, The University of Texas at Austin
   index: 2
 - name: Department of Mathematics, University of Houston
   index: 3
date: 28 May 2021
bibliography: paper.bib
---

# Summary
[`CLAIRE`](https://andreasmang.github.io/claire) [@claire-web] is a computational framework for **C**onstrained **LA**rge deformation diffeomorphic **I**mage **RE**gistration [@Mang:2019a]. It supports highly-optimized, parallel computational kernels for (multi-node) CPU [@Mang:2016a; @Gholami:2017a; @Mang:2019a] and (multi-node multi-)GPU architectures [@Brunn:2020a; @Brunn:2021a]. `CLAIRE` uses MPI for distributed-memory parallelism and can be scaled up to thousands of cores [@Mang:2019a; @Mang:2016a] and GPU devices [@Brunn:2020a]. The multi-GPU implementation uses device direct communication. The computational kernels are interpolation for semi-Lagrangian time integration, and a mixture of high-order finite difference operators and Fast-Fourier-Transforms (FFTs) for differentiation. `CLAIRE` uses a Newton--Krylov solver for numerical optimization [@Mang:2015a; @Mang:2017a]. It features various schemes for regularization of the control problem [@Mang:2016a] and different similarity measures. `CLAIRE` implements different preconditioners for the reduced space Hessian [@Brunn:2020a; @Mang:2019a] to optimize computational throughput and enable fast convergence. It uses `PETSc` [@petsc-web] for scalable and efficient linear algebra operations and solvers and `TAO` [@petsc-web; @Munson:2015a] for numerical optimization. `CLAIRE` can be downloaded at <https://github.com/andreasmang/claire>.

# Statement of Need
Image registration is required whenever images are taken at different points in time, from different viewpoints, and/or using different imaging modalities and these images need to be compared, combined, or integrated [@Fischer:2008a; @Modersitzki:2004a; @Modersitzki:2009a; @Sotiras:2013a]. Image registration is an inverse problem. The inputs to this inverse problem are two (or more) images $m_0(x)$ (the template image) and $m_1(x)$ (the reference image) of the same object. The task of image registration is to find a plausible map $y(x)$ that establishes spatial correspondences between the reference and template image, i.e., $m_0(x) \approx m_1(y(x))$. In `CLAIRE` the set of admissible spatial transformations $y$ is limited to diffeomorphisms, i.e., maps $y$ that are continuous, one-to-one, and have a smooth inverse. `CLAIRE` is related to a prominent class of formulations for these types of problems referred to as <em>large-deformation diffeomorphic metric mapping</em> [@Beg:2005a; @Trouve:1998a; @Younes:2010a].

Diffeomorphic image registration is an indispensable tool in medical image analysis [@Sotiras:2013a]. Computing diffeomorphisms that map one image to another is expensive. Deformable image registration is an infinite-dimensional problem that upon discretization leads to nonlinear optimality systems with millions or even billions of unknowns. For example, registering two typical medical imaging datasets of size $256^3$ necessitates solving for about 50 million unknowns (in our formulation). Additional complications are the ill-posedness and non-linearty of this inverse problem [@Fischer:2008a]. Consequently, image registration can take several minutes on multi-core high-end CPUs. Many of the available methods reduce the number of unknowns by using coarser resolutions either through parameterization or by solving the problem on coarser grids; they use simplified algorithms and deliver subpar registration quality. In the age of big data, clinical population studies that require thousands of registrations are incresingly common, and execution times of individual registrations become more critical. We provide technology that allows solving registration problems for clinical datasets in seconds. In addition, we have made available to the public a software that works on multi-node, multi-GPU architectures [@Brunn:2020a; @Brunn:2021a] that allows the registration of large-scale microscopic imaging data such as CLARITY imaging [@Tomer:2014a; @Kutten:2017a].

# Highlights
`CLAIRE` can be used to register images of $2048^3$ (25 B unknowns) on 64 nodes with 256 GPUs on TACC’s Longhorn system [@Brunn:2020a]. `CLAIRE` has been used for the registration of high resolution CLARITY imaging data [@Brunn:2020a]. The GPU version of `CLAIRE` can solve clinically relevant problems (50 M unknowns) in approximately 5 seconds on a single NVIDIA Tesla V100 [@Brunn:2020a]. `CLAIRE` has also been applied to hundreds of images in brain tumor imaging studies [@Bakas:2018a; @Mang:2017c; @Scheufele:2021a], and has been integrated with models for biophysics inversion [@Mang:2018a; @Mang:2020a; @Scheufele:2020a; @Scheufele:2019a; @Scheufele:2021a; @Subramanian:2020b] and Alzheimer's disease progression [@Scheufele:2020c]. `CLAIRE` uses highly optimized computational kernels and effective, state-of-the-art algorithms for time integration and numerical optimization. Our most recent version of `CLAIRE` features a Python interface to assist users in their applications.

We provide a detailed documentation on how to execute, compile, and install `CLAIRE` on various systems at our deployment page <https://andreasmang.github.io/claire>.

# Mathematics
`CLAIRE` uses an optimal control formulation. The diffeomorphism $y(x)$ is parameterized using a smooth, stationary velocity field $v(x)$. Given the template image $m_0(x)$ and the reference image $m_1(x)$, this velocity is found by solving the partial-differential equation constrained optimization problem of the form
$$
\operatorname{minimize}_{v,m} \operatorname{dist}(m(x,t=1),m_1) + \alpha\operatorname{reg}(v)
$$

subject to
$$
\begin{aligned}
\partial_t  m(x,t) + v(x) \cdot \nabla m(x,t) &= 0 \\
m(x,t=0) & = m_0(x)
\end{aligned}
$$

The first term in the objective functional measures the proximity of the deformed template image $m(x,t=1)$ and the reference image $m_1(x)$. The default option availble in `CLAIRE` is an $L^2$-distance. The second term controls the regularity of $v$. `CLAIRE` features different Sobolev norms. The default option is an $H^1$-seminorm. The constraint models the deformation the template image (i.e., the transport of the intensities of $m_0(x)$). `CLAIRE` also features additional hard constraints for controlling the divergence of $v(x)$ [@Mang:2016a]. For optimization, we use the method of Lagrange multipliers and solve the associated Karush--Kuhn--Tucker optimality system using a Newton--Krylov reduced space method [@Mang:2015a; @Mang:2015a].


# Acknowledgements
This work was partly supported by the National Science Foundation (DMS-1854853, DMS-2009923, DMS-2012825, CCF-1817048, CCF-1725743), the NVIDIA Corporation (NVIDIA GPU Grant Program), the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) under Germany’s Excellence Strategy-EXC 2075-390740016, by the U.S. Department of Energy, Office of Science, Office of Advanced Scientific Computing Research, Applied Mathematics program under Award Number DE-SC0019393; by the U.S. Air Force Office of Scientific Research award FA9550-17-1-0190; by the Portugal Foundation for Science and Technology and the UT Austin-Portugal program, and by NIH award 5R01NS042645-11A1. Any opinions, findings, and conclusions or recommendations expressed herein are those of the authors and do not necessarily reflect the views of the DFG, AFOSR, DOE, NIH, and NSF. Computing time on the Texas Advanced Computing Centers’ (TACC) systems was provided by an allocation from TACC and the NSF. This work was completed in part with resources provided by the Research Computing Data Core at the University of Houston.

# References
