# Change Log
All notable changes to this project will be documented in this file.

* This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
* The format of this log is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [Unreleased]

[Documentation](https://cr-sparse.readthedocs.io/en/latest/)

## [0.2.2] - 2021-12-02

[Documentation](https://cr-sparse.readthedocs.io/en/v0.2.2/)

### Improved

- Documentation
  - Introduction page revised
  - API docs improved
  - README revised
  - Algorithm page added
  - Quick start page revised

### Added

- JOSS Paper added and revised based on feedback from reviewers
- Linear Operators
  - input_shape, output_shape attributes introduced
  - fft, dot, norm estimate, reshape, scalar_mult, total variation, 


## [0.2.1] - 2021-11-01

[Documentation](https://cr-sparse.readthedocs.io/en/v0.2.1/)

### Added

- Compressive Sensing
  - 1 bit compressive sensing process
  - BIHT (Binary Iterative Hard Thresholding) algorithm for signal reconstruction from 1 bit measurements

## [0.2.0] - 2021-10-30

[Documentation](https://cr-sparse.readthedocs.io/en/v0.2.0/)

### Added

- Linear Operators
  - Convolution 1D, 2D, ND
  - Gram and Frame operators for a given linear operator
  - DWT 1D operator
  - DWT 2D operator
  - Block diagonal operator (by combining one or more operators)
- Sparse Linear Systems
  - Power iterations for computing the largest eigen value of a symmetric linear operator
  - LSQR solver for least squares problems with support for N-D data
  - ISTA: Iterative Shrinkage and Thresholding Algorithm
  - FISTA: Fast Iterative Shrinkage and Thresholding Algorithm
  - lanbpro, simple lansvd
- Geophysics
  - Ricker wavelet
  - Hard, soft and half thresholding operators for ND arrays (both absolute and percentile thresholds)
- Image Processing
  - Gaussian kernels
- Examples
  - Deconvolution
  - Image Deblurring
- Data generation
  - Random subspaces, uniform points on subspaces
  - two_subspaces_at_angle, three_subspaces_at_angle
  - multiple index_sets
  - sparse signals with bi-uniform non-zero values 
- Utilities
  - More functions for ND-arrays
  - Off diagonal elements in a matrix, min,max, mean
  - set_diagonal, abs_max_idx_cw, abs_max_idx_rw
- Linear Algebra
  - orth, row_space, null_space, left_null_space, effective_rank
  - subspaces: principal angles, is_in_subspace, project_to_subspace
  - mult_with_submatrix, solve_on_submatrix
  - lanbpro, simple lansvd
- Clustering
  - K-means clustering
  - Spectral clustering
  - Clustering error metrics
- Subspace clustering
  - OMP for sparse subspace clustering
  - Subspace preservation ratio metrics 

A paper is being prepared for JOSS.

### Improved

- Linear Operators
  - Ability to apply a 1D linear operator along a specific axis of input data
  - axis parameter added to various compressive sensing operators
- Code coverage
  - It is back to 90+% in the unit tests



## [0.1.6] - 2021-08-29

[Documentation](https://cr-sparse.readthedocs.io/en/v0.1.6/)

### Added

Wavelets
- CWT implementation based on PyWavelets: CMOR and MEXH
- integrate_wavelet, central_frequency, scale2frequency

Examples
- CoSaMP step by step
- Chirp CWT with Mexican Hat Wavelet
- Frequency Change Detection using DWT
- Cameraman Wavelet Decomposition


### Changed

Wavelets
- CWT API has been revised a bit.

### Updated

Examples
- Sparse recovery via ADMM

Signal Processing
- frequency_spectrum, power_spectrum

## [0.1.5] - 2021-08-22

[Documentation](https://cr-sparse.readthedocs.io/en/v0.1.5/)

### Added

Linear Operators
- Orthogonal basis operators: Cosine, Walsh Hadamard
- General Operators: FIR, Circulant, First Derivative, Second Derivative, Running average
- Operators: Partial Op
- DOT TEST for linear operators added.

Convex optimization algorithms
- Sparsifying basis support in yall1
- TNIPM (Truncated Newton Interior Points Method) implemented.

Wavelets
- Forward DWT
- Inverse DWT
- Padding modes: symmetric, reflect, constant, zero, periodic, periodization
- Wavelet families: HAAR, DB, MEYER, SYMMETRIC, COIFLET
- Full DWT/IDWT for periodization mode
- Filtering with Upsampling/Downsampling
- Quadrature Mirror Filters
- Forward and inverse DWT along a specific axis
- 2D Forward and inverse DWT for images
- Print wavelet info
- wavelist
- families
- build_wavelet
- wavefun
- CWT for Morlet and Ricker wavelets

Benchmarking
- Introduced airspeed velocity based benchmarks

General stuff
- Examples gallery introduced
- Unit test coverage is now back to 90%
- Documentation has been setup at ReadTheDocs.org also https://cr-sparse.readthedocs.io/en/latest/


## [0.1.4] - 2021-07-12

[Documentation](https://cr-sparse.readthedocs.io/en/v0.1.4/)

### Added

- A framework for linear operators
- ADMM based algorithms for l1 minimization
- Several greedy algorithms updated to support linear operators as well as plain matrices
- Hard thresholding pursuit added

## [0.1.3] - 2021-06-06
### Added
- Subspace Pursuit
- Iterative Hard Thresholding

## [0.1.0] - 2021-06-05

Initial release

[Unreleased]: https://github.com/carnotresearch/cr-sparse/compare/v0.2.2...HEAD
[0.2.2]: https://github.com/carnotresearch/cr-sparse/compare/v0.2.1...v0.2.2
[0.2.1]: https://github.com/carnotresearch/cr-sparse/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/carnotresearch/cr-sparse/compare/v0.1.6...v0.2.0
[0.1.6]: https://github.com/carnotresearch/cr-sparse/compare/v0.1.5...v0.1.6
[0.1.5]: https://github.com/carnotresearch/cr-sparse/compare/v0.1.4...v0.1.5
[0.1.4]: https://github.com/carnotresearch/cr-sparse/compare/0.1.3...v0.1.4
[0.1.3]: https://github.com/carnotresearch/cr-sparse/compare/v0.1...0.1.3
[0.1.0]: https://github.com/carnotresearch/cr-sparse/releases/tag/v0.1# Contributor Covenant Code of Conduct

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
reported by contacting the project team at contact@carnotresearch.com. All
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
# Guidelines for contributing code

## What do I need to know to help?

If you are looking to help to with a code contribution our project is primarily written in Python 
and uses Numpy, JAX and related technologies. 
If you don't feel ready to make a code contribution yet, no problem! 
You can also check out the documentation issues 
or the design issues that we have.

If you are interested in making a code contribution 
and would like to learn more about the technologies that we use, 
check out the list below.

* [JAX reference documentation](https://jax.readthedocs.io/en/latest/)
* [SciPy Lectures](https://scipy-lectures.org/)
* [CR-Sparse reference documentation](https://carnotresearch.github.io/cr-sparse/)

## How do I make a contribution?

Never made an open source contribution before? Wondering how contributions work in the in our project? Here's a quick rundown!

* Find an issue that you are interested in addressing or a feature that you would like to add.

* Fork the repository associated with the issue to your local GitHub organization. 
  This means that you will have a copy of the repository under your-GitHub-username/cr-sparse.

* Clone the repository to your local machine using 
  `git clone https://github.com/your-user-name/cr-sparse.git`.

* Create a new branch for your fix using `git checkout -b branch-name-here`.

* Make the appropriate changes for the issue you are trying to address or the feature that you want to add.

* Use `git add insert-paths-of-changed-files-here` to add the file contents of the changed files to the 
  "snapshot" git uses to manage the state of the project, also known as the index.

* Use `git commit -m "Insert a short message of the changes made here"` to store the contents of the index with a descriptive message.

* Push the changes to the remote repository using `git push origin branch-name-here`.

* Submit a pull request to the upstream repository.

* Title the pull request with a short description of the changes made and the issue or bug number associated with your change. 
  For example, you can title an issue like so "Added more log outputting to resolve #4352".

* In the description of the pull request, explain the changes that you made, any issues you think exist 
  with the pull request you made, and any questions you have for the maintainer. 
  It's OK if your pull request is not perfect (no pull request is), the reviewer will be able to help you 
  fix any problems and improve it!

* Wait for the pull request to be reviewed by a maintainer.

* Make changes to the pull request if the reviewing maintainer recommends them.

* Celebrate your success after your pull request is merged!

## Where can I go for help?

If you need help, you can ask questions on our [Discussions](https://github.com/carnotresearch/cr-sparse/discussions) forum for this project.

## What does the Code of Conduct mean for me?

Our Code of Conduct means that you are responsible for treating everyone 
on the project with respect and courtesy regardless of their identity. 
If you are the victim of any inappropriate behavior or 
comments as described in our Code of Conduct, 
we are here for you and will do the best to ensure that 
the abuser is reprimanded appropriately, per our code.
---
title: 'CR-Sparse: Hardware accelerated functional algorithms for sparse signal processing in Python using JAX'
tags:
    - Python
    - sparse and redundant representations
    - compressive sensing
    - wavelets
    - linear operators
    - sparse subspace clustering
    - functional programming
authors:
    - name: Shailesh Kumar
      affiliation: 1
      orcid: 0000-0003-2217-4768

affiliations:
    - name: Indian Institute of Technology, Delhi
      index: 1


date: 29 November 2021
bibliography: paper.bib
---

# Summary

We introduce [`CR-Sparse`](https://github.com/carnotresearch/cr-sparse), 
a Python library that enables to efficiently solve a wide variety of sparse representation based signal processing problems.
It is a cohesive collection of sub-libraries working together. Individual
sub-libraries provide functionalities for:
wavelets, linear operators, greedy and convex optimization 
based sparse recovery algorithms, subspace clustering, 
standard signal processing transforms,
and linear algebra subroutines for solving sparse linear systems. 
It has been built using Google JAX [@jax2018github], which enables the same high level
Python code to get efficiently compiled on CPU, GPU and TPU architectures
using XLA [@abadi2017computational]. 

![Sparse signal representations and compressive sensing](./srr_cs.png)

Traditional signal processing exploits the underlying structure in signals
by representing them using Fourier or wavelet orthonormal bases. 
In these representations,
most of the signal energy is concentrated in few coefficients allowing greater
flexibility in analysis and processing of signals. More flexibility can be
achieved by using overcomplete dictionaries [@mallat2008wavelet]
(e.g. unions of orthonormal bases). However, the construction of
sparse representations of signals in these overcomplete dictionaries 
is no longer straightforward and requires use of specialized sparse
coding algorithms like orthogonal matching pursuit [@pati1993orthogonal]
or basis pursuit [@chen2001atomic]. The key idea behind these algorithms 
is the fact that under-determined systems $A x = b$ can be solved efficiently
to provide sparse solutions $x$ if the matrix $A$ satisfies specific conditions
on its properties like coherence. Compressive sensing takes the same 
idea in the other direction and contends that signals having sparse representations
in suitable bases can be acquired by very few data-independent 
random measurements $y = \Phi x$ if the sensing or measurement system $\Phi$
satisfies certain conditions like restricted isometry property [@candes2008restricted].
The same sparse coding algorithms can be tailored for sparse signal recovery
from compressed measurements. 

A short mathematical introduction to compressive sensing and sparse representation problems 
is provided in [docs](https://cr-sparse.readthedocs.io/en/latest/intro.html).
For comprehensive introduction to sparse
representations and compressive sensing,
 please refer to excellent books [@mallat2008wavelet;@elad2010sparse;@foucart2013mathintro],
papers [@donoho2006compressed;@qaisar2013compressive;@marques2018review],
[Rice Compressive Sensing Resources](https://dsp.rice.edu/cs/) and references therein.

# Package Overview

The `cr.sparse.pursuit` package includes greedy and thresholding
based solvers for sparse recovery. It includes: 
`OMP`, `CoSaMP`, `HTP`, `IHT`, `SP` algorithms.
(provided in `cr.sparse.lop` package).
The `cr.sparse.cvx` package includes efficient solvers
for l1-minimization problems using convex optimization methods.
The `cr.sparse.sls` package provides JAX versions of
`LSQR`, `ISTA`, `FISTA`  algorithms for solving sparse linear 
systems.
These algorithms can work with unstructured random and dense sensing matrices
as well as structured sensing matrices represented as linear operators
The `cr.sparse.lop` package includes a collection of linear operators
influenced by `PyLops` [@ravasi2019pylops]. 
`cr.sparse.wt` package includes a JAX version of major functionality
from `PyWavelets` [@lee2019pywavelets] making it a first major pure 
Python wavelets implementation which can work across CPUs, GPUs and TPUs.


# Statement of need

Currently, there is no single Package which provides a 
comprehensive set of tools for solving sparse recovery problems
in one place. Individual researchers provide their codes
along with their research paper only for the algorithms they have
developed. Most of this work is available in the form of MATLAB [@MATLAB:2018]
libraries. E.g.: [`YALL1`](http://yall1.blogs.rice.edu) is the original MATLAB implementation of the ADMM based sparse recovery algorithms. 
[`L1-LS`](https://web.stanford.edu/~boyd/l1_ls/) 
is the original MATLAB implementation of the
Truncated Newton Interior Points Method for solving the l1-minimization problem.
[`Sparsify`](https://www.southampton.ac.uk/engineering/about/staff/tb1m08.page#software) 
provides the MATLAB implementations of IHT, NIHT, AIHT algorithms.
[`aaren/wavelets`](https://github.com/aaren/wavelets) 
is a CWT implementation following
[@torrence1998practical]. 
[`HTP`](https://github.com/foucart/HTP) provides implementation of Hard Thresholding
Pursuit in MATLAB.
[`WaveLab`](https://github.com/gregfreeman/wavelab850) is the 
reference open source wavelet implementation in MATLAB.
However, its API has largely been superseded by later libraries.
[`Sparse and Redundant Representations book code`](https://elad.cs.technion.ac.il/wp-content/uploads/2018/02/Matlab-Package-Book-1.zip) [@elad2010sparse]
provides basic implementations of a number of sparse recovery and related
algorithms.
Several of these libraries contain key 
performance critical sub-routines
in the form of C/C++ extensions making portability to GPUs harder. 

There are some Python libraries which focus on specific
areas however they are generally CPU based.
E.g., [`pyCSalgos`](https://github.com/nikcleju/pyCSalgos) is a
Python implementation of various Compressed Sensing algorithms.
[`spgl1`](https://github.com/drrelyea/spgl1) is a `NumPy` based
implementation of spectral projected gradient for L1 minimization.
`c-lasso` [@simpson2021classo] is a Python package for constrained sparse regression 
and classification. This is also CPU only. 
[`PyWavelets`](https://github.com/PyWavelets/pywt) is an excellent 
CPU only wavelets implementation in Python closely following the API
of Wavelet toolbox in MATLAB. The performance critical parts have been
written entirely in C. There are several attempts to port it on GPU
using `PyTorch` ([PyTorch-Wavelet-Toolbox](https://github.com/v0lta/PyTorch-Wavelet-Toolbox)) 
or `Tensorflow` ([tf-wavelets](https://github.com/UiO-CS/tf-wavelets)) backends.
[`PyLops`](https://github.com/PyLops/pylops) includes GPU support. 
They have built a [`backend.py`](https://github.com/PyLops/pylops/blob/master/pylops/utils/backend.py) 
layer to switch explicitly between
`NumPy` and [`CuPy`](https://cupy.dev/) for GPU support. 
In contrast, 
our use of JAX enables us to perform jit compilation with 
abstracted out end-to-end XLA optimization to multiple backend.

The algorithms in this package have a wide variety of applications. We list
a few: image denoising, deblurring, compression, inpainting, impulse noise removal,
super-resolution,
subspace clustering, dictionary learning, 
compressive imaging, medical imaging, compressive radar,
wireless sensor networks, astrophysical signals, cognitive radio,
sparse channel estimation, analog to information conversion, 
speech recognition, seismology, direction of arrival.


# Sparse signal processing problems and available solvers

We provide JAX based implementations for the following algorithms:

* `cr.sparse.pursuit.omp`: Orthogonal Matching Pursuit (OMP) [@pati1993orthogonal;@tropp2004greed;@davenport2010analysis] 
* `cr.sparse.pursuit.cosamp`: Compressive Sampling Matching Pursuit (CoSaMP) [@needell2009cosamp]
* `cr.sparse.pursuit.sp`: Subspace Pursuit (SP) [@dai2009subspace]
* `cr.sparse.pursuit.iht`: Iterative Hard Thresholding and its normalized version (IHT, NIHT) [@blumensath2009iterative;@blumensath2010normalized]
* `cr.sparse.pursuit.htp`: Hard Thresholding Pursuit and its normalized version (HTP, NHTP) [@foucart2011recovering]
* `cr.sparse.cvx.l1ls`: *Truncated Newton Interior Points Method* for solving the l1-minimization problem [@kim2007interior]
* `cr.sparse.cvx.admm`: Solvers for basis pursuit (BP), basis pursuit denoising (BPDN), basis pursuit with inequality constraints (BPIC),
  and their nonnegative variants based on ADMM [@yang2011alternating;@zhang2009user]
* `cr.sparse.sls.lsqr`: LSQR algorithm for sparse linear equations [@paige1982lsqr]
* `cr.sparse.sls.ista`: Iterative Shrinkage and Thresholding Algorithm (ISTA) [@daubechies2004iterative]
* `cr.sparse.sls.fista`: Fast Iterative Shrinkage Thresholding Algorithm (FISTA) [@beck2009fast]

The dictionaries and sensing matrices can be efficiently implemented 
using a pair of functions for the forward $A x$ and adjoint $A^H x$ operations.
`cr.sparse.lop` provides a collection of linear operators (similar to `PyLops` [@ravasi2019pylops])
which provide the forward and adjoint operation functions. 
These operators can be JIT compiled and used efficiently with the algorithms above.
Our 2D and ND operators accept 2D/ND arrays as input and return 2D/ND arrays as output.
The operators `+`, `-`, `@`, `**` etc. are overridden to provide operator calculus, 
i.e. ways to combine operators to generate new operators.

As an application area, the library includes an implementation of
sparse subspace clustering (SSC) by orthogonal matching pursuit
[@you2016scalable] in the `cr.sparse.cluster.ssc` package.
The `cr.sparse.cluster.spectral` package provides a custom implementation of spectral clustering
step of SSC.

# Experimental Results

We conducted a number of experiments to benchmark the runtime of 
`CR-Sparse` implementations viz. existing reference software
in Python or MATLAB. 
Jupyter notebooks to reproduce these micro-benchmarks
are available on the 
[`cr-sparse-companion`](https://github.com/carnotresearch/cr-sparse-companion) [@shailesh2021companion]
repository.

All Python based benchmarks have been run
on the machine configuration: 
Intel(R) Xeon(R) Gold 6130 CPU @ 2.10GHz, 16 Cores, 64 GB RAM, 
NVIDIA GeForce GTX 1060 6GB GPU, 
Ubuntu 18.04 64-Bit, Python 3.8.8, 
NVidia driver version 495.29.05,
CUDA version 11.5.

MATLAB based benchmarks were run on the machine configuration:
Intel(R) Core(TM) i7-10510U CPU @ 1.80GHz   2.30 GHz,
32 GB RAM, Windows 10 Pro, MATLAB R2020b.

The following table provides comparison of `CR-Sparse` against reference 
implementations on a set of representative problems:

\footnotesize

| Problem | Size | Ref tool | Ref time | Our time | Gain |
|:----:|:---:|:--:|:-:|:-:|:-:|
| Hard Thresholding Pursuit | M=2560, N=10240, K=200 | HTP (MATLAB) | 3.5687 s | 160 ms | 22x |  
| Orthogonal Matching Pursuit | M=2000, N=10000, K=100 | sckit-learn | 379 ms | 120 ms | 3.15x |  
| ADMM, BP | M=2000, N=20000, K=200 | YALL1 (MATLAB) | 1.542 sec | 445 ms | 3.46x |  
| ADMM, BPDN | M=2000, N=20000, K=200 | YALL1 (MATLAB) | 1.572.81 sec | 273 ms | 5.75x |  
| Image blurring | Image: 500x480, Kernel: 15x25 | Pylops | 6.63 ms | 1.64 ms | 4x |  
| Image deblurring using LSQR | Image: 500x480, Kernel: 15x25 | Pylops | 237 ms | 39.3 ms | 6x |  
| Image DWT2 | Image: 512x512 | PyWavelets | 4.48 ms | 656 µs | 6.83x |  
| Image IDWT2 | Image: 512x512 | PyWavelets | 3.4 ms | 614 µs | 5.54x |  
| OMP for SSC | 5 subspaces 50K points | SSCOMP_Code (MATLAB) | 52.5 s | 10.2 s | 4.6x |

\normalsize


We see significant though variable gains achieved by `CR-Sparse` on GPU. 
We have observed that gain tends to increase for larger problem sizes. 
GPUs tend to perform better when problem size increases as the matrix/vector 
products become bigger.
`vmap` and `pmap` tools provided by JAX can be used to easily 
parallelize the `CR-Sparse`  algorithms over multiple data and processors.

Following table compares the runtime of linear operators in `CR-Sparse` on GPU vs 
`PyLops` on CPU for large size problems. 
Timings are measured for both forward and adjoint operati`ons. 

\footnotesize

| Operator | Size | Fwd ref | Fwd our | Gain | Adj ref | Adj our | Gain |
|:--------:|:----:|:-:|:-:|:-:|:-:|:-:|:-:|
| Diagonal matrix mult| n=1M | 966 µs | 95.7 µs | 10x | 992 µs | 96.3 µs | 10x | 
| Matrix mult | (m,n)=(10K,10K) | 11 ms | 2.51 ms | 4.37x | 11.6 ms | 2.51 ms | 4.63x |
| First derivative | n=1M | 2.15 ms | 71.1 µs | 30.2x | 2.97 ms | 186 µs | 15.97x |
| HAAR DWT2, level=8 | in=(4K,4K) | 981 ms | 34.4 ms | 28.5x | 713 ms | 60.8 ms | 11.7x | 

\normalsize

# Limitations

Some of the limitations in the library come from the underlying 
JAX library. 
JAX is relatively new and still hasn't 
reached `1.0` level maturity. 
The programming model chosen by JAX places
several restrictions on expressing the program logic. For example,
JAX does not have support for dynamic or data dependent shapes
in their [JIT compiler](https://jax.readthedocs.io/en/latest/notebooks/thinking_in_jax.html#to-jit-or-not-to-jit).
Thus, any algorithm parameter which determines the size/shape
of individual arrays in an algorithm must be statically provided.
E.g. for the greedy algorithms like OMP, 
the sparsity level $K$ must be known in advance
and provided as a static parameter to the API as the size of
output array depends on $K$. 

The control flow primitives like `lax.while_loop`, `lax.fori_loop`
etc. in JAX require that the algorithm state flowing between iterations
must not change shape and size. This makes coding of algorithms
like OMP or SVT (singular value thresholding) very difficult.
An incremental QR or Cholesky decomposition based implementation of OMP requires
growing algorithm state. We ended up using a standard Python `for` loop
for now but the JIT compiler simply unrolls it and doesn't allow for tolerance
based early termination in them. 

1D convolutions are slow in JAX on CPU 
[#7961](https://github.com/google/jax/discussions/7961). 
This affects the performance of DWT/IDWT in `cr.sparse.dwt`. 
We are working on exploring ways of making it more efficient 
while keeping the API intact.

These restrictions necessitate good amount of creativity and a very
disciplined coding style so that efficient JIT friendly 
solvers can be developed. 

# Future Work

Currently, work is underway to provide a JAX based
implementation of [`TFOCS`](http://cvxr.com/tfocs/) [@becker2011templates]
in the dev branch.
This will help us increase the coverage to a wider set of
problems (like total variation minimization, Dantzig selector,
l1-analysis, nuclear norm minimization, etc.). As part of this
effort, we are expanding our collection of linear operators and
building a set of indicator and projector functions on to
convex sets and proximal operators [@parikh2014proximal].
This will enable us to cover other applications such as
SSC-L1 [@pourkamali2020efficient]. 
In future, we intend to increase the coverage in following areas:
More recovery algorithms (OLS, Split Bergmann, SPGL1, etc.) 
and specialized cases (partial known support, );
Bayesian Compressive Sensing;
Dictionary learning (K-SVD, MOD, etc.);
Subspace clustering;
Image denoising, compression, etc. problems using sparse representation principles;
Matrix completion problems;
Matrix factorization problems;
Model based / Structured compressive sensing problems;
Joint recovery problems from multiple measurement vectors.

# Acknowledgements

Shailesh would like to thank his Ph.D. supervisors Prof. Surendra 
Prasad and Prof. Brejesh Lall to inculcate his interest in this area
and support him over the years in his exploration. He would also
like to thank his employers, Interra Systems Inc., for allowing him
to pursue his research interests along with his day job.


# References

---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior:
1. Go to '...'
2. Click on '....'
3. Scroll down to '....'
4. See error

**Expected behavior**
A clear and concise description of what you expected to happen.


**Additional context**
Add any other context about the problem here.
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: ''
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
```{include} ../CHANGELOG.md
```