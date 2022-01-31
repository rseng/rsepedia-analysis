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
```Functional Models and Algorithms for Sparse Signal Processing   
==================================================================


|pypi| |license| |zenodo| |docs| |unit_tests| |coverage| |joss|


Introduction
-------------------


CR-Sparse is a Python library that enables efficiently solving
a wide variety of sparse representation based signal processing problems.
It is a cohesive collection of sub-libraries working together. Individual
sub-libraries provide functionalities for:
wavelets, linear operators, greedy and convex optimization 
based sparse recovery algorithms, subspace clustering, 
standard signal processing transforms,
and linear algebra subroutines for solving sparse linear systems. 
It has been built using `Google JAX <https://jax.readthedocs.io/en/latest/>`_, 
which enables the same high level
Python code to get efficiently compiled on CPU, GPU and TPU architectures
using `XLA <https://www.tensorflow.org/xla>`_. 

.. image:: docs/images/srr_cs.png

For detailed documentation and usage, please visit `online docs <https://cr-sparse.readthedocs.io/en/latest>`_.

Supported Platforms
----------------------

``CR-Sparse`` can run on any platform supported by ``JAX``. 
We have tested ``CR-Sparse`` on Mac and Linux platforms and Google Colaboratory.

``JAX`` is not officially supported on Windows platforms at the moment. 
Although, it is possible to build it from source using Windows Subsystems for Linux.

Installation
-------------------------------

Installation from PyPI:

.. code:: shell

    python -m pip install cr-sparse

Directly from our GITHUB repository:

.. code:: shell

    python -m pip install git+https://github.com/carnotresearch/cr-sparse.git



Examples/Usage
----------------

See the `examples gallery <https://cr-sparse.readthedocs.io/en/latest/gallery/index.html>`_ in the documentation.
Here is a small selection of examples:

* `Sparse recovery using Truncated Newton Interior Points Method <https://cr-sparse.readthedocs.io/en/latest/gallery/rec_l1/spikes_l1ls.html>`_ 
* `Sparse recovery with ADMM <https://cr-sparse.readthedocs.io/en/latest/gallery/rec_l1/partial_wh_sensor_cosine_basis.html>`_ 
* `Compressive sensing operators <https://cr-sparse.readthedocs.io/en/latest/gallery/lop/cs_operators.html>`_ 
* `Image deblurring with LSQR and FISTA algorithms <https://cr-sparse.readthedocs.io/en/latest/gallery/lop/deblurring.html>`_ 
* `Deconvolution of the effects of a Ricker wavelet <https://cr-sparse.readthedocs.io/en/latest/gallery/lop/deconvolution.html>`_ 
* `Wavelet transform operators <https://cr-sparse.readthedocs.io/en/latest/gallery/lop/wt_op.html>`_ 
* `CoSaMP step by step <https://cr-sparse.readthedocs.io/en/latest/gallery/pursuit/cosamp_step_by_step.html>`_ 


A more extensive collection of example notebooks is available in the `companion repository <https://github.com/carnotresearch/cr-sparse-companion>`_.
Some micro-benchmarks are reported `here <https://github.com/carnotresearch/cr-sparse/blob/master/paper/paper.md#runtime-comparisons>`_.


Contribution Guidelines/Code of Conduct
----------------------------------------

* `Contribution Guidelines <CONTRIBUTING.md>`_
* `Code of Conduct <CODE_OF_CONDUCT.md>`_

Citing CR-Sparse
------------------------


To cite this library:

.. code:: tex

    @article{Kumar2021,
      doi = {10.21105/joss.03917},
      url = {https://doi.org/10.21105/joss.03917},
      year = {2021},
      publisher = {The Open Journal},
      volume = {6},
      number = {68},
      pages = {3917},
      author = {Shailesh Kumar},
      title = {CR-Sparse: Hardware accelerated functional algorithms for sparse signal processing in Python using JAX},
      journal = {Journal of Open Source Software}
    }




`Documentation <https://carnotresearch.github.io/cr-sparse>`_ | 
`Code <https://github.com/carnotresearch/cr-sparse>`_ | 
`Issues <https://github.com/carnotresearch/cr-sparse/issues>`_ | 
`Discussions <https://github.com/carnotresearch/cr-sparse/discussions>`_ |


.. |docs| image:: https://readthedocs.org/projects/cr-sparse/badge/?version=latest
    :target: https://cr-sparse.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status
    :scale: 100%

.. |unit_tests| image:: https://github.com/carnotresearch/cr-sparse/actions/workflows/ci.yml/badge.svg
    :alt: Unit Tests
    :scale: 100%
    :target: https://github.com/carnotresearch/cr-sparse/actions/workflows/ci.yml


.. |pypi| image:: https://badge.fury.io/py/cr-sparse.svg
    :alt: PyPI cr-sparse
    :scale: 100%
    :target: https://badge.fury.io/py/cr-sparse

.. |coverage| image:: https://codecov.io/gh/carnotresearch/cr-sparse/branch/master/graph/badge.svg?token=JZQW6QU3S4
    :alt: Coverage
    :scale: 100%
    :target: https://codecov.io/gh/carnotresearch/cr-sparse


.. |license| image:: https://img.shields.io/badge/License-Apache%202.0-blue.svg
    :alt: License
    :scale: 100%
    :target: https://opensource.org/licenses/Apache-2.0

.. |codacy| image:: https://app.codacy.com/project/badge/Grade/36905009377e4a968124dabb6cd24aae
    :alt: Codacy Badge
    :scale: 100%
    :target: https://www.codacy.com/gh/carnotresearch/cr-sparse/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=carnotresearch/cr-sparse&amp;utm_campaign=Badge_Grade

.. |zenodo| image:: https://zenodo.org/badge/323566858.svg
    :alt: DOI
    :scale: 100%
    :target: https://zenodo.org/badge/latestdoi/323566858

.. |joss| image:: https://joss.theoj.org/papers/ebd4e5ca27a5db705b1dc382b64e0bed/status.svg
    :alt: JOSS
    :scale: 100%
    :target: https://joss.theoj.org/papers/ebd4e5ca27a5db705b1dc382b64e0bed
Organization of the test suite
================================
Example Notebooks
===========================


* `Browse all notebooks <https://nbviewer.jupyter.org/github/carnotresearch/cr-sparse/tree/master/notebooks/>`_


Sparse Signal Models
-------------------------

* `Solutions of underdetermined systems <https://nbviewer.jupyter.org/github/carnotresearch/cr-sparse/blob/master/notebooks/ssm/nb_underdetermined_systems.ipynb>`_


Greedy Algorithms
-----------------------

* `CoSaMP Step by Step <https://nbviewer.jupyter.org/github/carnotresearch/cr-sparse/blob/master/notebooks/pursuit/cosamp/cosamp_step_by_step.ipynb>`_



Experiments
===========================


* `Browse all experiment notebooks <https://nbviewer.jupyter.org/github/carnotresearch/cr-sparse/tree/master/notebooks/experiments/>`_


Sparse Recovery
-------------------------

* `Recovery Performance Comparison of Various Sparse Recovery Algorithms <https://nbviewer.jupyter.org/github/carnotresearch/cr-sparse/blob/master/notebooks/experiments/pursuit/comparison/nb_recovery_comparison_multiple_methods.ipynb>`_


Airspeed Velocity Benchmarks
====================================


This directory contains benchmarks written in the structure supported by airspeed velocity.

We assume that you have Airspeed Velocity installed. If not, try::

    pip install asv


To run the benchmarks (against the master), run::

    asv run

If you are running it for the first time on a machine, it will ask you a set of questions to capture the details of the machine.
We assume that that you have ``conda`` available. A ``conda`` environment is constructed under ``.asv/env`` directory to run the benchmarks.
The benchmarking results wlil be stored inside ``.asv/results`` directory. 

To see the results for the master, run::

    asv show master
.. _gallery:

Examples Gallery
===========================

.. _wavelets_gallery:

Wavelets
------------------------
.. _recovery_l1_gallery:

Sparse Recovery via L1 minimization
-----------------------------------------------

.. _cs_gallery:

Compressive Sensing
-----------------------------------
.. _cluster_gallery:

Data Clustering
------------------------
.. _lop_gallery:

Linear Operators
-----------------------------------------------
Development
=================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   dev/sourcecode
   dev/limitations
   benchmarks/index
   changelog.md
Acronyms
======================



.. list-table::
    :widths: 20 80
    :header-rows: 1

    * - Acronym
      - Full form
    * - ADC
      - Analog to Digital Conversion
    * - AIC
      - Analog to Information Conversion
    * - AIHT
      - Accelerated Iterative Hard Thresholding
    * - AMP
      - Approximate Message Passing
    * - AST
      - Affine Scaling Transformation
    * - BAOMP
      - Backtracking based Adaptive Orthogonal Matching Pursuit
    * - BCQP
      - Bound Constrained Quadratic Program
    * - BCS
      - Bayesian Compressive Sensing
    * - BOCP
      - Block Orthogonal Matching Pursuit
    * - BP
      - Basis Pursuit
    * - BPDN
      - Basis Pursuit DeNoising
    * - BPIC
      - Basis Pursuit with Inequality Constraints
    * - CG
      - Conjugate Gradient
    * - CoSaMP
      - Compressed Sampling Matching Pursuit
    * - CP
      - Chaining Pursuit
    * - CR
      - Cognitive Radio
    * - CS
      - Compressive Sensing, Compressed Sensing, Compressive Sampling
    * - DCT
      - Discrete Cosine Transform
    * - DOA
      - Direction Of Arrival
    * - D-OMP
      - Differential Orthogonal Matching Pursuit
    * - DS
      - Dantzig Selector
    * - EM
      - Expectation Maximization
    * - FBP
      - Forward Backward Pursuit
    * - FFT
      - Fast Fourier Transform
    * - FISTA
      - Fast Iterative Shrinkage Thresholding Algorithm
    * - FOCUSS
      - Focal Underdetermined System Solution
    * - GBP
      - Greedy Basis Pursuit
    * - GOAMP
      - Generalized Orthogonal Adaptive Matching Pursuit
    * - GP
      - Gradient Projection
    * - GPSR
      - Gradient Projection for Sparse Reconstruction
    * - HDTV
      - High Definition Television
    * - HHS
      - Heavy Hitters on Steroids
    * - HTP
      - Hard Thresholding Pursuit
    * - IHT
      - Iterative Hard Thresholding
    * - IoT
      - Internet of Things
    * - IP
      - Integer Programming, Interior Point methods
    * - IPM
      - Interior Points Method
    * - IRLS
      - Iterative Reweighted Least Squares
    * - IST
      - Iterative Shrinkage/Thresholding algorithms
    * - LARS
      - Least Angle Regression
    * - LASSO
      - Least Absolute Shrinkage and Selection Operator
    * - LP
      - Linear Programming
    * - LS
      - Least Squares
    * - LSP
      - Least Squares Program
    * - LSQR
      - Least Squares
    * - MM
      - Majorization Minimization
    * - MP
      - Matching Pursuit
    * - MRI
      - Magnetic Resonance Imaging
    * - NHTP
      - Normalized Hard Thresholding Pursuit
    * - NIHT
      - Normalized Iterative Hard Thresholding
    * - OLS
      - Orthogonal Least Squares
    * - OMP
      - Orthogonal Matching Pursuit
    * - PCG
      - Preconditioned Conjugate Gradient
    * - QCLP
      - Quadratically Constrained Linear Program
    * - QP
      - Quadratic Programming, Quadratic Program
    * - SOCP
      - Second Order Cone Program
    * - SP
      - Subspace Pursuit
    * - SSC
      - Sparse Subspace Clustering 
    * - TNIPM
      - Truncated Newton Interior Point Method
    * - WHT
      - Walsh Hadamard Transform
Quick Start
===================

|pypi| |license| |zenodo| |docs| |unit_tests| |coverage| |joss|


Platform Support
----------------------

``cr-sparse`` can run on any platform supported by ``JAX``. 
We have tested ``cr-sparse`` on Mac and Linux platforms and Google Colaboratory.

``JAX`` is not officially supported on Windows platforms at the moment. 
Although, it is possible to build it from source using Windows Subsystems for Linux.


Installation
-------------------------------

Installation from PyPI:

.. code:: shell

    python -m pip install cr-sparse



Directly from our GITHUB repository:

.. code:: shell

    python -m pip install git+https://github.com/carnotresearch/cr-sparse.git



Examples
----------------

* See the :ref:`examples gallery <gallery>`.
* A more extensive collection of example notebooks is available in the `companion repository <https://github.com/carnotresearch/cr-sparse-companion>`_.
* Some micro-benchmarks are reported `here <https://github.com/carnotresearch/cr-sparse/blob/master/paper/paper.md#runtime-comparisons>`_.


.. note::

    ``cr-sparse`` depends on its sister library `cr-nimble <https://github.com/carnotresearch/cr-nimble>`_.
    Normally, it would be installed automatically as a dependency. 
    You may want to install it directly from GITHUB if you need access to the latest code.

    .. code:: shell

        python -m pip install git+https://github.com/carnotresearch/cr-nimble.git


.. |docs| image:: https://readthedocs.org/projects/cr-sparse/badge/?version=latest
    :target: https://cr-sparse.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status
    :scale: 100%

.. |unit_tests| image:: https://github.com/carnotresearch/cr-sparse/actions/workflows/ci.yml/badge.svg
    :alt: Unit Tests
    :scale: 100%
    :target: https://github.com/carnotresearch/cr-sparse/actions/workflows/ci.yml


.. |pypi| image:: https://badge.fury.io/py/cr-sparse.svg
    :alt: PyPI cr-sparse
    :scale: 100%
    :target: https://badge.fury.io/py/cr-sparse

.. |coverage| image:: https://codecov.io/gh/carnotresearch/cr-sparse/branch/master/graph/badge.svg?token=JZQW6QU3S4
    :alt: Coverage
    :scale: 100%
    :target: https://codecov.io/gh/carnotresearch/cr-sparse


.. |license| image:: https://img.shields.io/badge/License-Apache%202.0-blue.svg
    :alt: License
    :scale: 100%
    :target: https://opensource.org/licenses/Apache-2.0

.. |codacy| image:: https://app.codacy.com/project/badge/Grade/36905009377e4a968124dabb6cd24aae
    :alt: Codacy Badge
    :scale: 100%
    :target: https://www.codacy.com/gh/carnotresearch/cr-sparse/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=carnotresearch/cr-sparse&amp;utm_campaign=Badge_Grade

.. |zenodo| image:: https://zenodo.org/badge/323566858.svg
    :alt: DOI
    :scale: 100%
    :target: https://zenodo.org/badge/latestdoi/323566858

.. |joss| image:: https://joss.theoj.org/papers/ebd4e5ca27a5db705b1dc382b64e0bed/status.svg
    :alt: JOSS
    :scale: 100%
    :target: https://joss.theoj.org/papers/ebd4e5ca27a5db705b1dc382b64e0bed
Theory
=================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   la/index
   sls/index
   ssc/index
   fwsp/index
   acronyms
CR-Sparse
=====================================

A JAX/XLA based library of accelerated models and algorithms for inverse problems in 
sparse representation and compressive sensing. 
`GITHUB <https://github.com/carnotresearch/cr-sparse>`_.


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   start
   intro
   tutorials/index
   source/index
   algorithms
   theory
   gallery/index
   zzzreference
   development



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
References
===================

.. bibliography::
    :cited:
Algorithms
=======================================

This section lists and organizes available and planned algorithms in the `CR-Sparse` package.


Sparse recovery algorithms
--------------------------------------------

See :cite:`marques2018review` for a review of sparse recovery algorithms.

.. rubric:: Convex relaxation algorithms

.. list-table::
    :widths: 70 10 10 10
    :header-rows: 1

    * - Algorithm
      - Acronym
      - Status
      - Docs
    * - Truncated Newton Interior Points Method
      - L1LS
      - done
      - :ref:`... <api:l1min:tnipm>`
    * - Basis Pursuit using ADMM
      - BP
      - done
      - :ref:`... <api:l1min:admmm>`
    * - Basis Pursuit Denoising using ADMM
      - BPDN
      - done
      - :ref:`... <api:l1min:admmm>`
    * - Basis Pursuit with Inequality Constraints using ADMM
      - BPIC
      - done
      - :ref:`... <api:l1min:admmm>`

.. rubric:: Greedy pursuit algorithms

.. list-table::
    :widths: 70 10 10 10
    :header-rows: 1

    * - Algorithm
      - Acronym
      - Status
      - Docs
    * - Orthogonal Matching Pursuit
      - OMP
      - done
      - :ref:`... <api:pursuit:matching>`
    * - Compressive Sampling Matching Pursuit
      - CoSaMP
      - done
      - :ref:`... <api:pursuit:matching>`
    * - Subspace Pursuit
      - CoSaMP
      - done
      - :ref:`... <api:pursuit:matching>`

.. rubric:: Shrinkage and thresholding algorithms

.. list-table::
    :widths: 70 10 10 10
    :header-rows: 1

    * - Algorithm
      - Acronym
      - Status
      - Docs
    * - Iterative Shrinkage Thresholding Algorithm
      - ISTA
      - done
      - :ref:`... <api:sls>`
    * - Fast Iterative Shrinkage Thresholding Algorithm
      - FISTA
      - done
      - :ref:`... <api:sls>`
    * - Iterative Hard Thresholding
      - IHT
      - done
      - :ref:`... <api:pursuit:ht>`
    * - Normalized Iterative Hard Thresholding
      - NIHT
      - done
      - :ref:`... <api:pursuit:ht>`
    * - Hard Thresholding Pursuit
      - HTP
      - done
      - :ref:`... <api:pursuit:ht>`
    * - Normalized Hard Thresholding Pursuit
      - NHTP
      - done
      - :ref:`... <api:pursuit:ht>`
Introduction
=====================

.. contents::
    :depth: 2
    :local:


This library aims to provide XLA/JAX based Python implementations for
various algorithms related to:

* Sparse approximation :cite:`mallat2008wavelet,elad2010sparse`
* Compressive sensing :cite:`donoho2006compressed,candes2006compressive,candes2008introduction,baraniuk2011introduction`
* Linear operators

.. image:: images/srr_cs.png

Bulk of this library is built using functional programming techniques
which is critical for the generation of efficient numerical codes for CPU
and GPU architectures.

Sparse approximation and recovery problems
------------------------------------------------

In the sparse approximation problems :cite:`mallat2008wavelet,elad2010sparse`, we have a 
dictionary of atoms designed for a class of signals
such that the dictionary enables us to construct
a sparse representation of the signal. The
sparse and redundant representation model is:

.. math::

    x = \mathcal{D} \alpha + \eta

where :math:`x \in \mathbb{R}^M` is a single from the given 
class of signals, :math:`\mathcal{D} \in \mathbb{R}^{M \times N}`
is a dictionary consisting of :math:`N` atoms (column vectors) chosen
specifically for the class of signals, :math:`\alpha`
is the sparse representation of :math:`x` in :math:`\mathcal{D}`
giving us an approximation :math:`\hat{x} = \mathcal{D} \alpha`
and :math:`\eta` is the approximation error. The
dictionary :math:`\mathcal{D}` is called the sparsifying dictionary.
The sparse approximation problem consists of finding
the best sparse :math:`\alpha` for a given :math:`x`.

In the compressed sensing (CS) setting,
a sparse signal :math:`x \in \mathbb{R}^N` is captured 
through :math:`M \ll N` linear measurements which are
sufficient to recover :math:`x` from the measurements.
The model is given by:

.. math::

    y = \Phi x + e

where :math:`y \in \mathbb{R}^M` is the vector of :math:`M` linear
measurements on :math:`x`, :math:`\Phi \in \mathbb{R}^{M \times N}` 
is the sensing matrix [or measurement matrix] whose
rows represent the linear functionals on :math:`x`, :math:`x \in \mathbb{R}^N`
is the sparse signal being measured and :math:`e` is the measurement
noise. Typically, :math:`x` by itself is not sparse but it has
a sparse representation in a sparsifying basis :math:`\Psi`
as :math:`x = \Psi \alpha`. The model then becomes:

.. math::

    y = \Phi \Psi \alpha + e.

Sparse recovery consists of finding :math:`\alpha` from
:math:`y` with minimum number of measurements possible.

Both sparse recovery and sparse approximation problems
can be addressed by same algorithms (though their 
performance analysis is different). To simplify the
notation, we will refer to :math:`\mathcal{D}` or :math:`\Phi` 
or :math:`\Phi \Psi` collectively as :math:`A` and attempt to
solve the under-determined system :math:`y = A x + e`
with the prior on the solution that very few entries
in :math:`x` are non-zero. In general, we assume that
:math:`A` is full rank, unless otherwise specified.

The indices of non-zero
entry of :math:`x` form the support of :math:`x`. Corresponding
columns in :math:`A` participate in the sparse
representation of :math:`y`. We can call these columns
also as the support of :math:`x`. 

.. math::

  \mathop{\mathrm{supp}}(x) \triangleq \{i : x_i \neq 0 \}.

Recovering the representation :math:`x`
involves identifying its support :math:`\Lambda = \mathop{\mathrm{supp}}(x)`
and identifying the non-zero entries over the support.
If the support has been 
correctly identified, a straight-forward
way to get the non-zero entries is to compute the
least squares solution :math:`A_{\Lambda}^{\dag} y`.
The :math:`\ell_0` norm of :math:`x` denoted by :math:`\| x\|_0` 
is the number of non-zero entries in :math:`x`.
A representation :math:`y = A x`
is sparse if :math:`\| x\|_0 \ll N`.
An algorithm which can
obtain such  a representation is called a *sparse coding
algorithm*.


:math:`\ell_0` problems
'''''''''''''''''''''''''''''''''

The :math:`K`-SPARSE approximation can be formally expressed as:

.. math::

  \begin{aligned}
    & \underset{x}{\text{minimize}} 
    & &  \| y - A x \|_2 \\
    & \text{subject to}
    & &  \| x \|_0 \leq K.
  \end{aligned}

If the measurements are noiseless, we are interested in 
exact recovery. 
The :math:`K`-EXACT-SPARSE approximation can be formally expressed as:

.. math::

  \begin{aligned}
    & \underset{x}{\text{minimize}} 
    & &  \| x \|_0 \\
    & \text{subject to}
    & &  y = \Phi x\\
    & \text{and}
    & &  \| x \|_0 \leq K
  \end{aligned}


We need to discover both the sparse support for :math:`x` and
the non-zero values over this support. A greedy algorithm
attempts to guess the support incrementally and solves
a smaller (typically least squares) subproblem to estimate
the nonzero values on this support. It then computes the
residual :math:`r = y - A x` and analyzes the correlation of :math:`r`
with the atoms in :math:`A`, via the vector :math:`h = A^T r`, to
improve its guess for the support and update :math:`x` accordingly.


:math:`\ell_1` problems
''''''''''''''''''''''''''''''''


We introduce the different :math:`\ell_1` minimization problems supported by the
``cr.sparse.cvx.admm`` package.

The :math:`\ell_0` problems are not convex. Obtaining a global minimizer 
is not feasible (NP hard). One way around is to use convex relaxation
where a cost function is replaced by its convex version. 
For :math:`\| x \|_0`, the closest convex function is :math:`\| x \|_1` 
or :math:`\ell_1` norm. With this, the exact-sparse recovery problem becomes

.. math::

  {\min}_{x} \| x\|_{1} \; \text{s.t.} \, A x = b


This problem is known as Basis Pursuit (BP) in literature. It can be shown 
that under appropriate conditions on :math:`A`, the basis pursuit solution
coincides with the exact sparse solution. In general, :math:`\ell_1`-norm
minimization problems tend to give sparse solutions.

If :math:`x` is sparse in an sparsifying basis :math:`\Psi` as :math:`x  = \Psi \alpha`
(i.e. :math:`\alpha` is sparse rather than :math:`x`), then we can adapt the
BP formulation as

.. math::

  {\min}_{x} \| W x\|_{1} \; \text{s.t.} \, A x = b

where :math:`W = \Psi^T` and :math:`A` is the sensing matrix :math:`\Phi`.

Finally, in specific problems, different atoms of :math:`\Psi` may 
have different importance. In this case, the :math:`\ell_1` norm
may be adapted to reflect this importance by a non-negative weight vector :math:`w`:

.. math::

  \| \alpha \|_{w,1} = \sum_{i=1}^{N} w_i | \alpha_i |.

This is known as the weighted :math:`\ell_1` semi-norm.

This gives us the general form of the basis pursuit problem

.. math::

  \tag{BP}
  {\min}_{x} \| W x\|_{w,1} \; \text{s.t.} \, A x = b


Usually, the measurement process introduces noise. Thus, 
a constraint :math:`A x = b` is too strict. We can relax this 
to allow for presence of noise as :math:`\| A x - b \|_2 \leq \delta`
where :math:`\delta` is an upper bound on the norm of the measurement noise
or approximation error. 
This gives us the Basis Pursuit with Inequality Constraints (BPIC) problem:

.. math::

  {\min}_{x} \| x\|_{1} \; \text{s.t.} \, \| A x - b \|_2 \leq \delta

The more general form is the L1 minimization problem with L2 constraints:

.. math::

  \tag{L1/L2con}
  {\min}_{x} \| W x\|_{w,1} \; \text{s.t.} \, \| A x - b \|_2 \leq \delta

The constrained BPIC problem can be transformed into an equivalent 
unconstrained convex problem:

.. math::

  {\min}_{x} \| x\|_{1} + \frac{1}{2\rho}\| A x - b \|_2^2.

This is known as Basis Pursuit Denoising (BPDN) in literature.
The more general form is the L1/L2 minimization:

.. math::

  \tag{L1/L2}
  {\min}_{x} \| W x\|_{w,1} + \frac{1}{2\rho}\| A x - b \|_2^2 

We also support corresponding non-negative counter-parts.
The nonnegative basis pursuit problem:

.. math::
  \tag{BP+}
  {\min}_{x} \| W x\|_{w,1} \; \text{s.t.} \, A x = b \, \, \text{and} \, x \succeq 0

The nonnegative L1/L2 minimization or basis pursuit denoising problem:

.. math::

  \tag{L1/L2+}
  {\min}_{x} \| W x\|_{w,1} + \frac{1}{2\rho}\| A x - b \|_2^2  \; \text{s.t.} \, x \succeq 0

The nonnegative L1 minimization problem with L2 constraints:

.. math::

  \tag{L1/L2con+}
  {\min}_{x} \| W x\|_{w,1} \; \text{s.t.} \, \| A x - b \|_2 \leq \delta \, \, \text{and} \, x \succeq 0


Functional Programming
---------------------------


Functional Programming is a programming paradigm where computer programs are constructed 
by applying and composing functions. Functions define a tree of expressions which 
map values to other values (akin to mathematical functions) rather than a sequence
of iterative statements. Some famous languages based on functional programming are
Haskell and Common Lisp.
A key idea in functional programming is a *pure function*. 
A pure function has following properties: 

* The return values are identical for identical arguments.
* The function has no side-effects (no mutation of local static variables, 
  non-local variables, etc.). 


XLA is a domain-specific compiler for linear algebra. 
XLA uses JIT (just-in-time) compilation techniques to analyze the structure of a 
numerical algorithm written using it.
It then specializes the algorithm for actual runtime dimensions and types of parameters involved,
fuses multiple operations together and emits efficient native machine code for
devices like CPUs, GPUs and custom accelerators (like Google TPUs).

JAX is a front-end for XLA and Autograd
with a NumPy inspired API.
Unlike NumPy, JAX arrays are always immutable. While ``x[0] = 10`` is perfectly fine
in NumPy as arrays are mutable, the equivalent functional code in JAX is
``x = x.at[0].set(10)``.


Linear Operators
-----------------------------------------

Efficient linear operator implementations provide much faster
computations compared to direct matrix vector multiplication.
PyLops :cite:`ravasi2019pylops` is a popular collection of
linear operators implemented in Python. 

A linear operator :math:`T : X \to Y` connects a model space :math:`X` 
to a data space :math:`Y`.

A linear operator satisfies following laws:

.. math::

    T (x + y) = T (x) + T (y)

and

.. math::

    T (\alpha x) = \alpha T(x)

Thus, for a general linear combination:

.. math::

    T (\alpha x + \beta y) = \alpha T (x) + \beta T (y)

We are concerned with linear operators :math:`T : \mathbb{F}^n \to \mathbb{F}^m`
where :math:`\mathbb{F}` is either the field of real numbers or 
complex numbers. 
:math:`X = \mathbb{F}^n` is the model space and 
:math:`Y = \mathbb{F}^m` is the data space.
Such a linear operator can be represented by a two dimensional matrix :math:`A`.
The forward operation is given by:

.. math::

    y = A x.

The corresponding adjoint operation is given by:

.. math::

    \hat{x} = A^H y

We represent a linear operator by a pair of functions ``times`` and ``trans``. 
The ``times`` function implements the forward operation while the ``trans``
function implements the adjoint operation.

An inverse problem consists of computing :math:`x` given :math:`y` and :math:`A`.

A framework for building and composing linear operators has been
provided in ``cr.sparse.lop``. Functionality includes:

* Basic operators: identity, matrix, diagonal, zero, flipud, 
  sum, pad_zeros, symmetrize, restriction, etc.
* Signal processing: fourier_basis_1d, dirac_fourier_basis_1d, etc.
* Random dictionaries: gaussian_dict, rademacher_dict, random_onb_dict, random_orthonormal_rows_dict, etc.
* Operator calculus: neg, scale, add, subtract, compose, transpose, hermitian, hcat, etc.
* Additional utilities



Greedy Sparse Recovery/Approximation Algorithms
------------------------------------------------

JAX based implementations for the following algorithms are included.

* Orthogonal Matching Pursuit :cite:`pati1993orthogonal,tropp2004greed`
* Compressive Sampling Matching Pursuit :cite:`needell2009cosamp`
* Subspace Pursuit :cite:`dai2009subspace`
* Iterative Hard Thresholding :cite:`blumensath2009iterative`
* Hard Thresholding Pursuit :cite:`foucart2011recovering`

Convex Optimization based Recovery Algorithms
-----------------------------------------------------

Convex optimization :cite:`boyd2004convex` based methods provide more 
reliable solutions to sparse recovery problems although they tend to be
computationally more complex. 
The first method appeared around 1998 as basis pursuit :cite:`chen1998atomic`.

Alternating directions :cite:`boyd2011distributed` based methods provide
simple yet efficient iterative solutions for sparse recovery. 

:cite:`yang2011alternating` describes inexact ADMM based solutions 
for a variety of :math:`\ell_1` minimization problems. The authors
provide a MATLAB package ``yall1`` :cite:`zhang2010user`. 
A port of ``yall1`` (Your algorithms for :math:`\ell_1`) has been provided.
It provides alternating directions method of multipliers based solutions for
basis pursuit, basis pursuit denoising, basis pursuit with inequality constraints,
their non-negative counterparts and other variants.



Evaluation Framework
--------------------------

The library also provides

* Various simple dictionaries and sensing matrices
* Sample data generation utilities
* Framework for evaluation of sparse recovery algorithms

.. highlight:: shell


Open Source Credits
-----------------------------

Major parts of this library are directly influenced by existing projects.
While the implementation in CR-Sparse is fresh (based on JAX), it has been
possible thanks to the extensive study of existing implementations. We list
here some of the major existing projects which have influenced the implementation
in CR-Sparse. Let us know if we missed anything. 

* `JAX <https://github.com/google/jax>`_ The overall project structure is heavily
  influenced by the conventions followed in JAX. We learned the functional programming
  techniques as applicable for linear algebra work by reading the source code of JAX.
* `SciPy <https://github.com/scipy/scipy>`_ JAX doesn't have all parts of SciPy ported
  yet. Some parts of SciPy have been adapted and re-written (in functional manner) 
  as per the needs of CR-Sparse. E.g. ``cr.sparse.dsp.signals``. The :cite:`torrence1998practical` version
  of CWT in ``cr.sparse.wt``.
* `OpTax <https://github.com/deepmind/optax>`_  This helped in understanding how to 
  use Named Tuples as states for iterative algorithms.  This was also useful 
  in conceptualizing the structure for ``cr.sparse.lop``. 
* `PyLops <https://github.com/PyLops/pylops>`_: The ``cr.sparse.lop`` library is 
  heavily influenced by it.
* `PyWavelets <https://github.com/PyWavelets/pywt>`_: The DWT and CWT implementations
  in ``cr.sparse.wt`` are largely derived from it. The filter coefficients for discrete
  wavelets have been ported from C to Python from here.
* `HTP <https://github.com/foucart/HTP>`_ Original implementation of Hard Thresholding
  Pursuit in MATLAB.
* `WaveLab <https://github.com/gregfreeman/wavelab850>`_ This MATLAB package helped a lot in
  initial understanding of DWT implementation.
* `YALL1 <http://yall1.blogs.rice.edu/>`_: This is the original MATLAB implementation of the
  ADMM based sparse recovery algorithm.
* `L1-LS <https://web.stanford.edu/~boyd/l1_ls/>`_ is the original MATLAB implementation of the
  Truncated Newton Interior Points Method for solving the l1-minimization problem.
* `Sparsify <https://www.southampton.ac.uk/engineering/about/staff/tb1m08.page#software>`_ provides
  the MATLAB implementations of IHT, NIHT, AIHT algorithms.
* `Sparse and Redundant Representations: <https://elad.cs.technion.ac.il/wp-content/uploads/2018/02/Matlab-Package-Book-1.zip>`_ 
  From Theory to Applications in Signal and Image Processing book code helped a lot in basic understanding
  of sparse representations.
* `aaren/wavelets <https://github.com/aaren/wavelets>`_ is a decent CWT implementation following
  :cite:`torrence1998practical`. Influenced: ``cr.sparse.wt``.
  

Further Reading
------------------
* `Functional programming <https://en.wikipedia.org/wiki/Functional_programming>`_
* `How to Think in JAX <https://jax.readthedocs.io/en/latest/notebooks/thinking_in_jax.html>`_
* `JAX - The Sharp Bits <https://jax.readthedocs.io/en/latest/notebooks/Common_Gotchas_in_JAX.html>`_


.. bibliography::
   :filter: docname in docnames


`Documentation <https://carnotresearch.github.io/cr-sparse>`_ | 
`Code <https://github.com/carnotresearch/cr-sparse>`_ | 
`Issues <https://github.com/carnotresearch/cr-sparse/issues>`_ | 
`Discussions <https://github.com/carnotresearch/cr-sparse/discussions>`_ |
`Sparse-Plex <https://sparse-plex.readthedocs.io>`_ 
Continuous Wavelet Transform
==================================


Complex Morlet Wavelets
-------------------------------------

There are several definitions of Complex Morlet Wavelets.

.. math::

    \psi(t) = \frac{1}{\sqrt[4]{\pi}} e^{j \omega_0t } e^{\frac{-t^2}{2}}


Its Fourier transform is:

.. math::

    \Psi(s \omega) = \frac{1}{\sqrt[4]{\pi}} H(\omega) e^{\frac{-(s\omega - \omega_0)^2}{2}}


where :math:`H(\omega)` is the Heaviside step function.

Second definition is more general and is based on two parameters:

- Central frequency: :math:`C`
- Bandwidth: :math:`B` 


.. math::

    \psi(t,B, C) = \frac{1}{\sqrt{\pi B}} \ e^{\frac{-t^2}{B}} \ e^{j2 \pi C t}


This is Gaussian modulated by a complex sinusoid with the standard deviation:

.. math::

    \sigma = \sqrt{\frac{T_p}{2}}

However, this definition doesn't have unit energy.Discrete Wavelet Transform
================================
Filter Banks
====================
Fourier and Wavelet Representations
==========================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   filter_banks
   dwt
   cwt
{{ fullname | escape | underline}}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}


   {% block attributes %}
   {% if attributes %}
   .. rubric:: {{ _('Attributes') }}

   .. autosummary::
   {% for item in attributes %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}Sparse Subspace Clustering
=================================================


.. currentmodule:: cr.sparse.cluster.ssc

SSC-OMP
--------------

.. autosummary::
  :toctree: _autosummary

  build_representation_omp
  build_representation_omp_jit
  batch_build_representation_omp
  batch_build_representation_omp_jit


.. rubric:: Utility functions

.. autosummary::
  :toctree: _autosummary

  sparse_to_full_rep
  sparse_to_bcoo
  bcoo_to_sparse
  bcoo_to_sparse_jit
  rep_to_affinity

Metrics for quality of sparse subspace clustering
------------------------------------------------------


.. autosummary::
  :toctree: _autosummary

  subspace_preservation_stats
  subspace_preservation_stats_jit
  sparse_subspace_preservation_stats
  sparse_subspace_preservation_stats_jit


Tools for analyzing data (with ground truth) 
------------------------------------------------------

.. autosummary::
  :toctree: _autosummary


  angles_between_points
  min_angles_inside_cluster
  min_angles_outside_cluster
  nearest_neighbors_inside_cluster
  nearest_neighbors_outside_cluster
  sorted_neighbors
  inn_positions


Examples
-----------------

* :ref:`gallery:cluster:ssc:omp`



Evaluation Framework
======================================

It is a set of tools to evaluate the performance of sparse recovery algorithms:

* Reconstruction quality of individual sparse recovery problems
* Success rates across multiple sparsity levels

Recovery performance of greedy solvers
-----------------------------------------

.. currentmodule:: cr.sparse

.. autosummary::
    :toctree: _autosummary

    RecoveryPerformance
    RecoveryTrialsAtFixed_M_N

First Order Conic Solvers
=================================

.. contents::
    :depth: 2
    :local:

.. currentmodule:: cr.sparse.focs

This module aims to implement the first order conic solvers 
for sparse signal recovery problems proposed in :cite:`becker2011templates`.
The implementation is adapted from TFOCS :cite:`becker2012tfocs`.

We consider problems of the form:

.. math::

    \text{minimize } \phi(x) = f( \AAA(x) + b) + h(x)

where:

* :math:`\AAA` is a linear operator from :math:`\RR^n \to \RR^m`.
* :math:`b` is a translation vector.
* :math:`f : \RR^m \to \RR` is a *smooth* convex function.
* :math:`h : \RR^n \to \RR` is a *prox-capable* convex function.

For a smooth function, its gradient :math:`g = \nabla f` must exist and be
easy to compute.

For a prox-capable function, it should have an efficient proximal operator:

.. math::

    p_f(x, t) = \text{arg} \min_{z \in \RR^n} f(x) + \frac{1}{2t} \| z - x \|_2^2

for any :math:`x \in \RR^n` and the step size :math:`t > 0`.

See the following sections for details:

* :ref:`api:lop`
* :ref:`api:opt:smooth`
* :ref:`api:opt:proximal`


The routine :py:func:`focs` provides the solver for the minimization
problem described above.

* Unconstrained smooth minimization problems can be handled by choosing
  :math:`h(x) = 0`. See :py:func:`cr.sparse.opt.prox_zero`.
* Convex constraints can be handled by adding their indicator functions 
  as part of :math:`h`.

For solving the minimization problem, an initial solution :math:`x_0` 
should be provided. If one is unsure of the initial solution, they can
provide :math:`0 \in \RR^n` as the initial solution.




Solvers
------------------

.. autosummary::
    :toctree: _autosummary

    focs
    l1rls
    l1rls_jit
    lasso
    lasso_jit
    owl1rls
    owl1rls_jit


Data types
------------------


.. autosummary::
    :toctree: _autosummary
    :nosignatures:
    :template: namedtuple.rst

    FOCSOptions
    FOCSState


Utilities
------------------

.. autosummary::
    :toctree: _autosummary

    matrix_affine_func



In the rest of the document, we will discuss how specific 
sparse signal recovery problems can be modeled and solved
using this module.

L1 regularized least square problem
-----------------------------------------


We consider the problem:

.. math::

    \text{minimize} \frac{1}{2} \| \AAA x - b \|_2^2 + \lambda \| x \|_1 

Choose:

- :math:`f(x) = \frac{1}{2}\| x \|_2^2`  see :func:`cr.sparse.opt.smooth_quad_matrix`
- :math:`h(x) = \| \lambda x \|_1` see :func:`cr.sparse.opt.prox_l1`
- :math:`\AAA` as the linear operator
- :math:`-b` as the translate input

With these choices, it is straight-forward to use :func:`focs` to solve 
the L1RLS problem.  This is implemented in the function :func:`l1rls`.


LASSO
------------------

LASSO (least absolute shrinkage and selection operator) 
s a regression analysis method that performs both variable selection 
and regularization in order to enhance the prediction accuracy 
and interpretability of the resulting statistical model.

We consider the problem: 

.. math::

    \begin{aligned}
    \underset{x}{\text{minimize}} \frac{1}{2} \| \AAA x - b \|_2^2\\
    \text{subject to } \| x \|_1 \leq \tau
    \end{aligned}


Choose:

- :math:`f(x) = \frac{1}{2}\| x \|_2^2`  see :func:`cr.sparse.opt.smooth_quad_matrix`
- :math:`h(x) = I_C(x)` as the indicator function for l1-ball :math:`C = \{x : \| x \|_1 \leq \tau\}`, 
  see :func:`cr.sparse.opt.prox_l1_ball`
- :math:`\AAA` as the linear operator
- :math:`-b` as the translate input

With these choices, it is straight-forward to use :func:`focs` to solve 
the LASSO problem.  This is implemented in the function :func:`lasso`.



Ordered weighted L1 regularized least square problem
-------------------------------------------------------------


We consider the problem:

.. math::

    \underset{x \in \RR^n}{\text{minimize}} \frac{1}{2} \| A x - b \|_2^2 + \sum_{i=1}^n \lambda_i | x |_{(i)} 

described in :cite:`lgorzata2013statistical`.


Choose:

- :math:`f(x) = \frac{1}{2}\| x \|_2^2`  see :func:`cr.sparse.opt.smooth_quad_matrix`
- :math:`h(x) = \sum_{i=1}^n \lambda_i | x |_{(i)}` see :func:`cr.sparse.opt.prox_owl1`
- :math:`\AAA` as the linear operator
- :math:`-b` as the translate input

With these choices, it is straight-forward to use :func:`focs` to solve 
the ordered weighted L1 regularized least square problem.  This is implemented in the function :func:`owl1rls`.
.. _api:opt:

Optimization
==================================

.. contents::
    :depth: 2
    :local:


This section includes several tools which 
form basic building blocks for other higher level
solvers.

We provide generators for:

* Indicator functions
* Projector function
* Proximal operators for functions which are prox capable
* Smooth functions for which gradients can be computed

Most of the descriptions focus on convex sets
:math:`C \subset \RR^n` and convex functions 
with signatures  :math:`\RR^n \to \RR`. However,
the implementations generalize well for the 
complex vector space :math:`\CC^n` with the
real inner product :math:`\langle x, y \rangle = \Re(x^H y)`.


Also, while we describe the math for vectors 
:math:`x \in \RR^n`, the memory layout of these
vectors could be in the form of multi-dimensional arrays too.


.. currentmodule:: cr.sparse.opt

.. _api:opt:indicators:

Indicator Function Generators
------------------------------

Let :math:`C` be a subset of :math:`\RR^n`. An indicator function 
:math:`I_C(x)` is defined as :math:`I_C(x) = 0` for every :math:`x \in C`.
For our purposes, its extended-value extension is relevant which
is defined as:

.. math::

  I_C(x) = \begin{cases} 
    0 & \text{if } x \in C \\
    \infty       & x \notin C
  \end{cases}

An *indicator function generator* creates an indicator function based
on the specification of the set :math:`C`. For example if :math:`C`
is a Euclidean ball, then the specification of :math:`C` is based on
its center and radius.

.. autosummary::
  :toctree: _autosummary
  
  indicator_zero
  indicator_singleton
  indicator_affine
  indicator_box
  indicator_box_affine
  indicator_conic
  indicator_l1_ball
  indicator_l2_ball

.. _api:opt:projectors:

Projection Function Generators
---------------------------------

A key idea in covex optimization is Projection on Convex Sets
(POCS). Given a convex set :math:`C \subset \RR^n` 
and a point :math:`x \in \RR^n`, the projection on to convex 
set :math:`C` is defined as :

.. math::

    P_C(x) =  \text{arg} \min_{v \in C} \| x - v \|_2.

i.e. find the point :math:`v` in :math:`C` which is 
closest to :math:`x`.  If :math:`x` is inside :math:`C`,
then the projection is :math:`x` itself. 

In general, computation of the projection for arbitrary
convex sets can be very hard. However, for specific 
classes of convex sets, the projection can be computed
efficiently. Here we provide projection function generators for
some classes of convex sets.

A *projection function generator* creates a projection function based
on the specification of the set :math:`C`.


.. autosummary::
  :toctree: _autosummary


  proj_identity
  proj_zero
  proj_singleton
  proj_affine
  proj_box
  proj_conic
  proj_l1_ball
  proj_l2_ball


.. _api:opt:smooth:


Smooth Function Generators
----------------------------------

The *smoothness* of a function is a property measured by the
number of continuous derivatives it has over some domain.

We are concerned with convex functions :math:`f : \RR^n \to \RR`
defined over convex sets :math:`C \in \RR^n`
and consider them to be smooth if the gradient :math:`g = \nabla f : \RR^n \to \RR^n`
exists over the interior of its domain. 

For a smooth convex function :math:`f`, we need to be able to

* Compute the value :math:`f(x)` at :math:`x \in \RR^n` using the extended-value extension
  (i.e. :math:`f(x) = \infty \text{if} x \notin C`.
* Compute the gradient :math:`g = \nabla f` at :math:`x`, :math:`g(x)`
* Compute the pair :math:`g(x), f(x)`  at :math:`x`

While JAX provides automatic gradient computation capability, it doesn't 
match exactly at the boundaries of the domain :math:`C` of :math:`f` for many of 
our functions of interest. Hence, we choose to provide our hand-coded 
gradient implementations.

Sometimes, the gradient computation :math:`g(x)`
and function computation :math:`f(x)` have common parts, 
and this can be exploited in improving the algorithm efficiency. 
Hence, we provide an ability to compute the pair too if an
algorithm desires to compute them together efficiently.

We represent smooth functions by the type :py:class:`cr.sparse.opt.SmoothFunction`.

.. autosummary::
    :toctree: _autosummary
    :nosignatures:
    :template: namedtuple.rst

    SmoothFunction

Let `op` be a variable of type :py:class:`cr.sparse.opt.SmoothFunction` 
which represents some smooth function :math:`f`. Then:

* `op.func(x)` returns the function value :math:`f(x)`.
* `op.grad(x)` returns the gradient of function :math:`g(x) = \nabla f(x)`.
* `op.grad_val(x)` returns the pair :math:`(g(x), f(x))`.

A *smooth function generator* creates 
an instance of :py:class:`cr.sparse.opt.SmoothFunction` based
on the specification of the smooth function :math:`f`.

Available smooth function generators:

.. autosummary::
  :toctree: _autosummary

  smooth_constant
  smooth_entropy
  smooth_huber
  smooth_linear
  smooth_entropy
  smooth_logdet
  smooth_quad_matrix

Operations on smooth functions:

.. autosummary::
  :toctree: _autosummary

  smooth_func_translate

We provide tools to build your own smooth functions. 

- You can provide just the definition of the smooth function and 
  we will use JAX to compute the gradient.
- You can provide definitions of the function :math:`f(x)` and the gradient :math:`g(x)`.
- You can provide definitions of :math:`f(x)`, :math:`g(x)` as well as an 
  efficient routine to compute the pair: :math:`(g(x), f(x))`.

.. autosummary::
  :toctree: _autosummary

  smooth_build
  smooth_build2
  smooth_build3

Utilities:

.. autosummary::
  :toctree: _autosummary

  build_grad_val_func
  smooth_value_grad

.. _api:opt:proximal:

Proximal Operator Generators
----------------------------------

The *proximal operator* for a function :math:`f` is defined as

.. math::

    p_f(x, t) = \text{arg} \min_{z \in \RR^n} f(x) + \frac{1}{2t} \| z - x \|_2^2

The proximal operator :math:`p` is a mapping from :math:`(\RR^n, \RR) \to \RR^n`
which maps a vector :math:`x`  to its proximal vector :math:`z = p_f(x,t)`
where :math:`t` can be thought of as a step size for computing the proximal vector.

The proximal operator can be thought of as a generalization of the projection operator. 
If :math:`f` is an indicator function for some convex set :math:`C`, then
the proximal operator is nothing but the projection function onto the set :math:`C`.
If we think of :math:`f` as a cost function over :math:`\RR^n` then indicator 
functions are possibly the toughest cost functions.

For smooth functions, the proximal operator reduces to the gradient step.

We informally call a function as *prox-capable* if its proximal operator 
can be computed efficiently.

For a prox-capable function :math:`f`, in typical proximal algorithms, we need to be able to

* Compute the value :math:`f(x)` at :math:`x \in \RR^n` using the extended-value extension
* Compute the proximal vector for the step size :math:`t` as :math:`p_f(x, t)`
* Compute the pair :math:`p_f(x, t), f(p_f(x, t))`

Note that in the last case, we first compute
:math:`z = p_f(x, t)` and then compute the value :math:`v  = f(z)`.


We represent *prox-capable* functions by the type :py:class:`cr.sparse.opt.ProxCapable`.

.. autosummary::
    :toctree: _autosummary
    :nosignatures:
    :template: namedtuple.rst

    ProxCapable

Let `op` be a variable of type :py:class:`cr.sparse.opt.ProxCapable` 
which represents some prox-capable function :math:`f`. Then:

* `op.func(x)` returns the function value :math:`f(x)`.
* `op.prox_op(x)` returns the proximal vector for a step size: :math:`z = p_f(x, t)`.
* `op.prox_vec_val(x)` returns the pair :math:`z,v = p_f(x, t), f(z)`.

A *proximal operator generator* creates 
an instance of :py:class:`cr.sparse.opt.ProxCapable` based
on the specification of the prox-capable function :math:`f`.

Available proximal operator generators:

.. autosummary::
  :toctree: _autosummary

  prox_zero
  prox_l1
  prox_l2
  prox_l1_pos
  prox_l1_ball
  prox_owl1

You can build your own :py:class:`cr.sparse.opt.ProxCapable` wrappers
by providing the definition of the function :math:`f(x)` and its 
proximal operator :math:`p_f(x, t)`.

.. autosummary::
  :toctree: _autosummary
  
  prox_build
  build_from_ind_proj





Simpler Projection Functions
----------------------------------

.. autosummary::
  :toctree: _autosummary

  project_to_ball
  project_to_box
  project_to_real_upper_limit



Shrinkage
---------------------

.. autosummary::
  :toctree: _autosummary

  shrink


Conjugate Gradient Methods
----------------------------------

.. rubric:: Normal Conjugate Gradients on Matrices

.. autosummary::
  :toctree: _autosummary

    cg.solve_from
    cg.solve_from_jit
    cg.solve
    cg.solve_jit


.. rubric:: Preconditioned Normal Conjugate Gradients on Linear Operators

These are more general purpose.

.. autosummary::
  :toctree: _autosummary

    pcg.solve_from
    pcg.solve_from_jit
    pcg.solve
    pcg.solve_jit

.. _api:pursuit:

Greedy Sparse Recovery
================================================

.. contents::
    :depth: 2
    :local:

.. rubric:: Algorithm versions

Several algorithms are available in multiple versions.

The library allows a dictionary or a sensing process to be represented as either:

* A matrix of real/complex values
* A linear operator with fast implementation (see ``cr.sparse.lop`` module)

E.g. a partial sensing sensing matrix can be implemented more efficiently using a linear operator consisting of 
applying the fourier transform followed by selecting a subset of fourier measurements. 

* A prefix ``matrix_`` indicates an implementation which accepts matrices as dictionaries or compressive sensors.
* A prefix ``operator_`` indicates an implementation which accepts linear operators described in ``cr.sparse.lop`` module as dictionaries or compressive sensors.

* A suffix ``_jit`` means that it is the JIT (Just In Time ) compiled version of the original implementation.
* A suffix ``_multi`` means that it is the version of implementation which can process multiple signals/measurement vectors 
  simultaneously. The recovery problem :math:`y = \Phi x + e` is extended to :math:`Y = \Phi X + E` such that:

  * Each column of :math:`Y` represents one signal/measurement vector
  * Each column of :math:`X` represents one representation vector to be recovered
  * Each column of :math:`E` representation corresponding measurement error/noise.


.. rubric:: Conditions on dictionaries/sensing matrices

Different algorithms have different requirements on the dictionaries or sensing matrices:

* Some algorithms accept overcomplete dictionaries/sensing matrices with unit norm columns 
* Some algorithms accept overcomplete dictionaries/sensing matrices with orthogonal rows
* Some algorithms accept any overcomplete dictionary



.. currentmodule:: cr.sparse.pursuit

.. _api:pursuit:matching:

Basic Matching Pursuit Based Algorithms
------------------------------------------

.. rubric:: Orthogonal Matching Pursuit

.. autosummary::
  :toctree: _autosummary

    omp.solve
    omp.matrix_solve
    omp.matrix_solve_jit
    omp.matrix_solve_multi

Compressive Sensing Matching Pursuit (CSMP) Algorithms

.. rubric:: Compressive Sampling Matching Pursuit

.. autosummary::
  :toctree: _autosummary

    cosamp.solve
    cosamp.matrix_solve
    cosamp.matrix_solve_jit
    cosamp.operator_solve
    cosamp.operator_solve_jit

.. rubric:: Subspace Pursuit

.. autosummary::
  :toctree: _autosummary

    sp.solve
    sp.matrix_solve
    sp.matrix_solve_jit
    sp.operator_solve
    sp.operator_solve_jit


.. _api:pursuit:ht:

Hard Thresholding Based Algorithms
-----------------------------------------

.. rubric:: Iterative Hard Thresholding

.. autosummary::
  :toctree: _autosummary

    iht.solve
    iht.matrix_solve
    iht.matrix_solve_jit
    iht.operator_solve
    iht.operator_solve_jit

.. rubric:: Hard Thresholding Pursuit

.. autosummary::
  :toctree: _autosummary

    htp.solve
    htp.matrix_solve
    htp.matrix_solve_jit
    htp.operator_solve
    htp.operator_solve_jit


Data Types
-------------------------------


.. autosummary::
  :nosignatures:
  :toctree: _autosummary
  :template: namedtuple.rst

    RecoverySolution

Utilities
-------------------------------

.. autosummary::
  :toctree: _autosummary

    abs_max_idx
    gram_chol_update



Using the greedy algorithms
-------------------------------

These algorithms solve the inverse problem :math:`y = \Phi x + e` where :math:`\Phi` and :math:`y`
are known and :math:`x` is desired with :math:`\Phi` being an overcomplete dictionary or a sensing matrix.

For sparse approximation problems, we require the following to invoke any of these algorithms:

* A sparsifying dictionary :math:`\Phi`.
* A signal :math:`y` which is expected to have a sparse or compressible representation :math:`x` in :math:`\Phi`.

For sparse recovery problems, we require the following to invoke any of these algorithms:

* A sensing matrix :math:`\Phi` with suitable RIP or other properties.
* A measurement vector :math:`y` generated by applying :math:`\Phi` to a sparse signal :math:`x`

.. rubric:: A synthetic example

Build a Gaussian dictionary/sensing matrix::

  from jax import random
  import cr.sparse.dict as crdict
  M = 128
  N = 256
  key = random.PRNGKey(0)
  Phi = crdict.gaussian_mtx(key, M,N)

Build a K-sparse signal with Gaussian non-zero entries::

  import cr.sparse.data as crdata
  import jax.numpy as jnp
  K = 16
  key, subkey = random.split(key)
  x, omega = crdata.sparse_normal_representations(key, N, K, 1)
  x = jnp.squeeze(x)

Build the measurement vector::

  y = Phi @ x

We have built the necessary inputs for a sparse recovery problem. It is time to run the solver.

Import a sparse recovery solver::

  from cr.sparse.pursuit import cosamp

Solve the recovery problem::

  solution =  cosamp.matrix_solve(Phi, y, K)

You can choose any other solver. 

The support for the non-zero entries in the solution is given by ``solution.I`` and
the values for non-zero entries are given by ``solution.x_I``. You can build the
sparse representation as follows::

  from cr.sparse import build_signal_from_indices_and_values
  x_hat = build_signal_from_indices_and_values(N, solution.I, solution.x_I)

Finally, you can use the utility to evaluate the quality of reconstruction::

  from cr.sparse.ef import RecoveryPerformance
  rp = RecoveryPerformance(Phi, y, x, x_hat)
  rp.print()

This would output something like::

  M: 128, N: 256, K: 16
  x_norm: 3.817, y_norm: 3.922
  x_hat_norm: 3.817, h_norm: 1.55e-06, r_norm: 1.72e-06
  recovery_snr: 127.83 dB, measurement_snr: 127.16 dB
  T0: [ 27  63  79  85  88 111 112 124 131 137 160 200 230 234 235 250]
  R0: [ 27  63  79  85  88 111 112 124 131 137 160 200 230 234 235 250]
  Overlap: [ 27  63  79  85  88 111 112 124 131 137 160 200 230 234 235 250], Correct: 16
  success: True

.. _api:lop:

Linear Operators
=======================

.. contents::
    :depth: 2
    :local:

We provide a collection of linear operators with efficient JAX based implementations 
that are
relevant in standard signal/image processing problems.
We also provide a bunch of utilities to combine and convert
linear operators.

This module is inspired by ``pylops`` although the implementation approach 
is different.

A linear operator :math:`T : X \to Y` connects a model space :math:`X` 
to a data space :math:`Y`.

A linear operator satisfies following laws:

.. math::

    T (x + y) = T (x) + T (y)


.. math::

    T (\alpha x) = \alpha T(x)


Thus, for a general linear combination:

.. math::

    T (\alpha x + \beta y) = \alpha T (x) + \beta T (y)

We are concerned with linear operators :math:`T : \mathbb{F}^n \to \mathbb{F}^m`
where :math:`\mathbb{F}` is either the field of real numbers or 
complex numbers. 
:math:`X = \mathbb{F}^n` is the model space and 
:math:`Y = \mathbb{F}^m` is the data space.

Such a linear operator can be represented by a two dimensional matrix :math:`A`.

The forward operation is given by:

.. math::

    y = A x.

The corresponding adjoint operation is given by:

.. math::

    \hat{x} = A^H y

We represent a linear operator by a pair of functions ``times`` and ``trans``. 
The ``times`` function implements the forward operation while the ``trans``
function implements the adjoint operation.

An inverse problem consists of computing :math:`x` given :math:`y` and :math:`A`.

..  rubric:: 1D, 2D, ND operators

* A 1D operator takes a 1D array as input and returns a 1D array as output.
  E.g. identity, pad, etc.
* A 2D operator takes a 2D array as input and returns a 2D array as output. 
  E.g. conv2D, dwt2D, etc.
* An ND operator takes an ND array as input and returns an ND array as output.

The vectors may be stored using multi-dimensional arrays in memory. 
E.g., images are usually stored in 2-3 dimensional arrays.
The operators themselves may work directly on multi-dimensions arrays. 
E.g. a 2D convolution operator can be applied directly
to an image to result in another image. 

In other words, the vectors from the model space as well as data space 
may be stored in memory using 1D,2D,...,ND array.  They should still
be treated as vectors for the purposes of this module.

This is a departure from pylops convention where input to a 
linear operator must be flattened into a 1D array and output 
needs to be reshaped again.
In this library, the input and output to a 2D linear operator 
would be a 2D array.

.. rubric:: axis parameter in a 1D linear operator

* A 1D linear operator may get an ND array as input.
* In this case, the axis parameter to the operator specifies the
  axis along which the linear operator is to be applied.
* The input ND array will be broken into slices of 1D arrays along the
  axis and the linear operator will be applied separately to each slice.
* Then the slices will be combined to generate the output ND array.
* E.g. if the input is a matrix then:
  
  * axis=0 means apply the linear operator over each column (along axis=0)
  * axis=1 means apply the linear operator over each row  (along axis=1)

This is based on the convention followed by ``numpy.apply_along_axis``.

.. currentmodule:: cr.sparse.lop


Data types
------------------


.. autosummary::
    :toctree: _autosummary
    :nosignatures:
    :template: namedtuple.rst

    Operator

Basic operators
------------------

.. autosummary::
    :toctree: _autosummary

    identity
    matrix
    real_matrix
    scalar_mult
    diagonal
    zero
    flipud
    sum
    dot
    pad_zeros
    symmetrize
    restriction
    reshape
    arr2vec

Operator calculus
------------------

It is possible to combine one or more linear operators
to create new linear operators. The functions in this
section provide different ways to combine linear operators.

.. autosummary::
    :toctree: _autosummary

    neg
    scale
    add
    subtract
    compose
    transpose
    adjoint
    hcat
    power
    block_diag

Signal processing operators
------------------------------------

.. autosummary::
    :toctree: _autosummary

    running_average
    fir_filter
    convolve
    convolve2D
    convolveND


Orthonormal transforms and bases
------------------------------------------------

.. autosummary::
    :toctree: _autosummary

    fft
    dwt
    dwt2D
    fourier_basis
    cosine_basis
    walsh_hadamard_basis

Unions of bases
--------------------

.. autosummary::
    :toctree: _autosummary

    dirac_fourier_basis


Random compressive sensing operators
--------------------------------------

.. autosummary::
    :toctree: _autosummary

    gaussian_dict
    rademacher_dict
    random_onb_dict
    random_orthonormal_rows_dict


Operators for special matrices
------------------------------------

.. autosummary::
    :toctree: _autosummary

    circulant


Derivatives (finite differences)
--------------------------------------

.. autosummary::
    :toctree: _autosummary

    first_derivative
    second_derivative
    tv
    tv2D

Convenience operators
-----------------------

These operators are technically not linear on :math:`\mathbb{F}^n \to \mathbb{F}^m`

.. autosummary::
    :toctree: _autosummary

    real


Operator parts
------------------

.. autosummary::
    :toctree: _autosummary

    column
    columns

Properties of a linear operator
--------------------------------------------

These are still experimental and not efficient.

.. autosummary::
    :toctree: _autosummary

    normest
    normest_jit
    upper_frame_bound


Utilities
-------------------

.. autosummary::
    :toctree: _autosummary

    jit
    to_matrix
    to_adjoint_matrix
    to_complex_matrix
.. _api:data:

Sample Data Generation Utilities
=====================================


Sparse Model Vectors
------------------------

.. currentmodule:: cr.sparse.data

.. autosummary::
  :toctree: _autosummary

  sparse_normal_representations
  sparse_spikes



Subspaces
------------------------

.. autosummary::
  :toctree: _autosummary

  random_subspaces
  uniform_points_on_subspaces
  two_subspaces_at_angle
  three_subspaces_at_angle
  
Utilities
==============================

.. contents::
    :depth: 2
    :local:


.. currentmodule:: cr.sparse



Sparse representations
------------------------------------

Following functions analyze or construct representation vectors which are known to be sparse.

.. autosummary::
  :toctree: _autosummary

    nonzero_values
    nonzero_indices
    support
    randomize_rows
    randomize_cols
    largest_indices
    largest_indices_by
    hard_threshold
    hard_threshold_sorted
    hard_threshold_by
    sparse_approximation
    build_signal_from_indices_and_values
    dynamic_range
    nonzero_dynamic_range


.. rubric:: Sparse representation matrices (row-wise)

.. autosummary::
  :toctree: _autosummary

    largest_indices_rw
    take_along_rows
    sparse_approximation_rw

.. rubric:: Sparse representation matrices (column-wise)

.. autosummary::
  :toctree: _autosummary

    largest_indices_cw
    take_along_cols
    sparse_approximation_cw




Basic Signal Information
---------------------------------------------

.. autosummary::
  :toctree: _autosummary

  frequency_spectrum
  power_spectrum
  energy

Basic Signal Processing
-------------------------------

.. autosummary::
  :toctree: _autosummary

  normalize
  interpft


Artificial Noise
-----------------------------------


.. autosummary::
  :toctree: _autosummary

  awgn_at_snr

.. _api:sls:

Sparse Linear Systems
=========================

The solvers in this module focus on 
traditional least square problems
for square or overdetermined linear systems
:math:`A x = b` 
where the matrix :math:`A` is sparse and is 
represented by a linear operator abstraction 
providing the matrix multiplication and adjoint 
multiplication functions.

.. currentmodule:: cr.sparse.sls

Solvers
-----------------

.. autosummary::
    :toctree: _autosummary

    lsqr
    lsqr_jit
    power_iterations
    power_iterations_jit
    ista
    ista_jit
    fista
    fista_jit

Data types
------------------


.. autosummary::
    :toctree: _autosummary
    :nosignatures:
    :template: namedtuple.rst

    LSQRSolution
    PowerIterSolution
    ISTAState
    FISTAState


.. _api:dict:

Sparsifying Dictionaries and Sensing Matrices
===================================================


.. currentmodule:: cr.sparse.dict


Functions for constructing sparsying dictionaries and sensing matrices
--------------------------------------------------------------------------


.. autosummary::
  :toctree: _autosummary


    gaussian_mtx
    rademacher_mtx
    random_onb
    hadamard
    hadamard_basis
    dirac_hadamard_basis
    cosine_basis
    dirac_cosine_basis
    dirac_hadamard_cosine_basis
    fourier_basis
    random_orthonormal_rows



Dictionary properties
-------------------------


.. autosummary::
  :toctree: _autosummary

    gram
    frame
    coherence_with_index
    coherence
    frame_bounds
    upper_frame_bound
    lower_frame_bound
    babel



Dictionary comparison
----------------------------

These functions are useful for comparing dictionaries
during the dictionary learning process. 

.. autosummary::
  :toctree: _autosummary

    mutual_coherence_with_index
    mutual_coherence
    matching_atoms_ratio.. _api:cs:

Compressive Sensing
===============================

.. contents::
    :depth: 2
    :local:


This package contains functions which are specific 
to Compressive Sensing applications.



1-Bit Compressive Sensing
----------------------------

.. currentmodule:: cr.sparse.cs.cs1bit


.. autosummary::
    :nosignatures:
    :toctree: _autosummary

    measure_1bit
    biht
    biht_jit

.. autosummary::
  :nosignatures:
  :toctree: _autosummary
  :template: namedtuple.rst

    BIHTState
.. _api:dsp:

Digital Signal Processing
===============================

.. contents::
    :depth: 2
    :local:

The ``CR-Sparse`` library has some handy digital signal processing routines
implemented in JAX. They are available as part of the ``cr.sparse.dsp``
package.


Utilities
-----------------------

.. currentmodule:: cr.sparse.dsp

.. autosummary::
    :nosignatures:
    :toctree: _autosummary

    time_values

Synthetic Signals
-----------------------

.. currentmodule:: cr.sparse.dsp.signals

.. autosummary::
    :nosignatures:
    :toctree: _autosummary

    pulse
    transient_sine_wave
    decaying_sine_wave
    chirp
    chirp_centered
    gaussian_pulse


.. currentmodule:: cr.sparse.dsp



Discrete Cosine Transform
-------------------------------------

.. autosummary::
    :nosignatures:
    :toctree: _autosummary

    dct
    idct
    orthonormal_dct
    orthonormal_idct

.. currentmodule:: cr.sparse.dsp

Fast Walsh Hadamard Transform
------------------------------

There is no separate Inverse Fast Walsh Hadamard Transform as FWHT is the inverse of
itself except for a normalization factor.
In other words,  ``x == fwht(fwht(x)) / n`` where n is the length of x.

.. autosummary::
    :nosignatures:
    :toctree: _autosummary

    fwht

.. _api:wavelets:

Wavelets
=====================

.. contents::
    :depth: 2
    :local:


``CR-Sparse`` provides support for both DWT (Discrete Wavelet Transform)
and CWT (Continuous Wavelet Transform).


The support for discrete wavelets is a partial port of 
wavelets functionality from the `PyWavelets <https://pywavelets.readthedocs.io/en/latest/>`_ project. 
The functionality has been written on top of JAX. While PyWavelets gets its performance from the C extensions in its 
implementation, we have built the functionality on top of JAX API. 

- API and implementation are both based on functional programming paradigm.
- There are no C extensions. All implementation is pure Python.
- The implementation takes advantage of XLA and can run easily on GPUs and TPUs.

Continuous Wavelet Transform have been implemented following :cite:`torrence1998practical`.
A reference implementation using NumPy is 
`here <https://github.com/aaren/wavelets>`_. 


The code examples in this section will assume following imports::

  import cr.sparse as crs
  import cr.sparse.wt as wt


Discrete Wavelets
--------------------------------------


API is available at two levels

- Functions which directly correspond to the high level API of PyWavelets.
- Lower level functions which are JIT compiled.

The high level functions involve handling a variety of use cases for
the arguments passed. For example, they can accept lists as well as
JAX nd-arrays. These functions cannot be JIT compiled. Lower level 
functions have been carefully designed to accept arguments which 
fit the JIT rules of JAX. They can be embedded in another JIT 
compiled function.

Current support is focused on discrete wavelet transforms. 
Following wavelets are supported.

    bior1.1 bior1.3 bior1.5 bior2.2 bior2.4 bior2.6 bior2.8 bior3.1 bior3.3 bior3.5 bior3.7 bior3.9 bior4.4 bior5.5 bior6.8 coif1 coif2 coif3 coif4 coif5 coif6 coif7 coif8 coif9 coif10 coif11 coif12 coif13 coif14 coif15 coif16 coif17 db1 db2 db3 db4 db5 db6 db7 db8 db9 db10 db11 db12 db13 db14 db15 db16 db17 db18 db19 db20 db21 db22 db23 db24 db25 db26 db27 db28 db29 db30 db31 db32 db33 db34 db35 db36 db37 db38 dmey haar rbio1.1 rbio1.3 rbio1.5 rbio2.2 rbio2.4 rbio2.6 rbio2.8 rbio3.1 rbio3.3 rbio3.5 rbio3.7 rbio3.9 rbio4.4 rbio5.5 rbio6.8 sym2 sym3 sym4 sym5 sym6 sym7 sym8 sym9 sym10 sym11 sym12 sym13 sym14 sym15 sym16 sym17 sym18 sym19 sym20


.. currentmodule:: cr.sparse.wt


High-level API
'''''''''''''''''''''''''''''''

.. rubric:: Data types

.. autosummary::
  :nosignatures:
  :toctree: _autosummary
  :template: namedtuple.rst

    FAMILY
    SYMMETRY

.. autosummary::
  :nosignatures:
  :toctree: _autosummary

    DiscreteWavelet

.. rubric:: Wavelets

.. autosummary::
  :toctree: _autosummary

  families
  build_wavelet
  wavelist
  is_discrete_wavelet
  wname_to_family_order
  build_discrete_wavelet



.. rubric:: Wavelet transforms

.. autosummary::
  :toctree: _autosummary

    dwt
    idwt
    dwt2
    idwt2
    downcoef
    upcoef
    wavedec
    waverec
    dwt_axis
    idwt_axis
    dwt_column
    dwt_row
    dwt_tube
    idwt_column
    idwt_row
    idwt_tube

.. rubric:: Utilities

.. autosummary::
  :toctree: _autosummary

    modes
    pad
    dwt_max_level
    dwt_coeff_len
    up_sample

Lower-level API
'''''''''''''''''''''''''''''


.. autosummary::
  :toctree: _autosummary

  dwt_
  idwt_ 
  downcoef_
  upcoef_
  dwt_axis_
  idwt_axis_


.. _ref-wt-modes:

Signal Extension Modes
''''''''''''''''''''''''''''''''

Real world signals are finite. They are typically stored in 
finite size arrays in computers. Computing the wavelet transform
of signal values around the boundary of the signal inevitably involves
assuming some form of signal extrapolation. A simple extrapolation
method is to extend the signal with zeros at the boundary. 
Reconstruction of the signal from its wavelet coefficients may introduce
boundary artifacts based on how the signal was extrapolated. A careful
choice of signal extension method is necessary based on actual 
application. 

We provide following signal extension modes at the moment.

zero
  Signal is extended by adding zeros::

    >>> wt.pad(jnp.array([1,2,4,-1,2,-1]), 2, 'zero')
    DeviceArray([ 0,  0,  1,  2,  4, -1,  2, -1,  0,  0], dtype=int64)


constant
  Border values of the signal are replicated::

    >>> wt.pad(jnp.array([1,2,4,-1,2,-1]), 2, 'constant')
    DeviceArray([ 1,  1,  1,  2,  4, -1,  2, -1, -1, -1], dtype=int64)


symmetric
  Signal is extended by mirroring the samples at the border in mirror form. 
  The border sample is also mirrored.::

    >>> wt.pad(jnp.array([1,2,4,-1,2,-1]), 2, 'symmetric')
    DeviceArray([ 2,  1,  1,  2,  4, -1,  2, -1, -1,  2], dtype=int64)


reflect
  Signal is extended by reflecting the samples around the border sample.
  Border sample is not copied in the extension.:: 

    >>> wt.pad(jnp.array([1,2,4,-1,2,-1]), 2, 'reflect')
    DeviceArray([ 4,  2,  1,  2,  4, -1,  2, -1,  2, -1], dtype=int64)

periodic
  Signal is extended periodically. The samples at the end repeat at the extension
  at the beginning. The samples at the beginning repeat at the extension at the end.::

    >>> wt.pad(jnp.array([1,2,4,-1,2,-1]), 2, 'periodic')
    DeviceArray([ 2, -1,  1,  2,  4, -1,  2, -1,  1,  2], dtype=int64)

periodization
  The signal is extended the same way as the periodic extension. The major difference is that
  the number of wavelet coefficients is identical to the length of the signal. All extra values
  are trimmed.

Many of the signal extension modes are similar to the padding modes supported by the
``jax.numpy.pad`` function. However, the naming convention is different and follows 
PyWavelets.


Continuous Wavelets
-----------------------------------


.. rubric:: Further Reading

.. bibliography::
   :filter: docname in docnamesAvailable Modules in CR-Sparse
====================================

.. contents::
    :depth: 2
    :local:

cr.sparse
-----------------

.. automodule:: cr.sparse
    :members:

cr.sparse.dict
-----------------

.. automodule:: cr.sparse.dict
    :members:


cr.sparse.data
-----------------

.. automodule:: cr.sparse.data
    :members:

cr.sparse.lop
-----------------

.. automodule:: cr.sparse.lop
    :members:

cr.sparse.wt
-----------------

.. automodule:: cr.sparse.wt
    :members:

cr.sparse.pursuit
-------------------

.. automodule:: cr.sparse.pursuit
    :members:
API Docs
=============================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   dsp
   wavelets
   lop
   sls
   dict
   pursuit
   cvx_recovery
   opt
   focs
   cs
   cluster
   ssc
   data
   util
   geo
   vision
   ef

.. ::

   modules



.. _api:geo:

Geophysical Signal Processing
===================================

This module contains some utilities for processing of 
geophysical signals. 

.. currentmodule:: cr.sparse.geo


Miscalleneous
-----------------

.. autosummary::
    :toctree: _autosummary

    ricker
.. _api:vision:

Computer Vision and Image Processing 
======================================

This module contains some utilities for image processing
and computer vision. 


.. note::

    CR-Sparse doesn't intend to be a full
    fledged image processing and computer vision solution.
    The utilities in this module are meant to complement 
    the functionality available in CR-Sparse to demonstrate
    how they can be used on images.


.. currentmodule:: cr.sparse.vision

Kernels
-----------------


.. autosummary::
    :toctree: _autosummary

    kernel_gaussian
    KERNEL_BOX_3X3
    KERNEL_BOX_5X5
    KERNEL_BOX_7X7
    KERNEL_BOX_21X21
    KERNEL_SHARPEN_3X3
    KERNEL_LAPLACIAN_3X3
    KERNEL_SOBEL_X
    KERNEL_SOBEL_Y

.. _api:l1min:

L1 Minimization
=============================================================================

.. contents::
    :depth: 2
    :local:

.. _api:l1min:tnipm:

Truncated Newton Interior Points Method (TNIPM)
---------------------------------------------------------

This method can be used to solve problems of type:

.. math::

  \min \| A x  - b \|_2^2 + \lambda \| x \|_1 

The method works as follows:

* The :math:`\ell_1` regularized LSP (Least Squares Program) above is 
  transformed into a convex quadratic program with linear inequality constraints.
* The logarithmic barrier for the quadratic program is constructed.
* The central path for this logarithmic barrier minimizer is identified.
* Truncated newton method is used to solve for points on the central path.

  * Compute the search direction as an approximate solution to the Newton system
  * Compute the step size by backtracking line search
  * Update the iterate
  * Construct a dual feasible point 
  * Evaluate the duality gap
  * Repeat till the relative duality gap is within a specified threshold

* The Newton system is solved approximately using preconditioned conjugate gradients.

The ``solve_from`` versions below are useful if the solution is known partially.
The ``solve`` versions are more directly applicable when solution is not known.
Both ``solve_from`` and ``solve`` are available as regular and JIT-compiled 
versions.


.. currentmodule:: cr.sparse.cvx

.. autosummary::
  :toctree: _autosummary

    l1ls.solve_from
    l1ls.solve_from_jit
    l1ls.solve
    l1ls.solve_jit

An example is provided in :ref:`recovery_l1_gallery`.


.. _api:l1min:admmm:

Alternating Directions Methods
-------------------------------------

This is based on :cite:`yang2011alternating`.

A tutorial has been provided to explore these 
methods in action. 
The ``yall1.solve`` method is an overall wrapper method 
for solving different types of :math:`\ell_1` minimization
problems. It in turn calls the lower level methods for solving
specific types of problems.

.. currentmodule:: cr.sparse.cvx.adm

.. autosummary::
  :toctree: _autosummary

    yall1.solve
    yall1.solve_bp
    yall1.solve_bp_jit
    yall1.solve_l1_l2
    yall1.solve_l1_l2_jit
    yall1.solve_l1_l2con
    yall1.solve_l1_l2con_jit
.. _api:cluster:

Data Clustering
==========================


Vector Quantization
------------------------

.. currentmodule:: cr.sparse.cluster.vq

.. autosummary::
  :toctree: _autosummary

  kmeans
  kmeans_jit
  kmeans_with_seed
  kmeans_with_seed_jit
  find_nearest
  find_nearest_jit
  find_assignment
  find_assignment_jit
  find_new_centroids
  find_new_centroids_jit
  
 
Spectral Clustering
------------------------

.. currentmodule:: cr.sparse.cluster.spectral

.. autosummary::
  :toctree: _autosummary

  unnormalized
  unnormalized_k
  unnormalized_k_jit
  normalized_random_walk


Data types
--------------------

.. currentmodule:: cr.sparse.cluster

.. autosummary::
  :nosignatures:
  :toctree: _autosummary
  :template: namedtuple.rst

  vq.KMeansState
  vq.KMeansSolution
  spectral.SpectralclusteringSolution

.. _sls:ista:

Iterative Shrinkage and Thresholding Algorithm
======================================================

.. contents::
    :depth: 2
    :local:

Our description is based on :cite:`daubechies2004iterative, elad2010sparse, zibulevsky2010l1`.

ISTA can be used to solve problems of the form:

.. math::
    :label: ista_l1_l2_minimization

    \widehat{x} = \text{arg} \min_{x} \frac{1}{2}\| b - A x \|_2^2 + \lambda \| x \|_1 

and its variants (with different regularizations).

Here the objective function is:

.. math::

    f(x) = \frac{1}{2}\| b - A x \|_2^2 + \lambda \| x \|_1

Derivation
------------------

For an identity :math:`A`, the problem reduces to:

.. math::

    \widehat{x} = \text{arg} \min_{x} \frac{1}{2}\| b - x \|_2^2 + \lambda \| x \|_1 

The term on the R.H.S. : 

.. math::

    f(x) = \frac{1}{2}\| b - x \|_2^2 + \lambda \| x \|_1  = \sum_{i=1}^n \left [\frac{1}{2} | b_i - x_i|^2  + \lambda |x_i| \right ]

is separable in the components of :math:`x`.

The scalar function

.. math::

    g(\tau) = \frac{1}{2} | \gamma - \tau|^2  + \lambda |\tau|

has a minimizer given by the soft thresholding function: 

.. math::

    \tau^* = \mathbf{ST}_{\lambda}(\gamma) = \begin{cases} 
     \gamma\left ( 1 -  \frac{\lambda}{|\gamma|} \right ) & \text{for} & |\gamma| \gt \lambda\\
     0 & \text{for} & |\gamma| \le \lambda
    \end{cases}

:cite:`zibulevsky2010l1` show that a similar solution emerges when :math:`A` is unitary too.

:cite:`daubechies2004iterative` introduced a surrogate term:

.. math::

    \text{dist}(x, x_0) = \frac{c}{2} \| x  - x_0 \|_2^2 - \frac{1}{2} \| A x - A x_0 \|_2^2 

For this function to be strictly convex w.r.t. :math:`x`, its Hessian must be positive definite. 
This holds if:

.. math::

    c > \| A^H A \|_2 = \lambda_{\max}(A^H A)

the largest eigen value :math:`A^H A`.

The new surrogate objective function is defined as:

.. math::

    \widetilde{f}(x, x_0) = f(x) + \text{dist}(x, x_0) 
    = \frac{1}{2}\| b - A x \|_2^2 + \lambda \| x \|_1 + \frac{c}{2} \| x  - x_0 \|_2^2 - \frac{1}{2} \| A x - A x_0 \|_2^2

Introducing:

.. math::

    \nu_0 = x_0 + \frac{1}{c} A^H (b - A x_0)

We can rewrite the surrogate objective as:

.. math::

    \widetilde{f}(x, x_0) = \text{Const} + \lambda \| x \|_1 + \frac{c}{2} \| x - \nu_0 \|_2^2

As discussed above, this surrogate objective can be minimized using the soft thresholding operator:

.. math::

    x^* = \mathbf{ST}_{\frac{\lambda}{c} }(\nu_0) = \mathbf{ST}_{\frac{\lambda}{c} }\left (x_0 + \frac{1}{c} A^H (b - A x_0) \right )


Starting with some :math:`x_0`, we can propose an itarative method over a sequence 
:math:`x_i = \{x_0, x_1, x_2, \dots \}` as:

.. math::
    :label: ista_iteration_soft_thresholding

    x_{i+1} = \mathbf{ST}_{\frac{\lambda}{c} }\left (x_i + \frac{1}{c} A^H (b - A x_i)\right )

This is the IST algorithm.

By changing the regularization in :eq:`ista_l1_l2_minimization`, we can derive different IST algorithms with different thresholding 
functions. The version below considers a generalized thresholding function which depends on the regularizer.

.. math::
    :label: ista_iteration_thresholding

    x_{i+1} = \mathbf{T}_{\frac{\lambda}{c} }\left (x_i + \frac{1}{c} A^H (b - A x_i)\right )

Sparsifying Basis
'''''''''''''''''''''''

Often, the signal :math:`x` (e.g. an image) may not be sparse or compressible 
but it has a sparse representation in some basis :math:`B`. We have 

.. math::

    \alpha  = B^H x 

as the representation of :math:`x` in the basis :math:`B`.

The regularization is then applied to the representation :math:`\alpha`.
:eq:`ista_l1_l2_minimization` becomes:

.. math::

    \widehat{x} = \text{arg} \min_{x} \frac{1}{2}\| b - A x \|_2^2 + \lambda \| B^H x \|_1 

We can rewrite this as:

.. math::

    \widehat{\alpha} = \text{arg} \min_{\alpha} \frac{1}{2}\| b - A B \alpha \|_2^2 + \lambda \| \alpha \|_1 

:eq:`ista_iteration_thresholding` changes to:

.. math::

    \alpha_{i+1} = \mathbf{T}_{\frac{\lambda}{c} }\left (\alpha_i + \frac{1}{c} B^H A^H (b - A B \alpha_i)\right )

By substituting :math:`\alpha = B^H x` and :math:`x = B \alpha`, we get:

.. math::

    \alpha_{i+1} = \mathbf{T}_{\frac{\lambda}{c} }\left (B^H x_i + \frac{1}{c} B^H A^H (b - A x_i)\right )

Simplifying further:

.. math::
    :label: ista_iteration_thresholding_basis

    x_{i+1} = B \mathbf{T}_{\frac{\lambda}{c} }\left ( B^H  \left (x_i + \frac{1}{c} A^H (b - A x_i)\right ) \right )

This is the version of IST algorithm with an operator :math:`A` and a basis :math:`B`.

Implementation
----------------------

We introduce the current residual:

.. math::

    r = b - A x

and a step size parameter :math:`\alpha = \frac{1}{c}`.
We also assume that a thresholding function :math:`\mathbf{T}` will be user defined.

This simplifies the iteration to:

.. math::

    x_{i+1} = B \mathbf{T}\left ( B^H \left (x_i + \alpha A^H r_i \right ) \right )


.. rubric:: Algorithm state

Current state of the algorithm is described by following quantities:

.. list-table::
    :header-rows: 1

    * - Term
      - Description
    * - ``x``
      - Current estimate of :math:`x`
    * - ``r``
      - Current residual :math:`r`
    * - ``r_norm_sqr``
      - Squared norm of the residual :math:`\| r \|_2^2`
    * - ``x_change_norm``
      - Change in the norm of :math:`x` given by :math:`\|x_{i+1} - x_{i} \|_2`
    * - ``iterations``
      - Number of iterations of algorithm run so far

.. rubric:: Algorithm initialization

We initialize the algorithm with:

* :math:`x \leftarrow x_0` an initial estimate of solution given by user.
  It can be 0.
* :math:`r \leftarrow b - A x_0`
* Compute ``r_norm_sqr`` 
* Give a very high value to ``x_change_norm`` (since there is no :math:`x_{-1}`).
* ``iterations = 0``


.. rubric:: Algorithm iteration

Following steps are involved in each iteration. These are
directly based on :eq:`ista_iteration_thresholding_basis`.

* Compute gradient:  :math:`g \leftarrow \alpha A^H r`
* Update estimate: :math:`x \leftarrow x + g`
* Transform :math:`\alpha \leftarrow B^H x`
* Threshold: :math:`\alpha \leftarrow \mathbf{T} (\alpha)`
* Inverse transform :math:`x \leftarrow B \alpha`
* Update residual: :math:`r \leftarrow b - A x`
Underdetermined Linear Systems
=================================


The discussion in this section is largely based on chapter 1 of 
:cite:`elad2010sparse`.

Consider a matrix :math:`\Phi \in \CC^{M \times N}` with :math:`M < N`. 

Define an under-determined system of linear equations:

.. math::

  \Phi x = y

where :math:`y \in \CC^M` is known and :math:`x \in \CC^N` is unknown. 

This system has :math:`N` unknowns and
:math:`M` linear equations. 
There are more unknowns than equations.

Let the columns of :math:`\Phi` be given by :math:`\phi_1, \phi_2, \dots, \phi_N`.

Column space of :math:`\Phi` (vector space spanned by all columns of :math:`\Phi`)  is denoted by :math:`\ColSpace(\Phi)`
i.e.

.. math::

  \ColSpace(\Phi) = \sum_{i=1}^{N} c_i \phi_i, \quad c_i \in \CC.


We know that :math:`\ColSpace(\Phi) \subset \CC^M`. 

Clearly :math:`\Phi x \in \ColSpace(\Phi)` for every :math:`x \in \CC^N`.  Thus if :math:`y \notin \ColSpace(\Phi)` then we have no solution. But, if :math:`y \in \ColSpace(\Phi)` then we have infinite number of solutions.

Let :math:`\NullSpace(\Phi)` represent the null space of :math:`\Phi` given by 

.. math::

  \NullSpace(\Phi) = \{ x \in \CC^N : \Phi x = 0\}.


Let :math:`\widehat{x}` be a solution of :math:`y = \Phi x`. And let :math:`z \in \NullSpace(\Phi)`. Then 

.. math::

  \Phi (\widehat{x} + z) = \Phi \widehat{x} + \Phi z = y + 0  = y.

Thus the set :math:`\widehat{x} + \NullSpace(\Phi)` forms the complete set of infinite solutions to the
problem :math:`y = \Phi x` where

.. math::

  \widehat{x} + \NullSpace(\Phi) = \{\widehat{x} + z \quad \Forall z \in \NullSpace(\Phi)\}.

Regularization
------------------------
One way to pick a specific solution from the set of all possible solutions is to introduce **regularization**. 

.. index:: Regularization

We define a cost function :math:`J(x) : \CC^N \to \RR` which defines the **desirability** of a given solution :math:`x` out
of infinitely possible solutions. The higher the cost, lower is the desirability of
the solution.

.. index:: Desirability

Thus the goal of the optimization problem is to find a desired :math:`x` with minimum possible cost.

We can write this optimization problem as
  
.. math::

  \begin{aligned}
    & \underset{x}{\text{minimize}} 
    & &  J(x) \\
    & \text{subject to}
    & &  y = \Phi x.
  \end{aligned}



If :math:`J(x)` is convex, then its possible to find a global minimum cost solution over the solution set.

If :math:`J(x)` is not convex, then it may not be possible to find a global minimum, we may have to
settle with a local minimum. 

A variety of such cost function based criteria can be considered. 

:math:`l_2` Regularization
--------------------------------

One of the most common criteria is to choose a solution with the smallest :math:`l_2` norm.

The problem can then be reformulated as an optimization problem 
  
.. math::

  \begin{aligned}
    & \underset{x}{\text{minimize}} 
    & &  \| x \|_2 \\
    & \text{subject to}
    & &  y = \Phi x.
  \end{aligned}


In fact minimizing :math:`\| x \|_2` is same as minimizing its square :math:`\| x \|_2^2 = x^H x`.

So an equivalent formulation is 

  
.. math::

  \begin{aligned}
    & \underset{x}{\text{minimize}} 
    & &  x^H x \\
    & \text{subject to}
    & &  y = \Phi x.
  \end{aligned}



A formal solution to :math:`l_2` norm minimization problem can be easily obtained using
Lagrange multipliers.

We define the Lagrangian
  
.. math::

  \mathcal{L}(x) = \|x\|_2^2 + \lambda^H (\Phi x  - y)

with :math:`\lambda \in \CC^M` being the Lagrange multipliers for the (equality) constraint set.

Differentiating :math:`\mathcal{L}(x)` w.r.t. :math:`x` we get
  
.. math::

  \frac{\partial \mathcal{L}(x)} {\partial x} = 2 x + \Phi^H \lambda.


By equating the derivative to :math:`0` we obtain the optimal value of :math:`x` as
  
.. math::

  x^* = - \frac{1}{2} \Phi^H \lambda.
  \label{eq:ssm:underdetermined_l2_optimal_value_expression_1}


Plugging this solution back into the constraint :math:`\Phi x= y` gives us
  
.. math::

  \Phi x^* = - \frac{1}{2} (\Phi \Phi^H) \lambda= y\implies \lambda = -2(\Phi \Phi^H)^{-1} y.


In above we are implicitly assuming that :math:`\Phi` is a full rank matrix thus, :math:`\Phi \Phi^H` is invertible
and positive definite.

Putting :math:`\lambda` back in above we obtain
the well known closed form least squares solution using pseudo-inverse solution
  
.. math::

  x^* = \Phi^H (\Phi \Phi^H)^{-1} y = \Phi^{\dag} y.

Convexity
------------------
Convex optimization
problems have a unique feature that it is possible to find the global optimal solution if
such a solution exists. 

The solution space  :math:`\Omega = \{x : \Phi x = y\}` is convex.
Thus the feasible set of solutions for the optimization problem
is also convex. All it remains is to make sure that we choose a cost function
:math:`J(x)` which happens to be convex. This will ensure that a global minimum can be found through
convex optimization techniques. Moreover, if :math:`J(x)` is strictly convex, then it is guaranteed
that the global minimum solution is *unique*. Thus even though, we may not have
a nice looking closed form expression for the solution of a strictly convex cost function minimization problem,
the guarantee of the existence and uniqueness of solution as well as well developed algorithms
for solving the problem make it very appealing to choose cost functions which are convex.

We recall that all :math:`l_p` norms with :math:`p \geq 1` are convex functions.
In particular :math:`l_{\infty}` and :math:`l_1` norms are very interesting and popular where
  
.. math::

  l_{\infty}(x) = \max(|x_i|), \, 1 \leq i \leq N

and
  
.. math::

  l_1(x) = \sum_{i=1}^{N} |x_i|.


In the following section we will attempt to find a unique solution to our 
optimization problem using :math:`l_1` norm.

:math:`l_1` Regularization
-----------------------------------

In this section we will restrict our attention to the
Euclidean space case where :math:`x \in \RR^N`,
:math:`\Phi \in \RR^{M \times N}` and :math:`y \in \RR^M`.

We choose our cost function :math:`J(x) = l_1(x)`.

The cost minimization problem can be reformulated as
  
.. math::

  \begin{aligned}
    & \underset{x}{\text{minimize}} 
    & &  \| x \|_1 \\
    & \text{subject to}
    & &  \Phi x = y.
  \end{aligned}



It's time to have a closer look at our cost function :math:`J(x) = \|x \|_1`. This function
is convex yet not strictly convex. 

For the :math:`l_1` norm minimization problem since :math:`J(x)` is not strictly convex,
hence a unique solution may not be guaranteed. In specific cases, there may be
infinitely many solutions. Yet what we can claim is

* these solutions are gathered in a set that is bounded and convex, and
* among these solutions, there exists at least one solution with at most
  :math:`M` non-zeros (as the number of constraints in :math:`\Phi x = y`).


.. theorem::

  Let :math:`S` denote the solution set of :math:`l_1` norm minimization problem.
  :math:`S`
  contains at least one solution :math:`\widehat{x}` with
  :math:`\| \widehat{x} \|_0 = M`.

See :cite:`elad2010sparse` for proof.


We thus note that :math:`l_1` norm has a tendency to prefer sparse solutions. This is a
well known and fundamental property of linear programming.

:math:`l_1` norm minimization problem as a linear programming problem
------------------------------------------------------------------------

We now show that :math:`l_1` norm minimization problem in :math:`\RR^N` 
is in fact a linear programming problem.

Recalling the problem:
    
.. math::
  :label: l1_norm_minimization_problem

  \begin{aligned}
    & \underset{x \in \RR^N}{\text{minimize}} 
    & &  \| x \|_1 \\
    & \text{subject to}
    & &  y = \Phi x.
  \end{aligned}


Let us write :math:`x` as :math:`u  - v`  where :math:`u, v \in \RR^N` are both non-negative vectors such that
:math:`u` takes all positive entries in :math:`x` while :math:`v` takes all the negative entries in :math:`x`.

We note here that by definition
    
.. math::

  \supp(u) \cap \supp(v) = \EmptySet

i.e. support of :math:`u` and :math:`v` do not overlap.

We now construct a vector
    
.. math::

  z = \begin{bmatrix}
  u \\ v
  \end{bmatrix} \in \RR^{2N}.


We can now verify that
    
.. math::

  \| x \|_1 = \|u\|_1 + \| v \|_1 = 1^T z.


And 
    
.. math::

  \Phi x = \Phi (u - v) = \Phi u - \Phi v = 
  \begin{bmatrix}
  \Phi & -\Phi
  \end{bmatrix}
  \begin{bmatrix}
  u \\ v
  \end{bmatrix}
  = \begin{bmatrix}
  \Phi & -\Phi
  \end{bmatrix} z 

where  :math:`z \succeq 0`.

Hence the optimization problem :eq:`l1_norm_minimization_problem` can be recast as
    
.. math::
  :label: l1_norm_minimization_problem_as_lp

  \begin{aligned}
    & \underset{z \in \RR^{2N}}{\text{minimize}} 
    & &  1^T z \\
    & \text{subject to}
    & &  \begin{bmatrix} \Phi & -\Phi \end{bmatrix} z = y\\
    & \text{and}
    & & z \succeq 0.
  \end{aligned}


This optimization problem has the classic Linear Programming structure since the
objective function is affine as well as constraints are affine.

Let :math:`z^* =\begin{bmatrix} u^* \\ v^* \end{bmatrix}` be an optimal solution to the
problem :eq:`l1_norm_minimization_problem_as_lp`.  

In order to show that the two optimization problems are equivalent, we need
to verify that our assumption about the decomposition of :math:`x` into positive entries in :math:`u` 
and negative entries in :math:`v` is indeed satisfied by the optimal solution :math:`u^*` and :math:`v^*`.
i.e. support of :math:`u^*` and :math:`v^*` do not overlap.

Since :math:`z \succeq 0` hence :math:`\langle u^* , v^* \rangle  \geq 0`. If support of :math:`u^*` and :math:`v^*` 
don't overlap, then we  have :math:`\langle u^* , v^* \rangle = 0`. And if they overlap then
:math:`\langle u^* , v^* \rangle > 0`.

Now for the sake of contradiction, let us assume that support of :math:`u^*` and :math:`v^*` do overlap for 
the optimal solution :math:`z^*`.

Let :math:`k` be one of the indices at which both :math:`u_k \neq 0` and :math:`v_k \neq 0`. Since :math:`z \succeq 0`, 
hence :math:`u_k > 0` and :math:`v_k > 0`.

Without loss of generality let us assume that :math:`u_k > v_k > 0`.

In the equality constraint
    
.. math::

  \begin{bmatrix} \Phi & -\Phi \end{bmatrix} \begin{bmatrix} u \\ v \end{bmatrix} = y


Both of these coefficients multiply the same column of :math:`\Phi` with opposite signs giving us a term
    
.. math::

  \phi_k (u_k - v_k). 


Now if we replace the two entries in :math:`z^*` by
    
.. math::

  u_k'  = u_k - v_k

and
    
.. math::

  v_k' = 0

to obtain an new vector :math:`z'`, 
we see that there is no impact in the equality constraint 
    
.. math::

  \begin{bmatrix} \Phi & -\Phi \end{bmatrix} z = y.

Also the positivity constraint
    
.. math::

  z \succeq 0

is satisfied. This means that :math:`z'` is a feasible solution.

On the other hand, the objective function :math:`1^T z` value reduces by :math:`2 v_k` for :math:`z'`. 
This contradicts our assumption that :math:`z^*` is the optimal solution.

Hence for the optimal solution of :eq:`l1_norm_minimization_problem_as_lp`
we have
    
.. math::

  \supp(u^*) \cap \supp(v^*) = \EmptySet

thus 
    
.. math::

  x^* = u^* - v^*

is indeed the desired solution for the optimization problem :eq:`l1_norm_minimization_problem`.


Sparse Linear Systems
==========================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   underdetermined
   thresholding
   ista
.. _sls:thresholding:

Thresholding
==========================

.. contents::
    :depth: 2
    :local:


For a real or complex scalar :math:`x`, the **hard thresholding** operator can be defined as:

.. math::

    \mathbf{HT}_{\gamma}(x) = \begin{cases} 
     x & \text{for} & |x| \gt \gamma\\
     0 & \text{for} & |x| \le \gamma
    \end{cases}

For a multi-dimensional array, we apply the same operator on each entry in the array.

The **soft thresholding** operator for real and complex scalars can be defined as:

.. math::

    \mathbf{ST}_{\gamma}(x) = \begin{cases} 
     x\left ( 1 -  \frac{\gamma}{|x|} \right ) & \text{for} & |x| \gt \gamma\\
     0 & \text{for} & |x| \le \gamma
    \end{cases}

For real numbers, this reduces to:

.. math::

    \mathbf{ST}_{\gamma}(x) = \begin{cases} 
     x - \gamma & \text{for} & x \gt \gamma \\
     x + \gamma & \text{for} & x \lt -\gamma \\
     0 & \text{for} & |x| \le \gamma
    \end{cases}


:cite:`chen2014irregular` consider the general minimization problem:

.. math::

    \widehat{x} = \text{arg} \min_{x} \| b - A x \|_2^2 + \mathbf{R}(x)

where :math:`x` is a vector in the model space, :math:`b` is a vector in the data space,
:math:`A` is a linear operator from model space to data space, and
:math:`\mathbf{R}(x)` is the regularization term on the model vector.
They specifically consider the regularizations

.. math:: 

    \mathbf{R}(x) =  \tau \|x\|_p^p
    
for :math:`p`-norms where :math:`0 \leq p \leq 1`. 
This can be solved by the IST (Iterative Shrinkage and Thresholding)
algorithm where the thresholding operator depends on the selection of :math:`p`-norm and the
regularization parameter :math:`\tau`.
They describe the IST iterations in terms of a more general thresholding operator :math:`\mathbf{T}_{\gamma(\tau, p)}(x)`:

.. math::

    x_{n+1} = \mathbf{T}_{\gamma(\tau, p)}\left [x_n + A^H (b - A x_n) \right ]


They provide the definition of the thresholding operator for:

* :math:`p=0` **hard thresholding**
* :math:`p=1` **soft thresholding**
* :math:`p=1/2` **half thresholding**


.. rubric:: Hard thresholding

Whe :math:`p=0`, we have:

.. math::

    \mathbf{R}(x) = \| x \|_0

.. math::

    \gamma(\tau, 0) = \sqrt{2 \tau}

The hard thresholding operator reduces to:

.. math::

    \mathbf{T}_{\gamma(\tau, 0)}(x) = \begin{cases} 
     x & \text{for} & |x| \gt \gamma (\tau, 0)\\
     0 & \text{for} & |x| \le \gamma  (\tau, 0)
    \end{cases}

.. rubric:: Soft thresholding

Whe :math:`p=1`, we have:

.. math::

    \mathbf{R}(x) = \| x \|_1

.. math::

    \gamma(\tau, 1) = \tau

The soft thresholding operator reduces to:

.. math::

    \mathbf{T}_{\gamma(\tau, 1)}(x) = \begin{cases} 
     x\left ( 1 -  \frac{\gamma}{|x|} \right ) & \text{for} & |x| \gt \gamma (\tau, 1)\\
     0 & \text{for} & |x| \le \gamma (\tau, 1)
    \end{cases}


.. rubric:: Half thresholding

Whe :math:`p=\frac{1}{2}`, we have:

.. math::

    \mathbf{R}(x) = \| x \|_{\frac{1}{2}}^{\frac{1}{2}}

.. math::

    \gamma(\tau, \frac{1}{2}) = \frac{3}{2} \tau^{\frac{2}{3}}


The half thresholding operator is more complicated:


.. math::

    \mathbf{T}_{\gamma(\tau, 1)}(x) = \begin{cases} 
     \frac{2}{3} x\left ( 1 +  \cos \left ( \frac{2}{3} \pi - \frac{2}{3} \arccos \left ( \frac{\tau}{8} \left (\frac{|x|}{3} \right)^{\frac{3}{2}}
         \right )   \right )  \right ) 
     & \text{for} & |x| \gt \gamma (\tau, \frac{1}{2})\\
     0 & \text{for} & |x| \le \gamma (\tau, \frac{1}{2})
    \end{cases}
Computation Time Comparison of Sparse Recovery Methods
===============================================================

Performance on CPU
-------------------------------------------------

.. rubric:: System configuration

* MacBook Pro 2019 Model
* Processor: 1.4 GHz Quad Core Intel Core i5
* Memory: 8 GB 2133 MHz LPDDR3

.. rubric:: Problem Specification

* Gaussian sensing matrices (normalized to unit norm columns)
* Sparse vectors with non-zero entries drawn from Gaussian distributions
* M, N, K have been chosen so that all algorithms under comparison are known to converge to successful 
  recovery.

.. rubric:: Remarks

* All algorithms have been benchmarked for both 32-bit and 64-bit floating point calculations. Benchmarks are separately presented for them.
* It was separately verified that sparse recovery results were identical for both with or without JIT acceleration.
* Python ``%timeit`` magic was used for benchmarking. 
* Every algorithm has been run several times on the given problem and the average time has been computed.
* Average times have been reported in jit_off and jit_on columns with milliseconds units.

.. rubric:: Algorithm structures

The table below highlights the differences
in the structure of different algorithms under consideration. These differencess are key reason for the computational
complexity.

.. list-table:: Comparison of algorithm structures
    :header-rows: 1

    * - method
      - Correlation with residual
      - Least squares
      - Hard thresholding
      - Step size
    * - OMP
      - Yes
      - Cholesky update
      - 1 atom
      - No
    * - SP
      - Yes
      - 2 (2K atoms and K atoms)
      - K atoms and K atoms
      - No
    * - CoSaMP
      - Yes
      - 1 (3K atoms)
      - 2K atoms and K atoms
      - No
    * - IHT
      - Yes
      - 0
      - K atoms
      - Fixed
    * - NIHT
      - Yes
      - 0
      - K atoms
      - Dynamic
    * - HTP
      - Yes
      - 1 (K atoms)
      - K atoms
      - Fixed
    * - NHTP
      - Yes
      - 1 (K atoms)
      - K atoms
      - Dynamic




.. rubric:: Benchmarks for 32-bit

.. list-table:: Average time (msec) and speedups due to JIT acceleration
    :header-rows: 1

    * - method
      - M
      - N
      - K
      - iterations
      - jit_off
      - jit_on
      - speedup
    * - OMP
      - 200
      - 1000
      - 20
      - 20
      - 105.78
      - 2.14
      - 49.48
    * - SP
      - 200
      - 1000
      - 20
      - 3
      - 1645.32
      - 2.73
      - 602.34
    * - CoSaMP
      - 200
      - 1000
      - 20
      - 4
      - 309.01
      - 6.20
      - 49.84
    * - IHT
      - 200
      - 1000
      - 20
      - 65
      - 232.99
      - 36.27
      - 6.42
    * - NIHT
      - 200
      - 1000
      - 20
      - 16
      - 240.96
      - 5.64
      - 42.72
    * - HTP
      - 200
      - 1000
      - 20
      - 5
      - 1491.00
      - 13.71
      - 108.76
    * - NHTP
      - 200
      - 1000
      - 20
      - 4
      - 1467.35
      - 1.98
      - 741.88


.. rubric:: Benchmarks for 64-bit

.. list-table:: Average time (msec) and speedups due to JIT acceleration
    :header-rows: 1

    * - method
      - M
      - N
      - K
      - iterations
      - jit_off
      - jit_on
      - speedup
    * - OMP
      - 200
      - 1000
      - 20
      - 20
      - 112.69
      - 2.43
      - 46.42
    * - SP
      - 200
      - 1000
      - 20
      - 4
      - 1324.79
      - 4.49
      - 295.02
    * - CoSaMP
      - 200
      - 1000
      - 20
      - 5
      - 293.50
      - 9.82
      - 29.90
    * - IHT
      - 200
      - 1000
      - 20
      - 77
      - 209.22
      - 48.81
      - 4.29
    * - NIHT
      - 200
      - 1000
      - 20
      - 19
      - 196.66
      - 7.23
      - 27.21
    * - HTP
      - 200
      - 1000
      - 20
      - 6
      - 1218.62
      - 18.96
      - 64.28
    * - NHTP
      - 200
      - 1000
      - 20
      - 5
      - 1238.37
      - 2.79
      - 443.68
      
​

Subjective Analysis 
----------------------------

.. rubric:: 64-bit vs 32-bit

* There are differences in number of iterations for convergence
* Every algorithm except OMP takes more iterations to converge with 64-bit compared to 32-bit floating point computations.
* In case of OMP, number of iterations is decided by sparsity. Hence, it is same for both 32-bit and 64-bit.
* It was separately established that success rates of these algorithms suffers somewhat for 32-bit floating point calculations.
* In other words, 32-bit computations are more aggressive and may be inaccurate.
* On CPUs, the floating point units are 64-bit. Hence, using 32-bit floating point computations doesn't give us much speedup. 
  32-bit computation would be more relevant for GPUs.
* The general trend of computation times (with JIT on) for both 32-bit and 64-bit are similar. i.e. algorithms which are
  slower for 32-bit are slower for 64-bit too.

Rest of the discussion is focused on the results for 64-bit sparse recovery. 
.. rubric:: All algorithms without JIT vs with JIT

* It is clear that all algorithms exhibit significant speedups with the introduction of 
  JIT acceleration.
* The speedup is as low as 4x for IHT and as high as 443x in NHTP.
* Before JIT, OMP is the fastest algorithm and SP is the slowest. 
* After JIT acceleration, OMP is the fastest algorithm while IHT is the slowest. NHTP comes as a close second. 
  Incidentally, NHTP is faster than OMP for 32-bit.
* NHTP and SP show significant speedups with JIT. HTP, OMP, CoSaMP and NIHT show modest gains. IHT doesn't seem to provide much 
  optimization opportunities.
* It appears that steps like dynamic step size computation (in NIHT, NHTP) and 
  least squares (in SP, CoSaMP, HTP, NHTP)
  tend to get aggressively optimized and lead to massive speed gains.

.. rubric:: OMP

* With JIT on, OMP is actually one of the fastest algorithms in the mix (for both 32-bit and 64-bit).
* In the current implementations, OMP is the only one in which the least squares step has
  been optimized using Cholesky updates. 
* This is possible as OMP structure allows for adding atoms one at a time to the mix.
* Other algorithms change several atoms [add / remove] in each iteration. Hence, such
  optimizations are not possible.
* The least squares steps in other algorithms can be accelerated using small number of conjugate gradients
  iterations. However, this hasn't been implemented yet.


.. rubric:: SP vs CoSaMP

* CoSaMP has one least squares step (on 3K indices) in each iteration.
* SP (Subspace Pursuit) has two least squares steps in each iteration.
* Without JIT, CoSaMP is 4x faster.
* With JIT, SP becomes 2x faster than CoSaMP.
* Thus, SP seems to provide more aggressive optimization opportunities.

.. rubric:: IHT vs NIHT

* IHT and NIHT are both simple algorithms. They don't involve a least squares step in their iterations.
* The main difference is that the step-size fixed for IHT and it is computed on every iteration in NIHT.
* The dynamic step size leads to reduction in the number of iterations for NIHT. From 77 to 19, 4x reduction.
* Without JIT, there is no significant difference between IHT and NIHT.
  Thus, step-size computation seems to contribute a lot to computation time without acceleration.
* With JIT, step-size computation seems to be aggressively optimized.
  NIHT after JIT is 6x faster than IHT even though the number of iterations reduces by only 4 times
  and there is extra overhead of computing the step size. This appears to be counter-intuitive.

.. rubric:: IHT vs HTP

* The major difference in the two algorithms is that HTP performs a least squares estimate
  on the current guess of signal support
* The number of iterations reduces 13 times due to the least squares step but it has its own extra overhead.
* Without JIT, HTP becomes much slower than IHT (6x slower). Thus, overhead of a least squares step is quite high.
* HTP is about 3x faster than IHT with JIT. This makes sense. The number of iterations reduced by 13
  times and the overhead of least squares was added.

.. rubric:: HTP vs NHTP

* Just like NIHT, NHTP also introduces computing the step size dynamically in every iteration.
* It helps in reducing the number of iterations from 6 to 5.
* In this case, the benefit of dynamic step size is not visible much in terms of iterations.
* Without JIT, NHTP is somewhat slower than HTP.
* However, with JIT, NHTP is 6x faster than HTP. This speedup is unusual as there is just
  20% reduction in number of iterations and there is the overhead of step size computation.

 

Benchmarks
=====================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   comparison
   ompOrthogonal Matching Pursuit
=============================



.. rubric:: Speed benchmarks for JAX implementation

Each row of the following table describes:

* problem type and configuration (M x N is dictionary size, K is sparsity level)
* Average time taken in CPU/GPU configurations
* Speed improvement ratios

.. rubric:: System used

* All benchmarks have been generated on Google Colab
* CPU and GPU configurations Google Colab have been used

.. list-table::
    :header-rows: 1

    * - M
      - N
      - K 
      - CPU
      - CPU + JIT
      - CPU / CPU + JIT
      - GPU 
      - GPU + JIT
      - GPU / GPU + JIT
      - CPU + JIT / GPU + JIT
    * - 256
      - 1024
      - 16
      - 148 ms
      - 8.27 ms
      - 17.9x
      - 139 ms
      - 1.28 ms
      - 108x
      - 6.46x

.. rubric:: Observations

* JIT (Just In Time) compilation seems to give significant performance improvements 
  in both CPU and GPU architectures
* Current implementation seems to be slower on GPU vs CPU with JIT. 
* GPU speed gain over CPU (with JIT on) is relatively meager. 
  On TensorFlow, people regularly report 30x improvements between CPU to GPU 
  for neural networks implemented using Keras. 


.. rubric:: Possible deficiencies

* There is opportunity to improve parallelization in the OMP implementation.
* Cholesky update based implements depends heavily on solving triangular systems.
* GPUs may not be great at solving triangular systems. 


 Sparse Subspace Clustering 
 ==================================

 .. toctree::
   :maxdepth: 2
   :caption: Contents:

   intro
Introduction to Sparse Subspace Clustering 
==============================================

Consider a dataset of :math:`S` points in the ambient data space
:math:`\mathbb{R}^M` which have been assembled in a matrix :math:`Y` of shape :math:`M \times S`.

In many applications, it often occurs that
if we *group* or *segment* the data set :math:`Y` into
multiple disjoint subsets (clusters): 
:math:`Y = Y_1 \cup \dots \cup Y_K`,
then each subset can be modeled sufficiently well by a low dimensional subspace
:math:`\mathbb{R}^D` where :math:`D \ll M`.
Some of the applications include:
motion segmentation :cite:`tomasi1991detection,tomasi1992shape, 
boult1991factorization,
poelman1997paraperspective,
gear1998multibody,
costeira1998multibody,
kanatani2001motion`, 
face clustering :cite:`basri2003lambertian, ho2003clustering, lee2005acquiring`
and handwritten digit recognition :cite:`zhang2012hybrid`.

*Subspace clustering* is a clustering framework which assumes
that the data-set can be segmented into clusters where points in
different clusters are drawn from different subspaces. Subspace clustering
algorithms are able to simultaneously segment the data into 
clusters corresponding to different subspaces as well as estimate
the subspaces from the data itself.
A comprehensive review of subspace clustering can be found in 
:cite:`vidal2010tutorial`.
Several state of the art algorithms are based on building
subspace preserving representations of individual data points
by treating the data set itself as a (self expressive) dictionary.
For creating subspace preserving representations, one resorts to
using sparse coding algorithms developed in sparse representations and 
compressive sensing literature. 

Two common algorithms are
*Sparse Subspace Clustering* using :math:`\ell_1` *regularization*
(SSC-:math:`\ell_1`):cite:`elhamifar2009sparse, elhamifar2013sparse` 
and *Sparse Subspace Clustering using Orthogonal
Matching Pursuit* (SSC-OMP) :cite:`dyer2013greedy, you2015sparse, you2016scalable`. 
While SSC-:math:`\ell_1` is guaranteed to give correct clustering under
broad conditions (arbitrary subspaces and corrupted data), it
requires solving a large scale convex optimization problem. On
the other hand, SSC-OMP 
is computationally efficient but its clustering accuracy is
poor (especially at low density of data points per subspace).

The dataset :math:`Y` is modeled as being sampled from a collection
or arrangement :math:`\mathcal{U}` of linear (or affine) subspaces
:math:`\mathcal{U}_k \subset \mathbb{R}^M` : 
:math:`\mathcal{U} = \{ \mathcal{U}_1  , \dots , \mathcal{U}_K \}`. 
The union of the subspaces
is denoted as
:math:`Z_{\mathcal{U}} = \mathcal{U}_1 \cup \dots \cup \mathcal{U}_K`.

Let the data set be :math:`\{ y_j  \in \mathbb{R}^M \}_{j=1}^S`
drawn from the union of subspaces under consideration.
:math:`S` is the total number of data points being analyzed
simultaneously.
We put the data points together in a *data matrix* as

.. math::

    Y  \triangleq \begin{bmatrix}
    y_1 & \dots & y_S
    \end{bmatrix}.

Let the vectors be drawn from a set of :math:`K` (linear or affine) subspaces, 
The subspaces are indexed by a variable :math:`k` with :math:`1 \leq k \leq K`.
The :math:`k`-th subspace is denoted by :math:`\mathcal{U}_k`. 
Let the (linear or affine) dimension
of :math:`k`-th subspace be :math:`\dim(\mathcal{U}_k) = D_k` with :math:`D_k \leq D \ll M`.

The vectors in :math:`Y` can be grouped (or segmented or clustered) 
as submatrices 
:math:`Y_1, Y_2, \dots, Y_K` such 
that all vectors in :math:`Y_k` are drawn from the subspace :math:`\mathcal{U}_k`. 
Thus, we can write

.. math::

    Y^* = Y \Gamma = \begin{bmatrix} y_1 & \dots & y_S \end{bmatrix} 
    \Gamma
    = \begin{bmatrix} Y_1 & \dots & Y_K \end{bmatrix} 

where :math:`\Gamma` is an :math:`S \times S` unknown permutation
matrix placing each vector to the right subspace. 

Let there be :math:`S_k` vectors in :math:`Y_k` with
:math:`S = S_1 + \dots + S_K`. 
Let :math:`Q_k` be an orthonormal basis for subspace :math:`\mathcal{U}_k`. Then,
the subspaces can be described as 

.. math::

    \mathcal{U}_k = \{ y \in \mathbb{R}^M : y = \mu_k + Q_k \alpha \}, \quad 1 \leq k \leq K 

For linear subspaces, :math:`\mu_k = 0`.

A dataset where each point can be expressed as a linear combination
of other points in the dataset is said to satisfy 
*self-expressiveness property*. The self-expressive 
representation of a point :math:`y_s` in :math:`Y` is given by 

.. math::

    y_s = Y c_s, \; c_{ss} = 0, \text{ or } Y = Y C, \quad \text{diag}(C) = 0

where :math:`C = \begin{bmatrix}c_1, \dots, c_S \end{bmatrix} \in \mathbb{R}^{S \times S}` 
is the matrix of representation coefficients. 

Let :math:`y_s` belong to :math:`k`-th subspace :math:`\mathcal{U}_k`. 
Let :math:`Y^{-s}` denote the dataset :math:`Y` excluding the point :math:`y_s` 
and  :math:`Y_k^{-s}` denote the
set of points in :math:`Y_k` excluding :math:`y_s`. If :math:`Y_k^{-s}` spans the subspace
:math:`\mathcal{U}_k`, then a representation of :math:`y_s` can be constructed entirely
from the points in :math:`Y_k^{-s}`. A representation is called 
*subspace preserving* if it consists of points within the same subspace.

If :math:`c_i` is a subspace preserving representation of :math:`y_i` and :math:`y_j`
belongs to a different subspace, then :math:`c_{ij} = 0`. Thus, if :math:`C` consists
entirely of subspace preserving representations, then :math:`C_{ij} = 0` whenever
:math:`y_i` and :math:`y_j` belong to different subspaces. 
In other words, if :math:`Y_{-k}` denotes the set of points from 
all subspaces excluding the subspace :math:`Y_k` corresponding
to the point :math:`y_i`, then points in :math:`Y_{-k}` do not
participate in the representation :math:`c_i`.

In the ``cr.sparse.cluster.ssc`` package, we provide a version of
OMP which can be used to construct the sparse self expressive representations 
:math:`C` of :math:`Y`. Once the representation has been constructed, we compute an
affinity matrix :math:`W = |C| + |C^T|`. 

We then apply spectral clustering on :math:`W` to complete SSC-OMP. 
For this, we have written a JAX version of spectral clustering
in ``cr.sparse.cluster.spectral`` package. In particular, it
uses our own Lanczos Bidiagonalization with Partial Orthogonalization (LANBPRO)
algorithm to compute the :math:`K` largest singular values of the
normalized affinity matrix
in as few iterations as possible. The intermediate variables 
:math:`C`, :math:`W` are maintained in the experimental sparse matrices
stored in BCOO format.
The LANBPRO algorithm also works on sparse matrices directly. 
Thus, even though :math:`C` is of size :math:`S \times S`, it can be stored 
efficiently in :math:`O(DS)` storage. This enables us to process 
hundreds of thousands of points efficiently. 
Limitations
======================


This section is a rough summary of some of the limitations 
encountered during the development of this package. 
Some of these limitations are due to not enough time 
spent on development. Others are due to lack of support
in JAX library.
Be aware that some of these limitations may be entirely 
due to my limited knowledge of JAX library. Hopefully, 
in future, with my better knowledge, or with better support
from JAX, some of these limitations will be alleviated.


.. rubric:: Key issues with JAX 

This is a list of key issues with JAX library which have
been acknowledged with JAX development team also.

* Lack of support for dynamic or data dependent shapes in JAX.
  See `#8042 <https://github.com/google/jax/discussions/8042>`_
* 1D convolution is slow in JAX (especially on CPU). 
  See `#7961 <https://github.com/google/jax/discussions/7961>`_.
* Support for sparse matrices is still under development in JAX.


Utilities
----------------------

- off_diagonal_elements cannot be jitted. 
  JAX arrays do not support boolean scalar indices.

Data clustering
-----------------------

- Spectral clustering assess the number of clusters 
  from data. It cannot be jitted. In turn, the k-means
  invocation from inside spectral clustering cannot be jitted either.
- Normalized spectral clustering could make use of 
  sparse matrices (identity and diagonal) whenever
  it becomes available.

Signal Processing
----------------------

- Walsh Hadamard transform implementation is technically correct 
  but possibly not vectorized enough. It is slow.

Linear Operators
---------------------

- Some operators like pad don't work well on matrix input yet.Source Code
===================



Working with the source code in development mode
-----------------------------------------------------


Clone the repository::

    git clone https://github.com/carnotresearch/cr-sparse.git


Change into the code::

    cd cr-sparse


Install the package in development mode::

    python -m pip install -e .
Vector Norms
================================

Our interest is in operators mapping vectors from a model space :math:`\mathbb{C}^n` to a data space :math:`\mathbb{C}^m`.

There are some simple and useful results on relationships between 
different  :math:`p`-norms listed in this section. We also discuss
some interesting properties of  :math:`l_1`-norm specifically.

.. _def:ssm:sign_vector:

.. definition::

     
    .. index:: Sign vector
    
    
    Let  :math:`v \in \mathbb{C}^n`. Let the entries in  :math:`v` be represented as
    
    .. math::
    
        v_i = r_i \exp (j \theta_i)
    
    where  :math:`r_i = | v_i |` with the convention that  :math:`\theta_i = 0` 
    whenever  :math:`r_i = 0`.
    
    The sign vector for  :math:`v` denoted by  :math:`\sgn(v)` is defined as
    
    .. math::
    
        \sgn(v)  = \begin{bmatrix}\sgn(v_1) \\ \vdots \\ \sgn(v_N)  \end{bmatrix}
    
    where
    
    .. math::
    
        \sgn(v_i) = \left\{
                \begin{array}{ll}
                    \exp (j \theta_i) & \mbox{if  $r_i \neq 0$ };\\
                    0 & \mbox{if  $r_i = 0$ }.
                \end{array}
              \right.
    


.. _res:ssm:l1_norm_as_inner_product_with_sign_vector:

.. lemma::


    
    For any  :math:`v \in \mathbb{C}^n` : 
    
    .. math::
    
        \| v \|_1 = \sgn(v)^H v = \langle v , \sgn(v) \rangle.
    


.. proof::

    
    .. math::
    
        \| v \|_1 = \sum_{i=1}^n r_i  = \sum_{i=1}^n \left [r_i e^{j \theta_i} \right ] e^{- j \theta_i} 
        = \sum_{i=1}^n v_i e^{- j \theta_i} = \sgn(v)^H v.
    
    Note that whenever  :math:`v_i = 0`, corresponding  :math:`0` entry in  :math:`\sgn(v)` has no effect on the sum.


.. _lem:ssm:l1_norm_l2_bounds:

.. lemma::


    
    Suppose  :math:`v \in \mathbb{C}^n`.  Then
    
    .. math::
    
         \| v \|_2 \leq \| v\|_1 \leq \sqrt{n} \| v \|_2.
    



.. proof::

    For the lower bound, we go as follows
    
    .. math::
    
        \| v \|_2^2 = \sum_{i=1}^n | v_i|^2  \leq \left ( \sum_{i=1}^n | v_i|^2  + 2 \sum_{i, j, i \neq j} | v_i | | v_j| \right )
        = \left ( \sum_{i=1}^n | v_i| \right )^2 = \| v \|_1^2.
    
    This gives us
    
    .. math::
    
        \| v \|_2 \leq \| v \|_1.
    
    
    We can write  :math:`l_1` norm as
    
    .. math::
    
        \| v \|_1 = \langle v, \sgn (v) \rangle.
    
    
    By Cauchy-Schwartz inequality we have
    
    .. math::
    
        \langle v, \sgn (v) \rangle \leq  \| v \|_2  \| \sgn (v) \|_2 
     
    
    Since  :math:`\sgn(v)` can have at most  :math:`n` non-zero values, each with magnitude 1,
    
    .. math::
    
        \| \sgn (v) \|_2^2 \leq n \implies \| \sgn (v) \|_2 \leq \sqrt{n}.
    
    Thus, we get
    
    .. math::
    
        \| v \|_1  \leq \sqrt{n} \| v \|_2.
    


.. _res:ssm:l2_upper_bound_max_norm:

.. lemma::


    
    Let  :math:`v \in \mathbb{C}^n`. Then
    
    .. math::
    
        \| v \|_2 \leq \sqrt{n} \| v \|_{\infty}
    


.. proof::

    
    .. math::
    
        \| v \|_2^2 = \sum_{i=1}^n | v_i |^2 \leq n \underset{1 \leq i \leq n}{\max} ( | v_i |^2) = n \| v \|_{\infty}^2.
    
    Thus
    
    .. math::
    
        \| v \|_2 \leq \sqrt{n} \| v \|_{\infty}.
    

.. _res:ssm:p_q_norm_bounds:

.. lemma::


    
    Let  :math:`v \in \mathbb{C}^n`. Let  :math:`1 \leq p, q \leq \infty`.
    Then
    
    .. math::
    
        \| v \|_q \leq \| v \|_p \text{ whenever } p \leq q.
    

.. _res:ssm:one_vec_l1_norm:

.. lemma::


    
    Let  :math:`\OneVec \in \mathbb{C}^n` be the vector of all ones i.e.  :math:`\OneVec = (1, \dots, 1)`.
    Let  :math:`v \in \mathbb{C}^n` be some arbitrary vector. Let  :math:`| v |` denote the vector of
    absolute values of entries in  :math:`v`. i.e.  :math:`|v|_i = |v_i| \Forall 1 \leq i \leq n`. Then
    
    .. math::
    
        \| v \|_1 = \OneVec^T | v | = \OneVec^H | v |.
     


.. proof::

    
    .. math::
    
        \OneVec^T | v | = \sum_{i=1}^n  | v |_i =   \sum_{i=1}^n  | v_i | = \| v \|_1.
    
    Finally since  :math:`\OneVec` consists only of real entries, hence its transpose and Hermitian 
    transpose are same.


.. _res:ssm:ones_matrix_l1_norm:

.. lemma::

    Let  :math:`\OneMat \in \CC^{n \times n}` be a square matrix of all ones. Let  :math:`v \in \mathbb{C}^n` 
    be some arbitrary vector. Then

    
    
    .. math::
    
        |v|^T \OneMat | v | = \| v \|_1^2.
    


.. proof::

    We know that
    
    .. math::
    
        \OneMat = \OneVec \OneVec^T
    
    Thus,
    
    .. math::
    
        |v|^T \OneMat | v |  = |v|^T  \OneVec \OneVec^T | v |  = (\OneVec^T | v | )^T \OneVec^T | v | =  \| v \|_1 \| v \|_1 = \| v \|_1^2.
    
    We used the fact that  :math:`\| v \|_1 = \OneVec^T | v |`.



.. _res:ssm:k_th_largest_entry_l1_norm:
.. _eq:ssm:k_th_largest_entry_l1_norm:

.. theorem::


    
     :math:`k`-th largest (magnitude) entry in a vector  :math:`x \in \mathbb{C}^n` denoted by  :math:`x_{(k)}` obeys

    
    .. math::
    
    
        
        | x_{(k)} | \leq  \frac{\| x \|_1}{k}
    


.. proof::

    Let  :math:`n_1, n_2, \dots, n_N` be a permutation of  :math:`\{ 1, 2, \dots, n \}` such that
    
    .. math::
    
        |x_{n_1} | \geq  | x_{n_2} | \geq \dots \geq  | x_{n_N} |.
    
    Thus, the  :math:`k`-th largest entry in  :math:`x` is  :math:`x_{n_k}`. It is clear that
    
    .. math::
    
        \| x \|_1 = \sum_{i=1}^n | x_i | = \sum_{i=1}^n |x_{n_i} |
    
    
    Obviously
    
    .. math::
    
        |x_{n_1} | \leq \sum_{i=1}^n |x_{n_i} | = \| x \|_1.
    
    Similarly
    
    .. math::
    
        k |x_{n_k} | = |x_{n_k} | + \dots + |x_{n_k} |  \leq |x_{n_1} | + \dots + |x_{n_k} | \leq \sum_{i=1}^n |x_{n_i} | \leq  \| x \|_1.
    
    Thus
    
    .. math::
    
        |x_{n_k} |  \leq \frac{\| x \|_1}{k}.
Linear Algebra
====================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   vector_norms
   sparse_vectors
Sparse Vectors
========================


In this section we explore some useful properties of  :math:`\Sigma_k`, the set of  :math:`k`-sparse signals in standard basis
for  :math:`\mathbb{C}^n`.

We recall that

.. math::

    \Sigma_k  = \{ x \in \mathbb{C}^n : \| x \|_0 \leq k \}.


This set is a union of  :math:`\binom{n}{k}` subspaces of  :math:`\mathbb{C}^n` each of which
is is constructed by an index set  :math:`\Lambda \subset \{1, \dots, n \}` with  :math:`| \Lambda | = k` choosing
:math:`k` specific dimensions of  :math:`\mathbb{C}^n`. 

We first present some lemmas which connect the  :math:`l_1`,  :math:`l_2` and  :math:`l_{\infty}` norms of vectors
in  :math:`\Sigma_k`.

.. _lem:u_sigma_k_norms:

.. lemma::


    
    Suppose  :math:`u \in \Sigma_k`.  Then
    
    .. math::
    
          \frac{\| u\|_1}{\sqrt{k}} \leq \| u \|_2 \leq \sqrt{k} \| u \|_{\infty}.
    



.. proof::

   We can write  :math:`l_1` norm as
    
    .. math::
    
        \| u \|_1 = \langle u, \sgn (u) \rangle.
    
    
    By Cauchy-Schwartz inequality we have
    
    .. math::
    
        \langle u, \sgn (u) \rangle \leq  \| u \|_2  \| \sgn (u) \|_2 
     
    
    Since  :math:`u \in \Sigma_k`,  :math:`\sgn(u)` can have at most  :math:`k` non-zero values each with magnitude 1.
    Thus, we have
    
    .. math::
    
        \| \sgn (u) \|_2^2 \leq k \implies \| \sgn (u) \|_2 \leq \sqrt{k}
    
    
    Thus we get the lower bound
    
    .. math::
    
        \| u \|_1 \leq \| u \|_2 \sqrt{k}
        \implies \frac{\| u \|_1}{\sqrt{k}} \leq \| u \|_2.
    
    
    Now  :math:`| u_i | \leq \max(| u_i |) = \| u \|_{\infty}`. So we have
      
    .. math::
    
          \| u \|_2^2 = \sum_{i= 1}^{n} | u_i |^2 \leq  k \| u \|_{\infty}^2
    
    since there are only  :math:`k` non-zero terms in the expansion of  :math:`\| u \|_2^2`.
    
    This establishes the upper bound:
    
    .. math::
    
          \| u \|_2 \leq \sqrt{k} \| u \|_{\infty}
    


Tutorials
=====================================

.. toctree::
   :maxdepth: 3
   :caption: Contents:

   dirac_cosine_dictionaries
   admm_l1
