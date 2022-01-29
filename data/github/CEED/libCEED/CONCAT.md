# libCEED: Efficient Extensible Discretization

[![GitHub Actions][github-badge]][github-link]
[![GitLab-CI][gitlab-badge]][gitlab-link]
[![Azure Pipelines][azure-badge]][azure-link]
[![Code coverage][codecov-badge]][codecov-link]
[![BSD-2-Clause][license-badge]][license-link]
[![Documentation][doc-badge]][doc-link]
[![JOSS paper][joss-badge]][joss-link]
[![Binder][binder-badge]][binder-link]

## Summary and Purpose

libCEED provides fast algebra for element-based discretizations, designed for
performance portability, run-time flexibility, and clean embedding in higher
level libraries and applications. It offers a C99 interface as well as bindings
for Fortran, Python, Julia, and Rust.
While our focus is on high-order finite elements, the approach is mostly
algebraic and thus applicable to other discretizations in factored form, as
explained in the [user manual](https://libceed.readthedocs.io/en/latest/) and
API implementation portion of the
[documentation](https://libceed.readthedocs.io/en/latest/api/).

One of the challenges with high-order methods is that a global sparse matrix is
no longer a good representation of a high-order linear operator, both with
respect to the FLOPs needed for its evaluation, as well as the memory transfer
needed for a matvec.  Thus, high-order methods require a new "format" that still
represents a linear (or more generally non-linear) operator, but not through a
sparse matrix.

The goal of libCEED is to propose such a format, as well as supporting
implementations and data structures, that enable efficient operator evaluation
on a variety of computational device types (CPUs, GPUs, etc.). This new operator
description is based on algebraically
[factored form](https://libceed.readthedocs.io/en/latest/libCEEDapi/#finite-element-operator-decomposition),
which is easy to incorporate in a wide variety of applications, without significant
refactoring of their own discretization infrastructure.

The repository is part of the
[CEED software suite](http://ceed.exascaleproject.org/software/), a collection of
software benchmarks, miniapps, libraries and APIs for efficient exascale
discretizations based on high-order finite element and spectral element methods.
See <http://github.com/ceed> for more information and source code availability.

The CEED research is supported by the
[Exascale Computing Project](https://exascaleproject.org/exascale-computing-project)
(17-SC-20-SC), a collaborative effort of two U.S. Department of Energy
organizations (Office of Science and the National Nuclear Security
Administration) responsible for the planning and preparation of a
[capable exascale ecosystem](https://exascaleproject.org/what-is-exascale), including
software, applications, hardware, advanced system engineering and early testbed
platforms, in support of the nation’s exascale computing imperative.

For more details on the CEED API see the [user manual](https://libceed.readthedocs.io/en/latest/).

% gettingstarted-inclusion-marker

## Building

The CEED library, `libceed`, is a C99 library with no required dependencies, and
with Fortran, Python, Julia, and Rust interfaces.  It can be built using:

```
make
```

or, with optimization flags:

```
make OPT='-O3 -march=skylake-avx512 -ffp-contract=fast'
```

These optimization flags are used by all languages (C, C++, Fortran) and this
makefile variable can also be set for testing and examples (below).

The library attempts to automatically detect support for the AVX
instruction set using gcc-style compiler options for the host.
Support may need to be manually specified via:

```
make AVX=1
```

or:

```
make AVX=0
```

if your compiler does not support gcc-style options, if you are cross
compiling, etc.

To enable CUDA support, add `CUDA_DIR=/opt/cuda` or an appropriate directory
to your `make` invocation. To enable HIP support, add `HIP_DIR=/opt/rocm` or
an appropriate directory. To store these or other arguments as defaults for
future invocations of `make`, use:

```
make configure CUDA_DIR=/usr/local/cuda HIP_DIR=/opt/rocm OPT='-O3 -march=znver2'
```

which stores these variables in `config.mk`.

## Additional Language Interfaces

The Fortran interface is built alongside the library automatically.

Python users can install using:

```
pip install libceed
```

or in a clone of the repository via `pip install .`.

Julia users can install using:

```
$ julia
julia> ]
pkg> add LibCEED
```

See the [LibCEED.jl documentation](http://ceed.exascaleproject.org/libCEED-julia-docs/dev/)
for more information.

Rust users can include libCEED via `Cargo.toml`:

```toml
[dependencies]
libceed = { git = "https://github.com/CEED/libCEED", branch = "main" }
```

See the [Cargo documentation](https://doc.rust-lang.org/cargo/reference/specifying-dependencies.html#specifying-dependencies-from-git-repositories) for details.

## Testing

The test suite produces [TAP](https://testanything.org) output and is run by:

```
make test
```

or, using the `prove` tool distributed with Perl (recommended):

```
make prove
```

## Backends

There are multiple supported backends, which can be selected at runtime in the examples:

| CEED resource              | Backend                                           | Deterministic Capable |
| :---                       | :---                                              | :---:                 |
||
| **CPU Native**             |
| `/cpu/self/ref/serial`     | Serial reference implementation                   | Yes                   |
| `/cpu/self/ref/blocked`    | Blocked reference implementation                  | Yes                   |
| `/cpu/self/opt/serial`     | Serial optimized C implementation                 | Yes                   |
| `/cpu/self/opt/blocked`    | Blocked optimized C implementation                | Yes                   |
| `/cpu/self/avx/serial`     | Serial AVX implementation                         | Yes                   |
| `/cpu/self/avx/blocked`    | Blocked AVX implementation                        | Yes                   |
||
| **CPU Valgrind**           |
| `/cpu/self/memcheck/*`     | Memcheck backends, undefined value checks         | Yes                   |
||
| **CPU LIBXSMM**            |
| `/cpu/self/xsmm/serial`    | Serial LIBXSMM implementation                     | Yes                   |
| `/cpu/self/xsmm/blocked`   | Blocked LIBXSMM implementation                    | Yes                   |
||
| **CUDA Native**            |
| `/gpu/cuda/ref`            | Reference pure CUDA kernels                       | Yes                   |
| `/gpu/cuda/shared`         | Optimized pure CUDA kernels using shared memory   | Yes                   |
| `/gpu/cuda/gen`            | Optimized pure CUDA kernels using code generation | No                    |
||
| **HIP Native**             |
| `/gpu/hip/ref`             | Reference pure HIP kernels                        | Yes                   |
| `/gpu/hip/shared`          | Optimized pure HIP kernels using shared memory    | Yes                   |
| `/gpu/hip/gen`             | Optimized pure HIP kernels using code generation  | No                    |
||
| **MAGMA**                  |
| `/gpu/cuda/magma`          | CUDA MAGMA kernels                                | No                    |
| `/gpu/cuda/magma/det`      | CUDA MAGMA kernels                                | Yes                   |
| `/gpu/hip/magma`           | HIP MAGMA kernels                                 | No                    |
| `/gpu/hip/magma/det`       | HIP MAGMA kernels                                 | Yes                   |
||
| **OCCA**                   |
| `/*/occa`                  | Selects backend based on available OCCA modes     | Yes                   |
| `/cpu/self/occa`           | OCCA backend with serial CPU kernels              | Yes                   |
| `/cpu/openmp/occa`         | OCCA backend with OpenMP kernels                  | Yes                   |
| `/gpu/cuda/occa`           | OCCA backend with CUDA kernels                    | Yes                   |
| `/gpu/hip/occa`~           | OCCA backend with HIP kernels                     | Yes                   |

The `/cpu/self/*/serial` backends process one element at a time and are intended for meshes
with a smaller number of high order elements. The `/cpu/self/*/blocked` backends process
blocked batches of eight interlaced elements and are intended for meshes with higher numbers
of elements.

The `/cpu/self/ref/*` backends are written in pure C and provide basic functionality.

The `/cpu/self/opt/*` backends are written in pure C and use partial e-vectors to improve performance.

The `/cpu/self/avx/*` backends rely upon AVX instructions to provide vectorized CPU performance.

The `/cpu/self/memcheck/*` backends rely upon the [Valgrind](http://valgrind.org/) Memcheck tool
to help verify that user QFunctions have no undefined values. To use, run your code with
Valgrind and the Memcheck backends, e.g. `valgrind ./build/ex1 -ceed /cpu/self/ref/memcheck`. A
'development' or 'debugging' version of Valgrind with headers is required to use this backend.
This backend can be run in serial or blocked mode and defaults to running in the serial mode
if `/cpu/self/memcheck` is selected at runtime.

The `/cpu/self/xsmm/*` backends rely upon the [LIBXSMM](http://github.com/hfp/libxsmm) package
to provide vectorized CPU performance. If linking MKL and LIBXSMM is desired but
the Makefile is not detecting `MKLROOT`, linking libCEED against MKL can be
forced by setting the environment variable `MKL=1`.

The `/gpu/cuda/*` backends provide GPU performance strictly using CUDA.

The `/gpu/hip/*` backends provide GPU performance strictly using HIP. They are based on
the `/gpu/cuda/*` backends.  ROCm version 3.6 or newer is required.

The `/gpu/*/magma/*` backends rely upon the [MAGMA](https://bitbucket.org/icl/magma) package.
To enable the MAGMA backends, the environment variable `MAGMA_DIR` must point to the top-level
MAGMA directory, with the MAGMA library located in `$(MAGMA_DIR)/lib/`.
By default, `MAGMA_DIR` is set to `../magma`; to build the MAGMA backends
with a MAGMA installation located elsewhere, create a link to `magma/` in libCEED's parent
directory, or set `MAGMA_DIR` to the proper location.  MAGMA version 2.5.0 or newer is required.
Currently, each MAGMA library installation is only built for either CUDA or HIP.  The corresponding
set of libCEED backends (`/gpu/cuda/magma/*` or `/gpu/hip/magma/*`) will automatically be built
for the version of the MAGMA library found in `MAGMA_DIR`.

Users can specify a device for all CUDA, HIP, and MAGMA backends through adding `:device_id=#`
after the resource name.  For example:

> - `/gpu/cuda/gen:device_id=1`

The `/*/occa` backends rely upon the [OCCA](http://github.com/libocca/occa) package to provide
cross platform performance. To enable the OCCA backend, the environment variable `OCCA_DIR` must point
to the top-level OCCA directory, with the OCCA library located in the `${OCCA_DIR}/lib` (By default,
`OCCA_DIR` is set to `../occa`).

Additionally, users can pass specific OCCA device properties after setting the CEED resource.
For example:

> - `"/*/occa:mode='CUDA',device_id=0"`

Bit-for-bit reproducibility is important in some applications.
However, some libCEED backends use non-deterministic operations, such as `atomicAdd` for increased performance.
The backends which are capable of generating reproducible results, with the proper compilation options, are highlighted in the list above.

## Examples

libCEED comes with several examples of its usage, ranging from standalone C
codes in the `/examples/ceed` directory to examples based on external packages,
such as MFEM, PETSc, and Nek5000. Nek5000 v18.0 or greater is required.

To build the examples, set the `MFEM_DIR`, `PETSC_DIR`, and
`NEK5K_DIR` variables and run:

```
cd examples/
```

% running-examples-inclusion-marker

```console
# libCEED examples on CPU and GPU
cd ceed/
make
./ex1-volume -ceed /cpu/self
./ex1-volume -ceed /gpu/cuda
./ex2-surface -ceed /cpu/self
./ex2-surface -ceed /gpu/cuda
cd ..

# MFEM+libCEED examples on CPU and GPU
cd mfem/
make
./bp1 -ceed /cpu/self -no-vis
./bp3 -ceed /gpu/cuda -no-vis
cd ..

# Nek5000+libCEED examples on CPU and GPU
cd nek/
make
./nek-examples.sh -e bp1 -ceed /cpu/self -b 3
./nek-examples.sh -e bp3 -ceed /gpu/cuda -b 3
cd ..

# PETSc+libCEED examples on CPU and GPU
cd petsc/
make
./bps -problem bp1 -ceed /cpu/self
./bps -problem bp2 -ceed /gpu/cuda
./bps -problem bp3 -ceed /cpu/self
./bps -problem bp4 -ceed /gpu/cuda
./bps -problem bp5 -ceed /cpu/self
./bps -problem bp6 -ceed /gpu/cuda
cd ..

cd petsc/
make
./bpsraw -problem bp1 -ceed /cpu/self
./bpsraw -problem bp2 -ceed /gpu/cuda
./bpsraw -problem bp3 -ceed /cpu/self
./bpsraw -problem bp4 -ceed /gpu/cuda
./bpsraw -problem bp5 -ceed /cpu/self
./bpsraw -problem bp6 -ceed /gpu/cuda
cd ..

cd petsc/
make
./bpssphere -problem bp1 -ceed /cpu/self
./bpssphere -problem bp2 -ceed /gpu/cuda
./bpssphere -problem bp3 -ceed /cpu/self
./bpssphere -problem bp4 -ceed /gpu/cuda
./bpssphere -problem bp5 -ceed /cpu/self
./bpssphere -problem bp6 -ceed /gpu/cuda
cd ..

cd petsc/
make
./area -problem cube -ceed /cpu/self -degree 3
./area -problem cube -ceed /gpu/cuda -degree 3
./area -problem sphere -ceed /cpu/self -degree 3 -dm_refine 2
./area -problem sphere -ceed /gpu/cuda -degree 3 -dm_refine 2

cd fluids/
make
./navierstokes -ceed /cpu/self -degree 1
./navierstokes -ceed /gpu/cuda -degree 1
cd ..

cd solids/
make
./elasticity -ceed /cpu/self -mesh [.exo file] -degree 2 -E 1 -nu 0.3 -problem Linear -forcing mms
./elasticity -ceed /gpu/cuda -mesh [.exo file] -degree 2 -E 1 -nu 0.3 -problem Linear -forcing mms
cd ..
```

For the last example shown, sample meshes to be used in place of
`[.exo file]` can be found at <https://github.com/jeremylt/ceedSampleMeshes>

The above code assumes a GPU-capable machine with the OCCA backend
enabled. Depending on the available backends, other CEED resource
specifiers can be provided with the `-ceed` option. Other command line
arguments can be found in [examples/petsc](https://github.com/CEED/libCEED/blob/main/examples/petsc/README.md).

% benchmarks-marker

## Benchmarks

A sequence of benchmarks for all enabled backends can be run using:

```
make benchmarks
```

The results from the benchmarks are stored inside the `benchmarks/` directory
and they can be viewed using the commands (requires python with matplotlib):

```
cd benchmarks
python postprocess-plot.py petsc-bps-bp1-*-output.txt
python postprocess-plot.py petsc-bps-bp3-*-output.txt
```

Using the `benchmarks` target runs a comprehensive set of benchmarks which may
take some time to run. Subsets of the benchmarks can be run using the scripts in the `benchmarks` folder.

For more details about the benchmarks, see the `benchmarks/README.md` file.

## Install

To install libCEED, run:

```
make install prefix=/usr/local
```

or (e.g., if creating packages):

```
make install prefix=/usr DESTDIR=/packaging/path
```

The usual variables like `CC` and `CFLAGS` are used, and optimization flags
for all languages can be set using the likes of `OPT='-O3 -march=native'`. Use
`STATIC=1` to build static libraries (`libceed.a`).

To install libCEED for Python, run:

```
pip install libceed
```

with the desired setuptools options, such as `--user`.

### pkg-config

In addition to library and header, libCEED provides a [pkg-config](https://en.wikipedia.org/wiki/Pkg-config)
file that can be used to easily compile and link.
[For example](https://people.freedesktop.org/~dbn/pkg-config-guide.html#faq), if
`$prefix` is a standard location or you set the environment variable
`PKG_CONFIG_PATH`:

```
cc `pkg-config --cflags --libs ceed` -o myapp myapp.c
```

will build `myapp` with libCEED.  This can be used with the source or
installed directories.  Most build systems have support for pkg-config.

## Contact

You can reach the libCEED team by emailing [ceed-users@llnl.gov](mailto:ceed-users@llnl.gov)
or by leaving a comment in the [issue tracker](https://github.com/CEED/libCEED/issues).

## How to Cite

If you utilize libCEED please cite:

```
@article{libceed-joss-paper,
  author       = {Jed Brown and Ahmad Abdelfattah and Valeria Barra and Natalie Beams and Jean Sylvain Camier and Veselin Dobrev and Yohann Dudouit and Leila Ghaffari and Tzanio Kolev and David Medina and Will Pazner and Thilina Ratnayaka and Jeremy Thompson and Stan Tomov},
  title        = {{libCEED}: Fast algebra for high-order element-based discretizations},
  journal      = {Journal of Open Source Software},
  year         = {2021},
  publisher    = {The Open Journal},
  volume       = {6},
  number       = {63},
  pages        = {2945},
  doi          = {10.21105/joss.02945}
}

@misc{libceed-user-manual,
  author       = {Abdelfattah, Ahmad and
                  Barra, Valeria and
                  Beams, Natalie and
                  Brown, Jed and
                  Camier, Jean-Sylvain and
                  Dobrev, Veselin and
                  Dudouit, Yohann and
                  Ghaffari, Leila and
                  Kolev, Tzanio and
                  Medina, David and
                  Pazner, Will and
                  Ratnayaka, Thilina and
                  Thompson, Jeremy L and
                  Tomov, Stanimire},
  title        = {{libCEED} User Manual},
  month        = jul,
  year         = 2021,
  publisher    = {Zenodo},
  version      = {0.9.0},
  doi          = {10.5281/zenodo.5077489}
}
```

For libCEED's Python interface please cite:

```
@InProceedings{libceed-paper-proc-scipy-2020,
  author    = {{V}aleria {B}arra and {J}ed {B}rown and {J}eremy {T}hompson and {Y}ohann {D}udouit},
  title     = {{H}igh-performance operator evaluations with ease of use: lib{C}{E}{E}{D}'s {P}ython interface},
  booktitle = {{P}roceedings of the 19th {P}ython in {S}cience {C}onference},
  pages     = {85 - 90},
  year      = {2020},
  editor    = {{M}eghann {A}garwal and {C}hris {C}alloway and {D}illon {N}iederhut and {D}avid {S}hupe},
  doi       = {10.25080/Majora-342d178e-00c}
}
```

The BiBTeX entries for these references can be found in the
`doc/bib/references.bib` file.

## Copyright

The following copyright applies to each file in the CEED software suite, unless
otherwise stated in the file:

> Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at the
> Lawrence Livermore National Laboratory. LLNL-CODE-734707. All Rights reserved.

See files LICENSE and NOTICE for details.

[github-badge]: https://github.com/CEED/libCEED/workflows/C/Fortran/badge.svg
[github-link]: https://github.com/CEED/libCEED/actions
[gitlab-badge]: https://gitlab.com/libceed/libCEED/badges/main/pipeline.svg?key_text=GitLab-CI
[gitlab-link]: https://gitlab.com/libceed/libCEED/-/pipelines?page=1&scope=all&ref=main
[azure-badge]: https://dev.azure.com/CEED-ECP/libCEED/_apis/build/status/CEED.libCEED?branchName=main
[azure-link]: https://dev.azure.com/CEED-ECP/libCEED/_build?definitionId=2
[codecov-badge]: https://codecov.io/gh/CEED/libCEED/branch/main/graphs/badge.svg
[codecov-link]: https://codecov.io/gh/CEED/libCEED/
[license-badge]: https://img.shields.io/badge/License-BSD%202--Clause-orange.svg
[license-link]: https://opensource.org/licenses/BSD-2-Clause
[doc-badge]: https://readthedocs.org/projects/libceed/badge/?version=latest
[doc-link]: https://libceed.readthedocs.io/en/latest/?badge=latest
[joss-badge]: https://joss.theoj.org/papers/10.21105/joss.02945/status.svg
[joss-link]: https://doi.org/10.21105/joss.02945
[binder-badge]: http://mybinder.org/badge_logo.svg
[binder-link]: https://mybinder.org/v2/gh/CEED/libCEED/main?urlpath=lab/tree/examples/python/tutorial-0-ceed.ipynb
# libCEED Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, religion, or sexual identity
and orientation.

We pledge to act and interact in ways that contribute to an open, welcoming,
diverse, inclusive, and healthy community.

## Our Standards

Examples of behavior that contributes to a positive environment for our
community include:

* Demonstrating empathy and kindness toward other people
* Being respectful of differing opinions, viewpoints, and experiences
* Giving and gracefully accepting constructive feedback
* Accepting responsibility and apologizing to those affected by our mistakes,
  and learning from the experience
* Focusing on what is best not just for us as individuals, but for the
  overall community

Examples of unacceptable behavior include:

* The use of sexualized language or imagery, and sexual attention or
  advances of any kind
* Trolling, insulting or derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or email
  address, without their explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

## Enforcement Responsibilities

Community leaders are responsible for clarifying and enforcing our standards of
acceptable behavior and will take appropriate and fair corrective action in
response to any behavior that they deem inappropriate, threatening, offensive,
or harmful.

Community leaders have the right and responsibility to remove, edit, or reject
comments, commits, code, wiki edits, issues, and other contributions that are
not aligned to this Code of Conduct, and will communicate reasons for moderation
decisions when appropriate.

## Scope

This Code of Conduct applies within all community spaces, and also applies when
an individual is officially representing the community in public spaces.
Examples of representing our community include using an official e-mail address,
posting via an official social media account, or acting as an appointed
representative at an online or offline event.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported to the community leaders responsible for enforcement at
jed@jedbrown.org, valeria.barra@colorado.edu, or tzanio@llnl.gov.
All complaints will be reviewed and investigated promptly and fairly.

All community leaders are obligated to respect the privacy and security of the
reporter of any incident.

## Enforcement Guidelines

Community leaders will follow these Community Impact Guidelines in determining
the consequences for any action they deem in violation of this Code of Conduct:

### 1. Correction

**Community Impact**: Use of inappropriate language or other behavior deemed
unprofessional or unwelcome in the community.

**Consequence**: A private, written warning from community leaders, providing
clarity around the nature of the violation and an explanation of why the
behavior was inappropriate. A public apology may be requested.

### 2. Warning

**Community Impact**: A violation through a single incident or series
of actions.

**Consequence**: A warning with consequences for continued behavior. No
interaction with the people involved, including unsolicited interaction with
those enforcing the Code of Conduct, for a specified period of time. This
includes avoiding interactions in community spaces as well as external channels
like social media. Violating these terms may lead to a temporary or
permanent ban.

### 3. Temporary Ban

**Community Impact**: A serious violation of community standards, including
sustained inappropriate behavior.

**Consequence**: A temporary ban from any sort of interaction or public
communication with the community for a specified period of time. No public or
private interaction with the people involved, including unsolicited interaction
with those enforcing the Code of Conduct, is allowed during this period.
Violating these terms may lead to a permanent ban.

### 4. Permanent Ban

**Community Impact**: Demonstrating a pattern of violation of community
standards, including sustained inappropriate behavior,  harassment of an
individual, or aggression toward or disparagement of classes of individuals.

**Consequence**: A permanent ban from any sort of public interaction within
the community.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 2.0, available at
https://www.contributor-covenant.org/version/2/0/code_of_conduct.html.

Community Impact Guidelines were inspired by [Mozilla's code of conduct
enforcement ladder](https://github.com/mozilla/diversity).

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see the FAQ at
https://www.contributor-covenant.org/faq. Translations are available at
https://www.contributor-covenant.org/translations.
# libCEED: How to Contribute

Contributions to libCEED are encouraged.
<!---
Please use a pull request to the appropriate branch ('stable' for
backward-compatible bug fixes for the last stable release, main' for
new features and everything else).
-->
Please make your commits well-organized and
[atomic](https://en.wikipedia.org/wiki/Atomic_commit#Atomic_commit_convention),
using `git rebase --interactive` as needed.  Check that tests
(including "examples") pass using `make prove-all`.  If adding a new
feature, please add or extend a test so that your new feature is
tested.

In typical development, every commit should compile, be covered by the
test suite, and pass all tests.  This improves the efficiency of
reviewing and facilitates use of
[`git bisect`](https://git-scm.com/docs/git-bisect).

Open an issue or RFC (request for comments) pull request to discuss
any significant changes before investing time.  It is useful to create
a WIP (work in progress) pull request for any long-running development
so that others can be aware of your work and help to avoid creating
merge conflicts.

Write commit messages for a reviewer of your pull request and for a
future developer (maybe you) that bisects and finds that a bug was
introduced in your commit.  The assumptions that are clear in your
mind while committing are likely not in the mind of whomever (possibly
you) needs to understand it in the future.

Give credit where credit is due using tags such as `Reported-by:
Helpful User <helpful@example.com>` or
[`Co-authored-by: Snippet Mentor <code.by@comment.com>`](https://help.github.com/en/github/committing-changes-to-your-project/creating-a-commit-with-multiple-authors#creating-co-authored-commits-on-the-command-line).
Please use a real name and email for your author information (`git
config user.name` and `user.email`).  If your author information or
email becomes inconsistent (look at `git shortlog -se`), please edit
`.mailmap` to obtain your preferred name and email address.

When contributors make a major contribution and support it, their names
are included in the automatically generated user-manual documentation.

Please avoid "merging from upstream" (like merging 'main' into your
feature branch) unless there is a specific reason to do so, in which
case you should explain why in the merge commit.
[Rationale](https://lwn.net/Articles/328436/) from
[Junio](https://gitster.livejournal.com/42247.html) and
[Linus](http://yarchive.net/comp/linux/git_merges_from_upstream.html).

You can use `make style` to help conform to coding conventions of the
project, but try to avoid mixing whitespace or formatting changes with
content changes (see atomicity above).

By submitting a pull request, you are affirming the following.

## [Developer's Certificate of Origin 1.1](https://developercertificate.org/)

By making a contribution to this project, I certify that:

(a) The contribution was created in whole or in part by me and I
    have the right to submit it under the open source license
    indicated in the file; or

(b) The contribution is based upon previous work that, to the best
    of my knowledge, is covered under an appropriate open source
    license and I have the right under that license to submit that
    work with modifications, whether created in whole or in part
    by me, under the same open source license (unless I am
    permitted to submit under a different license), as indicated
    in the file; or

(c) The contribution was provided directly to me by some other
    person who certified (a), (b) or (c) and I have not modified
    it.

(d) I understand and agree that this project and the contribution
    are public and that a record of the contribution (including all
    personal information I submit with it, including my sign-off) is
    maintained indefinitely and may be redistributed consistent with
    this project or the open source license(s) involved.

## Authorship

libCEED contains components authored by many individuals.  It is
important that contributors receive appropriate recognition through
informal and academically-recognized credit systems such as
publications.  Status as a named author on the users manual and
libCEED software publications will be granted for those who

1. make significant contributions to libCEED (in implementation,
  documentation, conceptualization, review, etc.) and
2. maintain and support those contributions.

Maintainers will do their best to notice when contributions reach this
level and add your name to `AUTHORS`, but please email or create an
issue if you believe your contributions have met these criteria and
haven't yet been acknowledged.

Authors of publications about libCEED as a whole, including
DOI-bearing archives, shall offer co-authorship to all individuals
listed in the `AUTHORS` file.  Authors of publications claiming
specific libCEED contributions shall evaluate those listed in
`AUTHORS` and offer co-authorship to those who made significant
intellectual contributions to the work.

Note that there is no co-authorship expectation for those publishing
about use of libCEED (versus creation of new features in libCEED), but
see the [citing section](https://libceed.readthedocs.io/en/latest/gettingstarted/#how-to-cite)
and use your judgment regarding significance of support/advice you may have
received in developing your use case and interpreting results.
# Release Procedures

*These notes are meant for a maintainer to create official releases.*

In preparing a release, create a branch to hold pre-release commits. We ideally want all release mechanics (for all languages) to be in one commit, which will then be tagged. (This will change if/when we stop synchronizing releases across all language bindings.)

## Core C library

Some minor bookkeeping updates are needed when releasing a new version of the core library.

The version number must be updated in

* `include/ceed/ceed.h`
* `ceed.pc.template`
* `Doxyfile`
* `CITATION.cff`

as well as `include/ceed/ceed.h` (`CEED_VERSION_MAJOR`, `CEED_VERSION_MINOR`).

Additionally, the release notes in `doc/sphinx/source/releasenotes.rst` should be updated. Use `git log --first-parent v0.7..` to get a sense of the pull requests that have been merged and thus might warrant emphasizing in the release notes. While doing this, gather a couple sentences for key features to highlight on [GitHub releases](https://github.com/CEED/libCEED/releases). The "Current Main" heading needs to be named for the release.

Use `make doc-latexpdf` to build a PDF users manual and inspect it for missing references or formatting problems (e.g., with images that were converted to PDF). This contains the same content as the website, but will be archived on Zenodo.

### Quality control and good citizenry

1. If making a minor release, check for API and ABI changes that could break [semantic versioning](https://semver.org/). The [ABI compliance checker](https://github.com/lvc/abi-compliance-checker) is a useful tool, as is `nm -D libceed.so` and checking for public symbols (capital letters like `T` and `D` that are not namespaced).

2. Double check testing on any architectures that may not be exercised in continuous integration (e.g., HPC facilities) and with users of libCEED, such as MFEM and PETSc applications. While unsupported changes do not prevent release, it's polite to make a PR to support the new release, and it's good for quality to test before taggin a libCEED release.

3. Update and test all the language bindings (see below) within the branch.

4. Check that `spack install libceed@develop` works prior to tagging. The Spack `libceed/package.py` file should be updated immediately after tagging a release.

### Tagging and releasing on GitHub

0. Confirm all the steps above, including all language bindings.
1. `git commit -am'libCEED 0.8.1'`
More frequently, this is amending the commit message on an in-progress commit, after rebasing if applicable on latest `main`.
2. `git push` updates the PR holding release; opportunity for others to review
3. `git switch main && git merge --ff-only HEAD@{1}` fast-forward merge into `main`
4. `git tag --sign -m'libCEED 0.8.1'`
5. `git push origin main v0.8.1`
6. Draft a [new release on GitHub](https://github.com/CEED/libCEED/releases), using a few sentences gathered from the release notes.
7. Submit a PR to Spack.
8. Publish Julia, Python, and Rust packages.

### Archive Users Manual on Zenodo

Generate the PDF using `make doc-latexpdf`, click "New version" on the [Zenodo
record](https://zenodo.org/record/4302737) and upload. Update author info if applicable (new
authors, or existing authors changing institutions). Make a new PR to update the version
number and DOI in `README.rst` and `doc/bib/references.bib`.

## Julia

libCEED's Julia interface (LibCEED.jl) has two components:

* LibCEED.jl, the user-facing package that contains the Julia interface.
* libCEED_jll, a binary wrapper package ("jll package") that contains prebuilt binaries of the
  libCEED library for various architectures.

When there is a new release of libCEED, both of these components need to be updated. First,
libCEED_jll is updated, and then LibCEED.jl.

### Updating libCEED_jll

The binary wrapper package libCEED_jll is updated by making a pull request against
[Yggdrasil](https://github.com/JuliaPackaging/Yggdrasil), the Julia community build tree. In this
PR, the file `L/libCEED/build_tarballs.jl` should be changed to update version number and change the
hash of the libCEED commit to use to build the binaries, similar to the following diff:
```diff
diff --git a/L/libCEED/build_tarballs.jl b/L/libCEED/build_tarballs.jl
--- a/L/libCEED/build_tarballs.jl
+++ b/L/libCEED/build_tarballs.jl
@@ -3,11 +3,11 @@
 using BinaryBuilder, Pkg

 name = "libCEED"
-version = v"0.7.0"
+version = v"0.8.0"

 # Collection of sources required to complete build
 sources = [
-    GitSource("https://github.com/CEED/libCEED.git", "06988bf74cc6ac18eacafe7930f080803395ba29")
+    GitSource("https://github.com/CEED/libCEED.git", "e8f234590eddcce2220edb1d6e979af7a3c35f82")
 ]
```
After the PR is merged into Yggdrasil, the new version of libCEED_jll will automatically be
registered, and then we can proceed to update LibCEED.jl.

### Updating LibCEED.jl

After the binary wrapper package libCEED_jll has been updated, we are ready to update the main Julia
interface LibCEED.jl. This requires updating the file `julia/LibCEED.jl/Project.toml` in the libCEED
repository. The version number should be incremented, and the dependency on the updated version of
libCEED_jll should be listed:
```diff
diff --git a/julia/LibCEED.jl/Project.toml b/julia/LibCEED.jl/Project.toml
--- a/julia/LibCEED.jl/Project.toml
+++ b/julia/LibCEED.jl/Project.toml
@@ -1,7 +1,7 @@
 name = "LibCEED"
 uuid = "2cd74e05-b976-4426-91fa-5f1011f8952b"
-version = "0.1.0"
+version = "0.1.1"

 [deps]
 CEnum = "fa961155-64e5-5f13-b03f-caf6b980ea82"
@@ -26,4 +26,4 @@ Cassette = "0.3"
 Requires = "1"
 StaticArrays = "0.12"
 UnsafeArrays = "1"
-libCEED_jll = "0.7"
+libCEED_jll = "0.8"
```
Once this change is merged into libCEED's `main` branch, the updated package version can be
registered using the GitHub registrator bot by commenting on the commit:

> @JuliaRegistrator register branch=main subdir=julia/LibCEED.jl

At this point, the bot should create a PR against the [general Julia
registry](https://github.com/JuliaRegistries/General), which should be merged automatically after a
short delay.

### Moving development tests to release tests

LibCEED.jl has both _development_ and _release_ unit tests. The _release_ tests are run both with
the current build of libCEED, and with the most recent release of libCEED_jll. The _development_
tests may use features which were not available in the most recent release, and so they are only run
with the current build of libCEED.

Upon release, the development tests may be moved to the release tests, so that these features will
be tested against the most recent release of libCEED_jll. The release tests are found in the file
`julia/LibCEED.jl/test/runtests.jl` and the development tests are found in
`julia/LibCEED.jl/test/rundevtests.jl`.

## Python

The Python package gets its version from `ceed.pc.template` so there are no file modifications necessary.

1. `make wheel` builds and tests the wheels using Docker. See the [manylinux repo](https://github.com/pypa/manylinux) for source and usage inforamtion. If this succeeds, the completed wheels are in `wheelhouse/libceed-0.8-cp39-cp39-manylinux2010_x86_64.whl`.
2. Manually test one or more of the wheels by creating a virtualenv and using `pip install wheelhouse/libceed-0.8-cp39-cp39-manylinux2010_x86_64.whl`, then `python -c 'import libceed'` or otherwise running tests.
3. Create a `~/.pypirc` with entries for `testpypi` (`https://test.pypi.org/legacy/`) and the real `pypi`.
4. Upload to `testpypi` using
```console
$ twine upload --repository testpypi wheelhouse/libceed-0.8-cp39-cp39-manylinux2010_x86_64.whl
```
5. Test installing on another machine/in a virtualenv:
```console
$ pip install --index-url https://test.pypi.org/simple --extra-index-url https://pypi.org/simple libceed
```
The `--extra-index-url` argument allows dependencies like `cffi` and `numpy` from being fetched from the non-test repository.
6. Do it live:
```console
$ twine upload --repository pypi wheelhouse/libceed-0.8-cp39-cp39-manylinux2010_x86_64.whl
```
Note that this cannot be amended.

## Rust

The Rust crates for libCEED are split into
1. [`libceed-sys`](https://crates.io/crates/libceed-sys), which handles building/finding the `libceed.so` or `libceed.a` library and providing unsafe Rust bindings (one to one with the C interface, using C FFI datatypes)
2. [`libceed`](https://crates.io/crates/libceed) containing the safe and idiomatic Rust bindings.

We currently apply the same version number across both of these crates. There are some tests for version strings matching, but in short, one needs to update the following locations.

```console
$ git grep '0\.8' -- rust/
rust/libceed-sys/Cargo.toml:version = "0.8.0"
rust/libceed-sys/README.md:libceed-sys = "0.8.0"
rust/libceed-sys/build.rs:        .atleast_version("0.8")
rust/libceed/Cargo.toml:version = "0.8.0"
rust/libceed/Cargo.toml:libceed-sys = { version = "0.8", path = "../libceed-sys" }
rust/libceed/README.md:libceed = "0.8.0"
```

After doing this,

1. `cargo package --list` to see that the file list makes sense.
2. `cargo package` to build crates locally
3. `cargo publish` to publish the crates to https://crates.io
# LibCEED.jl: Julia Interface for [libCEED](https://github.com/CEED/libCEED)

Please see the [LibCEED.jl
documentation](http://ceed.exascaleproject.org/libCEED-julia-docs/dev/) for
usage and API documentation.

## Installation

The LibCEED.jl package can be installed with Julia's package manager by running
`] add LibCEED`. This will automatically install a pre-built binary of the
libCEED library. If you require features of a specific build of libCEED (e.g.
CUDA/GPU support, specific compiler flags, etc.) then you should compile your
own version of the libCEED library, and configure LibCEED.jl to use this binary
as described in the [Configuring LibCEED.jl](#configuring-libceedjl) section.

**Warning:** the pre-built libCEED binaries do not support CUDA backends

The pre-built binaries automatically installed by LibCEED.jl (through the
[libCEED_jll](https://juliahub.com/ui/Packages/libCEED_jll/LB2fn) package) are
not built with CUDA support. If you want to run libCEED on the GPU, you will
have to build libCEED from source and configure LibCEED.jl as described in the
[Configuring LibCEED.jl](#configuring-libceedjl) section.

### Configuring LibCEED.jl

By default, LibCEED.jl will use the pre-built libCEED binaries provided by the
[libCEED_jll](https://juliahub.com/ui/Packages/libCEED_jll/LB2fn) package. If
you wish to use a different libCEED binary (e.g. one built from source),
LibCEED.jl can be configured using Julia's _preferences_ mechanism. Note that
this preference will be set for the currently active Julia environment, and can
be different between different environments. The Julia session must be restarted
for changes to take effect.

```julia
julia> using LibCEED
julia> set_libceed_path!("/path/to/libceed.so")
[ Info: Setting the libCEED library path to /path/to/libceed.so.
[ Info: Restart the Julia session for changes to take effect.
```

See [Preferences.jl](https://github.com/JuliaPackaging/Preferences.jl) for more
information.
# Linear Algebra

User Q-functions often perform small (1x1, 2x2, or 3x3) linear algebra
operations (determinant, matrix-vector product, etc.) at every Q-point. For good
performance, it is important to use specialized versions of these operations for
the given size.

If the matrix or vector is given in a statically sized container (e.g. using
[StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl/)) then this
happens automatically. However, if the matrix is not statically sized, and
instead is given as, for example, a view into a larger array, then LibCEED.jl
provides some convenient specialized functions.

In order to allow for generic code, the [`CeedDim`](@ref) struct is used for
dispatch. An object `D = CeedDim(dim)` can be created, and passed as a second
argument to functions like `det` to choose the specialized implementations. In
this case, `dim` should be known as a compile-time constant, otherwise it will
result in a type instability, and give poor performance.

For example:
```@repl
using LibCEED, LinearAlgebra

dim = 3;
J = rand(dim, dim);

det(J) # Slow!
det(J, CeedDim(dim)) # Fast!
```

```@docs
CeedDim
det(J, ::CeedDim{1})
setvoigt
getvoigt
```
# LibCEED.jl Docs

Documentation for the LibCEED.jl Julia interface to the
[libCEED](https://github.com/ceed/libceed) library.

For further information, see also the [libCEED
documentation](https://libceed.readthedocs.io/).

Several [short examples](Examples.md) are included to demonstrate the
functionality.

### Installation

The LibCEED.jl package can be installed with Julia's package manager by running
`] add LibCEED`. This will automatically install a pre-built binary of the
libCEED library. If you require features of a specific build of libCEED (e.g.
CUDA/GPU support, specific compiler flags, etc.) then you should compile your
own version of the libCEED library, and configure LibCEED.jl to use this binary
as described in the [Configuring LibCEED.jl](@ref) section.

!!! warning "The pre-built libCEED binaries do not support CUDA backends"
    The pre-built binaries automatically installed by LibCEED.jl (through the
    [libCEED_jll](https://juliahub.com/ui/Packages/libCEED_jll/LB2fn) package)
    are not built with CUDA support. If you want to run libCEED on the GPU, you
    will have to build libCEED from source and configure LibCEED.jl as described
    in the [Configuring LibCEED.jl](@ref) section.

#### Configuring LibCEED.jl

By default, LibCEED.jl will use the pre-built libCEED binaries provided by the
[libCEED_jll](https://juliahub.com/ui/Packages/libCEED_jll/LB2fn) package. If
you wish to use a different libCEED binary (e.g. one built from source),
LibCEED.jl can be configured using Julia's _preferences_ mechanism. Note that
this preference will be set for the currently active Julia environemnt, and can
be different between different environments. The Julia session must be restarted
for changes to take effect.

```julia
julia> using LibCEED
julia> set_libceed_path!("/path/to/libceed.so")
[ Info: Setting the libCEED library path to /path/to/libceed.so.
[ Info: Restart the Julia session for changes to take effect.
```

See the [library configuration documentation](LibCEED.md) for more details. For
information on Julia's preferences system, see
[Preferences.jl](https://github.com/JuliaPackaging/Preferences.jl).

### Features of the high-level interface for libCEED

#### User Q-functions

With LibCEED.jl, it is much easier to write dimension-independent user-defined
Q-functions that automatically work on the GPU. See the [related
documentation](UserQFunctions.md) for more information.

#### Safe access to CeedVector objects

When accessing [`CeedVector`](@ref) objects, the C interface requires the user
to manually call `CeedVectorGetArray`, paired with `CeedVectorRestoreArray`. If
the user wants read-only access, then the user must call
`CeedVectorGetArrayRead`, paired with `CeedVectorRestoreArrayRead`. This can
possibly be bug-prone, because the user may forget to restore the array, or may
match the `Read` version to get the array with non-`Read` version to restore the
array (or vice versa).

In LibCEED.jl, this difficulty is mitigated using the [`witharray`](@ref)
function and [`@witharray`](@ref) macro. There are also read-only versions,
[`witharray_read`](@ref) and [`@witharray_read`](@ref). When using this
functionality, it is impossible to forget to restore the array, and the correct
version is always paired properly.

For example, in `ex1-volume`, the following C code
```c
// Compute and print the sum of the entries of 'v' giving the mesh volume.
const CeedScalar *v_host;
CeedVectorGetArrayRead(v, CEED_MEM_HOST, &v_host);
CeedScalar vol = 0.;
for (CeedInt i = 0; i < sol_size; i++) {
  vol += v_host[i];
}
CeedVectorRestoreArrayRead(v, &v_host);
```
is replaced with the following equivalent Julia code
```julia
# Compute and print the sum of the entries of 'v' giving the mesh volume.
vol = witharray_read(sum, v, MEM_HOST)
```

In `ex2-surface`, the following C code
```c
// Initialize 'u' with sum of coordinates, x+y+z.
CeedScalar *u_host;
const CeedScalar *x_host;
CeedVectorGetArray(u, CEED_MEM_HOST, &u_host);
CeedVectorGetArrayRead(mesh_coords, CEED_MEM_HOST, &x_host);
for (CeedInt i = 0; i < sol_size; i++) {
  u_host[i] = 0;
  for (CeedInt d = 0; d < dim; d++)
    u_host[i] += x_host[i+d*sol_size];
}
CeedVectorRestoreArray(u, &u_host);
CeedVectorRestoreArrayRead(mesh_coords, &x_host);
```
is replaced with the following equivalent Julia code
```julia
@witharray_read(x_host=mesh_coords, size=(mesh_size÷dim, dim),
    @witharray(u_host=u, size=(sol_size,1),
        sum!(u_host, x_host)))
```
The macro version can provide better performance if a closure is required, and
allow for convenient reshaping of the vector into equivalently sized matrices
or tensors.

### Library configuration
```@contents
Pages = [
  "LibCEED.md",
]
```

### Ceed objects
```@contents
Pages = [
   "Ceed.md",
   "CeedVector.md",
   "ElemRestriction.md",
   "Basis.md",
   "QFunction.md",
   "Operator.md",
]
```

### Utilities
```@contents
Pages = [
   "Misc.md",
   "Globals.md",
   "Quadrature.md",
]
```

### C interface
```@contents
Pages = ["C.md"]
```
# Context

```@docs
Context
```
# Library configuration

By default, LibCEED.jl uses a "basic version" of the libCEED library that is
bundled as a pre-built binary. In order to access more advanced features (CUDA
support, architecture-specific compiler flags, etc.), users can use LibCEED.jl
with a other versions of the libCEED library (e.g. compiled from source).

This is achieved by by calling [`set_libceed_path!`](@ref) with the path to the
library file. The choice of library file is stored as a per-environment
preference. For changes to take effect, the Julia session must be restarted. The
library currently being used by LibCEED.jl can be queried using
[`get_libceed_path`](@ref).

The version number of the currently loaded libCEED library can also be queried
using [`ceedversion`](@ref).

```@docs
ceedversion
isrelease
set_libceed_path!
use_prebuilt_libceed!
get_libceed_path
```
# ElemRestriction

```@docs
ElemRestriction
ElemRestrictionNone
create_elem_restriction
create_elem_restriction_strided
apply!(r::ElemRestriction, u::CeedVector, ru::CeedVector; tmode=NOTRANSPOSE, request=RequestImmediate())
apply(r::ElemRestriction, u::AbstractVector; tmode=NOTRANSPOSE)
create_evector
create_lvector
create_vectors
getcompstride
getnumelements
getelementsize
getlvectorsize
getnumcomponents(r::ElemRestriction)
getmultiplicity!
getmultiplicity
```
# Examples

LibCEED.jl includes three short examples, which are analogues of the two
examples in `libCEED/examples/ceed`.

These examples are:
- `ex1-volume-c.jl`, an almost one-to-one translation of `ex1-volume.c`, using
  the low-level C interface. This example uses low-level user Q-functions
  defined in `ex1-function-c.jl`.
- `ex1-volume.jl`, a higher-level more idiomatic version of `ex1-volume.c`,
  using user Q-functions defined using [`@interior_qf`](@ref).
- `ex2-surface.jl`, a higher-level, idiomatic version of `ex2-surface.c`.
# Basis

!!! info "Column-major vs. row-major storage"
    libCEED internally uses row-major (C convention) storage of matrices,
    while Julia uses column-major (Fortran convention) storage.

    LibCEED.jl will typically handle the conversion between these formats by
    transposing or permuting the dimensions of the input and output matrices
    and tensors.

```@docs
Basis
BasisCollocated
create_tensor_h1_lagrange_basis
create_tensor_h1_basis
create_h1_basis
apply!(b::Basis, nelem, tmode::TransposeMode, emode::EvalMode, u::LibCEED.AbstractCeedVector, v::LibCEED.AbstractCeedVector)
apply(b::Basis, u::AbstractVector; nelem=1, tmode=NOTRANSPOSE, emode=EVAL_INTERP)
getdimension
gettopology
getnumcomponents(b::Basis)
getnumnodes
getnumnodes1d
getnumqpts
getnumqpts1d
getqref
getqweights
getinterp
getinterp1d
getgrad
getgrad1d
```
# Constants and Enumerations

```@docs
CeedScalar
CeedInt
QuadMode
MemType
CopyMode
EvalMode
TransposeMode
NormType
Topology
STRIDES_BACKEND
```
# CeedVector

```@docs
CeedVector
setvalue!
Base.setindex!(v::CeedVector, v2::CeedScalar)
Base.setindex!(v::CeedVector, v2::AbstractArray)
Base.Vector(v::CeedVector)
LinearAlgebra.norm(v::CeedVector, n::NormType)
LinearAlgebra.norm(v::CeedVector, p::Real)
@witharray
@witharray_read
witharray
witharray_read
setarray!
syncarray!
takearray!
scale!
LinearAlgebra.axpy!(a::Real, x::CeedVector, y::CeedVector)
pointwisemult!
```
# Defining User Q-Functions

An important feature of LibCEED.jl is the ability to define [user
Q-functions](https://libceed.readthedocs.io/en/latest/libCEEDapi/#gallery-of-qfunctions)
natively in Julia. These user Q-functions work with both the CPU and CUDA
backends.

User Q-functions describe the action of the $D$ operator at quadrature points
(see [libCEED's theoretical
framework](https://libceed.readthedocs.io/en/latest/libCEEDapi/#theoretical-framework)).
Since the Q-functions are invoked at every quadrature point, efficiency is
very important.

## Apply mass Q-function in C

Before describing how to define user Q-functions in Julia, we will briefly given
an example of a user Q-function defined in C. This is the "apply mass"
Q-function from `ex1-volume.c`, which computes the action of the mass operator.
The mass operator on each element can be written as $B^\intercal D B$, where $B$
is the basis operator, and $D$ represents multiplication by quadrature weights
and geometric factors (i.e. the determinant of the mesh transformation Jacobian
at each qudarture point). It is the action of $D$ that the Q-function must
implement. The C source of the Q-function is:

```c
/// libCEED Q-function for applying a mass operator
CEED_QFUNCTION(f_apply_mass)(void *ctx, const CeedInt Q,
                             const CeedScalar *const *in,
                             CeedScalar *const *out) {
  const CeedScalar *u = in[0], *qdata = in[1];
  CeedScalar *v = out[0];
  // Quadrature Point Loop
  CeedPragmaSIMD
  for (CeedInt i=0; i<Q; i++) {
    v[i] = qdata[i] * u[i];
  } // End of Quadrature Point Loop
  return 0;
}
```

From this example, we see that a user Q-function is a C callback that takes a
"data context" pointer, a number of quadrature points, and two arrays of arrays,
one for inputs, and one for outputs.

In this example, the first input array is `u`, which is the value of the trial
function evaluated at each quadrature point. The second input array is `qdata`,
which contains the precomputed geometric factors. There is only one output
array, `v`, which will store the pointwise product of `u` and `data`. Given the
definition of this Q-function, the `CeedQFunction` object is created by
```c
CeedQFunctionCreateInterior(ceed, 1, f_apply_mass, f_apply_mass_loc, &apply_qfunc);
CeedQFunctionAddInput(apply_qfunc, "u", 1, CEED_EVAL_INTERP);
CeedQFunctionAddInput(apply_qfunc, "qdata", 1, CEED_EVAL_NONE);
CeedQFunctionAddOutput(apply_qfunc, "v", 1, CEED_EVAL_INTERP);
```
When adding the inputs and outputs, `CEED_EVAL_INTERP` indicates that the $B$
basis operator should be used to interpolate the trial and test functions from
nodal points to quadrature points, and `CEED_EVAL_NONE` indicates that the
`qdata` is already precomputed at quadrature points, and no interpolation is
requried.

## Apply mass Q-function in Julia

We now replicate this Q-function in Julia. The main way of defining user
Q-functions in Julia is using the [`@interior_qf`](@ref) macro. The above C code
(both the definition of the Q-function, its creation, and adding the inputs and
outputs) is analogous to the following Julia code:

```julia
@interior_qf apply_qfunc = (
    ceed, Q,
    (u, :in, EVAL_INTERP, Q), (qdata, :in, EVAL_NONE, Q),
    (v, :out, EVAL_INTERP, Q),
    @inbounds @simd for i=1:Q
        v[i] = qdata[i]*u[i]
    end
)
```

This creates a [`QFunction`](@ref) object named `apply_qfunc`. The Q-function is
defined by the tuple on the right-hand side. `ceed` is the name of the
[`Ceed`](@ref) object where the Q-function will be created, and the second
argument, `Q`, is the name of that variable that will contain the number of
quadrature points. The next three arguments are specifications of the input and
output fields:
```julia
    (u, :in, EVAL_INTERP, Q),
    (qdata, :in, EVAL_NONE, Q),
    (v, :out, EVAL_INTERP, Q),
```
Each input or output field specification is a tuple, where the first entry is
the name of the array, and the second entry is either `:in` or `:out`, according
to whether the array is an input or output array. The third entry is the
[`EvalMode`](@ref) of the field. The remaining entries are the dimensions of the
array. The first dimension is always equal to the number of quadrature points.
In this case, all the arrays are simply vectors whose size is equal to the
number of quadrature points, but in more sophisticated examples (e.g. the [apply
diffusion Q-function](@ref applydiff)) these arrays could consists of vectors or
matrices at each quadrature point. After providing all of the array
specifications, the body of the Q-function is provided.

## [Apply diffusion Q-function in Julia](@id applydiff)

For a more sophisticated example of a Q-function, we consider the "apply
diffusion" Q-function, used in `ex2-surface`. This Q-function computes the
action of the diffusion operator. When written in the form $B^\intercal D B$, in
this case $B$ represents the basis gradient matrix, and $D$ represents
multiplication by $w \det(J) J^{-\intercal} J^{-1}$, where $J$ is the mesh
transformation Jacobian, and $w$ is the quadrature weight.

This Q-function is implemented in Julia as follows:
```julia
@interior_qf apply_qfunc = (
    ceed, Q, dim=dim,
    (du, :in, EVAL_GRAD, Q, dim),
    (qdata, :in, EVAL_NONE, Q, dim*(dim+1)÷2),
    (dv, :out, EVAL_GRAD, Q, dim),
    @inbounds @simd for i=1:Q
        dXdxdXdxT = getvoigt(@view(qdata[i,:]), CeedDim(dim))
        dui = SVector{dim}(@view(du[i,:]))
        dv[i,:] .= dXdxdXdxT*dui
    end
)
```
In contrast to the previous example, before the field specifications, this
Q-function includes a _constant definition_ `dim=dim`. The
[`@interior_qf`](@ref) macro allows for any number of constant definitions,
which make the specified values available within the body of the Q-function as
compile-time constants.

In this example, `dim` is either 1, 2, or 3 according to the spatial dimension
of the problem. When the user Q-function is defined, LibCEED.jl will JIT compile
the body of the Q-function and make it available to libCEED as a C callback. In
the body of this Q-function, `dim` will be available, and its value will be a
compile-time constant, allowing for (static) dispatch based on the value of
`dim`, and eliminating branching.

Note that `dim` is also available for use in the field specifications. In this
example, the field specifications are slightly more involved that in the
previous example. The arrays are given by
```julia
    (du, :in, EVAL_GRAD, Q, dim),
    (qdata, :in, EVAL_NONE, Q, dim*(dim+1)÷2),
    (dv, :out, EVAL_GRAD, Q, dim),
```
Note that the input array `du` has [`EvalMode`](@ref) `EVAL_GRAD`, meaning that
this array stores the gradient of the trial function at each quadrature point.
Therefore, at each quadrature point, `du` stores a vector of length `dim`, and
so the shape of `du` is `(Q, dim)`. Similarly, the action of $D$ is given by
$w \det(J) J^{-\intercal} J^{-1} \nabla u$, which is also a vector of length `dim` at
each quadrature point. This means that the output array `dv` also has shape
`(Q, dim)`.

The geometric factors stored in `qdata` represent the symmetric matrix $w
\det(J) J^{-\intercal} J^{-1}$ evaluated at every quadrature point. In order to
reduce data usage, instead of storing this data as a $d \times d$ matrix, we use
the fact that we know it is symmetric to only store $d(d+1)/2$ entries, and the
remaining entries we infer by symmetry. These entries are stored using the
[Voigt convention](https://en.wikipedia.org/wiki/Voigt_notation). LibCEED.jl
provides some [utilities](Misc.md#LibCEED.getvoigt) for storing and extracting
symmetric matrices stored in this fashion.

After the field specifications, we have the body of the Q-function:
```julia
@inbounds @simd for i=1:Q
    dXdxdXdxT = getvoigt(@view(qdata[i,:]), CeedDim(dim))
    dui = SVector{dim}(@view(du[i,:]))
    dv[i,:] .= dXdxdXdxT*dui
end
```
First, the matrix $w \det(J) J^{-\intercal} J^{-1}$ is stored in the variable
`dXdxdXdxT`. The symmetric entries of this matrix are accesed using
`@view(qdata[i,:])`, which avoids allocations. [`getvoigt`](@ref) is used to
convert from Voigt notation to a symmetric matrix, which returns a statically
sized `SMatrix`. The version for the correct spatial dimension is selected using
`CeedDim(dim)`, which allows for compile-time dispatch, since `dim` is a
constant whose value is known as a constant when the Q-function is JIT compiled.

Then, the gradient of $u$ at the given quadrature point is loaded as a
fixed-size `SVector`. The result is placed into the output array, where the
[`StaticArrays.jl`](https://github.com/JuliaArrays/StaticArrays.jl) package
evaluates `dXdxdXdxT*dui` using an optimized matrix-vector product for small
matrices (since their sizes are known statically).

## GPU Kernels

If the `Ceed` resource uses a CUDA backend, then the user Q-functions defined
using [`@interior_qf`](@ref) are automatically compiled as CUDA kernels using
[`CUDA.jl`](https://github.com/JuliaGPU/CUDA.jl). Some Julia features are not
available in GPU code (for example, dynamic dispatch), so if the Q-function is
intended to be run on the GPU, the user should take care when defining the body
of the user Q-function.
# Low-level C interface

The low-level interface (provided in the `LibCEED.C` module) is in one-to-one
correspondence with the C libCEED iterface, and is automatically generated (with
some minor manual modifications) using the Julia package
[Clang.jl](https://github.com/JuliaInterop/Clang.jl/). The code used to generate
bindings is available in `generate_bindings.jl`.

With the low-level interface, the user is responsible for freeing all allocated
memory (calling the appropriate `Ceed*Destroy` functions). This interface is not
type-safe, and calling functions with the wrong arguments can cause libCEED to
crash.

It is generally recommended for users to use the Julia interface exported from
the `LibCEED` module, unless other specific low-level functionality is required.
# Operator

```@docs
Operator
create_composite_operator
apply!(op::Operator, vin::LibCEED.AbstractCeedVector, vout::LibCEED.AbstractCeedVector; request)
apply_add!
assemble_diagonal!
assemble_add_diagonal!
```
# Ceed

```@docs
Ceed
getresource
get_preferred_memtype
isdeterministic
iscuda
```
# Quadrature

```@docs
gauss_quadrature
lobatto_quadrature
```
# QFunction

```@docs
QFunction
@interior_qf
create_interior_qfunction(::Ceed, ::AbstractString)
create_identity_qfunction
set_context!
apply!(qf::QFunction, Q, vin, vout)
```
# libCEED: Tests

This page provides a brief description of the tests for the libCEED library.

The tests are organized by API object, and some tests are further organized,
as required.

0. Ceed Tests
1. CeedVector Tests  
    1.0. CeedVector general tests  
    1.1. CeedVector error tests
2. CeedElemRestriction Tests
3. CeedBasis Tests  
    3.0. CeedBasis utility tests  
    3.1. CeedBasis tensor basis tests  
    3.2. CeedBasis simplex basis tests
4. CeedQFunction Tests  
    4.0. CeedQFunction user code tests  
    4.1. CeedQFunction gallery code tests
5. CeedOperator Tests  
    5.0. CeedOperator with tensor bases tests  
    5.1. CeedOperator with simplex bases tests  
    5.2. CeedOperator with operator composition tests  
    5.3. CeedOperator and CeedQFunction assembly tests  
    5.4. CeedOperator element inverse tests  
    5.5. CeedOperator multigrid level tests
# libCEED: Benchmarks

This directory contains benchmark problems for performance evaluation of libCEED
backends.

## Running the Benchmarks

Example:
```sh
benchmark.sh -c /cpu/self -r petsc-bpsraw.sh -b bp1 -n 16 -p 16
```
where the option `-c <specs-list>` specifies a list of libCEED specs to
benchmark, `-b <bp-list>` specifies a list of CEED benchmark problems to run,
`-n 16` is the total number of processors and `-p 16` is the number
of processors per node.

Multiple backends, benchmark problems, and processor configurations can be
benchmarked with:
```sh
benchmark.sh -c "/cpu/self/ref/serial /cpu/self/ref/blocked" -r petsc-bpsraw.sh -b "bp1 bp3" -n "16 32 64" -p "16 32 64"
```

The results from the benchmarks are written to files named `*-output.txt`.

For a short help message, use the option `-h`.

When running the tests `petsc-bpsraw.sh`, the following
variables can be set on the command line:
* `max_dofs_node=<number>`, e.g. `max_dofs_node=1000000` - this sets the upper
  bound of the problem sizes, per compute node; the default value is 3*2^20.
* `max_p=<number>`, e.g. `max_p=12` - this sets the highest degree for which the
  tests will be run (the lowest degree is 1); the default value is 8.

## Post-processing the results

After generating the results, use the `postprocess-plot.py` script (which
requires the python package matplotlib) or the `postprocess-table.py` script,
e.g.:
```sh
python postprocess-plot.py petsc-bpsraw-bp1-*-output.txt
```
The plot ranges and some other options can be adjusted by editing the values
in the beginning of the script `postprocess-plot.py`.

Note that the `postprocess-*.py` scripts can read multiple files at a time just
by listing them on the command line and also read the standard input if no files
were specified on the command line.
# libCEED for Python

## Install

To install libCEED for Python, run

    pip install libceed

or in a clone of the repository via `pip install .`

## Examples and Tutorials

For examples and short tutorials see the folder `examples/tutorials`. It
contains some [Jupyter](https://jupyter.org/) notebooks using Python and C.
Jupyter can be installed locally so that users can edit and interact with these
notebook.

`tutorial-0`-`tutorial-5` illustrate libCEED for Python, each one focusing on one
class of objects.

`tutorial-6` shows a standalone libCEED C example.
# libCEED Python Tests

These files provide libCEED for Python tests. Full examples of finite element
operators can be found in the file `test-5-operator.py`.

## Testing

To run the tests, first build the user QFunctions file by running

  python setup-qfunctions.py build

Then, to run the test suite use the command

  pytest test-*.py --ceed /cpu/self/ref/serial

## Building QFunctions

To build user defined QFunctions, modify `libceed-qfunctions.c` to include
the apropriate QFunction single source file and run
`python setup-qfunctions.py build`. The files `test-4-qfunction.py` and
`test-5-operator.py` both contain the example function `load_qfs_so()` for
loading the user defined QFunctions so the QFunction pointers can be passed to
libCEED.
(examples)=

# Examples

This section contains a mathematical description of all examples provided with libCEED
in the {file}`examples/` directory.
These examples are meant to demonstrate use of libCEED from standalone definition of operators to integration with external libraries such as [PETSc](https://www.mcs.anl.gov/petsc), [MFEM](https://mfem.org), and [Nek5000](https://nek5000.mcs.anl.gov/), as well as more substantial mini-apps.

```{toctree}
:maxdepth: 2

notation.md
ceed/index.md
bps.md
petsc/index.md
fluids/index.md
solids/index.md
```
(bps)=

# CEED Bakeoff Problems

```{include} ./README.md
:start-after: bps-inclusion-marker
:end-before: bps-exclusion-marker
```

(mass-operator)=

## Mass Operator

The Mass Operator used in BP1 and BP2 is defined via the $L^2$ projection
problem, posed as a weak form on a Hilbert space $V^p \subset H^1$, i.e.,
find $u \in V^p$ such that for all $v \in V^p$

$$
\langle v,u \rangle = \langle v,f \rangle ,
$$ (eq-general-weak-form)

where $\langle v,u\rangle$ and $\langle v,f\rangle$ express the continuous
bilinear and linear forms, respectively, defined on $V^p$, and, for sufficiently
regular $u$, $v$, and $f$, we have:

$$
\begin{aligned} \langle v,u \rangle &:= \int_{\Omega} \, v \, u \, dV ,\\ \langle v,f \rangle &:= \int_{\Omega} \, v \, f \, dV . \end{aligned}
$$

Following the standard finite/spectral element approach, we formally
expand all functions in terms of basis functions, such as

$$
\begin{aligned}
u(\bm x) &= \sum_{j=1}^n u_j \, \phi_j(\bm x) ,\\
v(\bm x) &= \sum_{i=1}^n v_i \, \phi_i(\bm x) .
\end{aligned}
$$ (eq-nodal-values)

The coefficients $\{u_j\}$ and $\{v_i\}$ are the nodal values of $u$
and $v$, respectively. Inserting the expressions {eq}`eq-nodal-values`
into {eq}`eq-general-weak-form`, we obtain the inner-products

$$
\langle v,u \rangle = \bm v^T M \bm u , \qquad  \langle v,f\rangle =  \bm v^T \bm b \,.
$$ (eq-inner-prods)

Here, we have introduced the mass matrix, $M$, and the right-hand side,
$\bm b$,

$$
M_{ij} :=  (\phi_i,\phi_j), \;\; \qquad b_{i} :=  \langle \phi_i, f \rangle,
$$

each defined for index sets $i,j \; \in \; \{1,\dots,n\}$.

(laplace-operator)=

## Laplace's Operator

The Laplace's operator used in BP3-BP6 is defined via the following variational
formulation, i.e., find $u \in V^p$ such that for all $v \in V^p$

$$
a(v,u) = \langle v,f \rangle , \,
$$

where now $a (v,u)$ expresses the continuous bilinear form defined on
$V^p$ for sufficiently regular $u$, $v$, and $f$, that is:

$$
\begin{aligned} a(v,u) &:= \int_{\Omega}\nabla v \, \cdot \, \nabla u \, dV ,\\ \langle v,f \rangle &:= \int_{\Omega} \, v \, f \, dV . \end{aligned}
$$

After substituting the same formulations provided in {eq}`eq-nodal-values`,
we obtain

$$
a(v,u) = \bm v^T K \bm u ,
$$

in which we have introduced the stiffness (diffusion) matrix, $K$, defined as

$$
K_{ij} = a(\phi_i,\phi_j),
$$

for index sets $i,j \; \in \; \{1,\dots,n\}$.
# libCEED: Examples

This page provides a brief description of the examples for the libCEED
library.

## Basic libCEED Examples

Two examples that rely only upon libCEED without any external libraries are provided in the [ceed/](./ceed) folder. For more details, please see the dedicated [documentation section](https://libceed.readthedocs.io/en/latest/examples/ceed/index.html).

## Bakeoff Problems

% bps-inclusion-marker

The Center for Efficient Exascale Discretizations (CEED) uses Bakeoff Problems (BPs)
to test and compare the performance of high-order finite element implementations. The
definitions of the problems are given on the ceed
[website](https://ceed.exascaleproject.org/bps/). Each of the following bakeoff
problems that use external discretization libraries (such as MFEM, PETSc, and Nek5000)
are located in the subdirectories `mfem/`, `petsc/`, and
`nek5000/`, respectively.

Here we provide a short summary:

:::{list-table}
:header-rows: 1
:widths: auto
* - User code
  - Supported BPs
* - `mfem`
  - * BP1 (scalar mass operator) with $Q=P+1$
    * BP3 (scalar Laplace operator) with $Q=P+1$
* - `petsc`
  - * BP1 (scalar mass operator) with $Q=P+1$
    * BP2 (vector mass operator) with $Q=P+1$
    * BP3 (scalar Laplace operator) with $Q=P+1$
    * BP4 (vector Laplace operator) with $Q=P+1$
    * BP5 (collocated scalar Laplace operator) with $Q=P$
    * BP6 (collocated vector Laplace operator) with $Q=P$
* - `nek5000`
  - * BP1 (scalar mass operator) with $Q=P+1$
    * BP3 (scalar Laplace operator) with $Q=P+1$
:::

These are all **T-vector**-to-**T-vector** and include parallel scatter, element
scatter, element evaluation kernel, element gather, and parallel gather (with the
parallel gathers/scatters done externally to libCEED).

BP1 and BP2 are $L^2$ projections, and thus have no boundary condition.
The rest of the BPs have homogeneous Dirichlet boundary conditions.

The BPs are parametrized by the number $P$ of Gauss-Legendre-Lobatto nodal points
(with $P=p+1$, and $p$ the degree of the basis polynomial) for the Lagrange
polynomials, as well as the number of quadrature points, $Q$.
A $Q$-point Gauss-Legendre quadrature is used for all BPs except BP5 and BP6,
which choose $Q = P$ and Gauss-Legendre-Lobatto quadrature to collocate with the
interpolation nodes. This latter choice is popular in applications that use spectral
element methods because it produces a diagonal mass matrix (enabling easy explicit
time integration) and significantly reduces the number of floating point operations
to apply the operator.

% bps-exclusion-marker

For a more detailed description of the operators employed in the BPs, please see the dedicated [BPs documentation section](https://libceed.readthedocs.io/en/latest/examples/bps.html).

## PETSc+libCEED Navier-Stokes Solver

The Navier-Stokes problem solves the compressible Navier-Stokes
equations using an explicit or implicit time integration. A more detailed
description of the problem formulation can be found in the
[fluids/](./fluids) folder and the corresponding [fluids documentation page](https://libceed.readthedocs.io/en/latest/examples/fluids/index.html).

## PETSc+libCEED Solid mechanics elasticity mini-app

This example solves the steady-state static momentum balance equations using unstructured high-order finite/spectral element spatial discretizations. A more detailed
description of the problem formulation can be found in the
[solids/](./solids) folder and the corresponding [solids documentation page](https://libceed.readthedocs.io/en/latest/examples/solids/index.html).

## PETSc+libCEED Surface Area Examples

These examples, located in the [petsc/](./petsc) folder, use the mass operator to compute the surface area of a
cube or a discrete cubed-sphere, using PETSc. For a detailed description, please see the corresponding [area documentation page](https://libceed.readthedocs.io/en/latest/examples/petsc/index.html#area).

## PETSc+libCEED Bakeoff Problems on the Cubed-Sphere

These examples, located in the [petsc/](./petsc) folder, reproduce the Bakeoff Problems 1-6 on a discrete
cubed-sphere, using PETSc. For a detailed description, please see the corresponding [problems on the cubed-sphere documentation page](https://libceed.readthedocs.io/en/latest/examples/petsc/index.html#bakeoff-problems-on-the-cubed-sphere).

## Running Examples

To build the examples, set the `MFEM_DIR`, `PETSC_DIR`, and
`NEK5K_DIR` variables and, from the `examples/` directory, run

```{include} ../README.md
:start-after: running-examples-inclusion-marker
:end-before: benchmarks-marker
```
(common-notation)=

# Common notation

For most of our examples, the spatial discretization
uses high-order finite elements/spectral elements, namely, the high-order Lagrange
polynomials defined over $P$ non-uniformly spaced nodes, the
Gauss-Legendre-Lobatto (GLL) points, and quadrature points $\{q_i\}_{i=1}^Q$, with
corresponding weights $\{w_i\}_{i=1}^Q$ (typically the ones given by Gauss
or Gauss-Lobatto quadratures, that are built in the library).

We discretize the domain, $\Omega \subset \mathbb{R}^d$ (with $d=1,2,3$,
typically) by letting $\Omega = \bigcup_{e=1}^{N_e}\Omega_e$, with $N_e$
disjoint elements. For most examples we use unstructured meshes for which the elements
are hexahedra (although this is not a requirement in libCEED).

The physical coordinates are denoted by
$\bm{x}=(x,y,z) \equiv (x_0,x_1,x_2) \in\Omega_e$,
while the reference coordinates are represented as
$\bm{X}=(X,Y,Z) \equiv (X_0,X_1,X_2) \in \textrm{I}=[-1,1]^3$
(for $d=3$).
# Standalone libCEED

The following two examples have no dependencies, and are designed to be self-contained.
For additional examples that use external discretization libraries (MFEM, PETSc, Nek5000
etc.) see the subdirectories in {file}`examples/`.

(ex1-volume)=

## Ex1-Volume

This example is located in the subdirectory {file}`examples/ceed`. It illustrates a
simple usage of libCEED to compute the volume of a given body using a matrix-free
application of the mass operator. Arbitrary mesh and solution orders in 1D, 2D, and 3D
are supported from the same code.

This example shows how to compute line/surface/volume integrals of a 1D, 2D, or 3D
domain $\Omega$ respectively, by applying the mass operator to a vector of
$1$s. It computes:

$$
I = \int_{\Omega} 1 \, dV .
$$ (eq-ex1-volume)

Using the same notation as in {ref}`theoretical-framework`, we write here the vector
$u(x)\equiv 1$ in the Galerkin approximation,
and find the volume of $\Omega$ as

$$
\sum_e \int_{\Omega_e} v(x) 1 \, dV
$$ (volume-sum)

with $v(x) \in \mathcal{V}_p = \{ v \in H^{1}(\Omega_e) \,|\, v \in P_p(\bm{I}), e=1,\ldots,N_e \}$,
the test functions.

(ex2-surface)=

## Ex2-Surface

This example is located in the subdirectory {file}`examples/ceed`. It computes the
surface area of a given body using matrix-free application of a diffusion operator.
Similar to {ref}`Ex1-Volume`, arbitrary mesh and solution orders in 1D, 2D, and 3D
are supported from the same code. It computes:

$$
I = \int_{\partial \Omega} 1 \, dS ,
$$ (eq-ex2-surface)

by applying the divergence theorem.
In particular, we select $u(\bm x) = x_0 + x_1 + x_2$, for which $\nabla u = [1, 1, 1]^T$, and thus $\nabla u \cdot \hat{\bm n} = 1$.

Given Laplace's equation,

$$
\nabla \cdot \nabla u = 0, \textrm{ for  } \bm{x} \in \Omega ,
$$

let us multiply by a test function $v$ and integrate by parts to obtain

$$
\int_\Omega \nabla v \cdot \nabla u \, dV - \int_{\partial \Omega} v \nabla u \cdot \hat{\bm n}\, dS = 0 .
$$

Since we have chosen $u$ such that $\nabla u \cdot \hat{\bm n} = 1$, the boundary integrand is $v 1 \equiv v$. Hence, similar to {eq}`volume-sum`, we can evaluate the surface integral by applying the volumetric Laplacian as follows

$$
\int_\Omega \nabla v \cdot \nabla u \, dV \approx \sum_e \int_{\partial \Omega_e} v(x) 1 \, dS .
$$
## libCEED: Basic Examples

Two examples are provided that rely only upon libCEED without any external
libraries.

### Example 1: ex1-volume

This example uses the mass matrix to compute the length, area, or volume of a
region, depending upon runtime parameters.

### Example 2: ex2-surface

This example uses the diffusion matrix to compute the surface area of a region,
in 1D, 2D or 3D, depending upon runtime parameters.
(example-petsc-navier-stokes)=

# Compressible Navier-Stokes mini-app

This example is located in the subdirectory {file}`examples/fluids`.
It solves the time-dependent Navier-Stokes equations of compressible gas dynamics in a static Eulerian three-dimensional frame using unstructured high-order finite/spectral element spatial discretizations and explicit or implicit high-order time-stepping (available in PETSc).
Moreover, the Navier-Stokes example has been developed using PETSc, so that the pointwise physics (defined at quadrature points) is separated from the parallelization and meshing concerns.

## Running the mini-app

```{include} README.md
:start-after: inclusion-fluids-marker
```
## The Navier-Stokes equations

The mathematical formulation (from {cite}`giraldoetal2010`, cf. SE3) is given in what follows.
The compressible Navier-Stokes equations in conservative form are

$$
\begin{aligned}
\frac{\partial \rho}{\partial t} + \nabla \cdot \bm{U} &= 0 \\
\frac{\partial \bm{U}}{\partial t} + \nabla \cdot \left( \frac{\bm{U} \otimes \bm{U}}{\rho} + P \bm{I}_3 -\bm\sigma \right) + \rho g \bm{\hat k} &= 0 \\
\frac{\partial E}{\partial t} + \nabla \cdot \left( \frac{(E + P)\bm{U}}{\rho} -\bm{u} \cdot \bm{\sigma} - k \nabla T \right) &= 0 \, , \\
\end{aligned}
$$ (eq-ns)

where $\bm{\sigma} = \mu(\nabla \bm{u} + (\nabla \bm{u})^T + \lambda (\nabla \cdot \bm{u})\bm{I}_3)$ is the Cauchy (symmetric) stress tensor, with $\mu$ the dynamic viscosity coefficient, and $\lambda = - 2/3$ the Stokes hypothesis constant.
In equations {eq}`eq-ns`, $\rho$ represents the volume mass density, $U$ the momentum density (defined as $\bm{U}=\rho \bm{u}$, where $\bm{u}$ is the vector velocity field), $E$ the total energy density (defined as $E = \rho e$, where $e$ is the total energy), $\bm{I}_3$ represents the $3 \times 3$ identity matrix, $g$ the gravitational acceleration constant, $\bm{\hat{k}}$ the unit vector in the $z$ direction, $k$ the thermal conductivity constant, $T$ represents the temperature, and $P$ the pressure, given by the following equation of state

$$
P = \left( {c_p}/{c_v} -1\right) \left( E - {\bm{U}\cdot\bm{U}}/{(2 \rho)} - \rho g z \right) \, ,
$$ (eq-state)

where $c_p$ is the specific heat at constant pressure and $c_v$ is the specific heat at constant volume (that define $\gamma = c_p / c_v$, the specific heat ratio).

The system {eq}`eq-ns` can be rewritten in vector form

$$
\frac{\partial \bm{q}}{\partial t} + \nabla \cdot \bm{F}(\bm{q}) -S(\bm{q}) = 0 \, ,
$$ (eq-vector-ns)

for the state variables 5-dimensional vector

$$
\bm{q} =        \begin{pmatrix}            \rho \\            \bm{U} \equiv \rho \bm{ u }\\            E \equiv \rho e        \end{pmatrix}        \begin{array}{l}            \leftarrow\textrm{ volume mass density}\\            \leftarrow\textrm{ momentum density}\\            \leftarrow\textrm{ energy density}        \end{array}
$$

where the flux and the source terms, respectively, are given by

$$
\begin{aligned}
\bm{F}(\bm{q}) &=
\underbrace{\begin{pmatrix}
    \bm{U}\\
    {(\bm{U} \otimes \bm{U})}/{\rho} + P \bm{I}_3 \\
    {(E + P)\bm{U}}/{\rho}
\end{pmatrix}}_{\bm F_{\text{adv}}} +
\underbrace{\begin{pmatrix}
0 \\
-  \bm{\sigma} \\
 - \bm{u}  \cdot \bm{\sigma} - k \nabla T
\end{pmatrix}}_{\bm F_{\text{diff}}},\\
S(\bm{q}) &=
- \begin{pmatrix}
    0\\
    \rho g \bm{\hat{k}}\\
    0
\end{pmatrix}.
\end{aligned}
$$ (eq-ns-flux)

Let the discrete solution be

$$
\bm{q}_N (\bm{x},t)^{(e)} = \sum_{k=1}^{P}\psi_k (\bm{x})\bm{q}_k^{(e)}
$$

with $P=p+1$ the number of nodes in the element $e$.
We use tensor-product bases $\psi_{kji} = h_i(X_0)h_j(X_1)h_k(X_2)$.

For the time discretization, we use two types of time stepping schemes.

- Explicit time-stepping method

  The following explicit formulation is solved with the adaptive Runge-Kutta-Fehlberg (RKF4-5) method by default (any explicit time-stepping scheme available in PETSc can be chosen at runtime)

  $$
  \bm{q}_N^{n+1} = \bm{q}_N^n + \Delta t \sum_{i=1}^{s} b_i k_i \, ,
  $$

  where

  $$
  \begin{aligned}
     k_1 &= f(t^n, \bm{q}_N^n)\\
     k_2 &= f(t^n + c_2 \Delta t, \bm{q}_N^n + \Delta t (a_{21} k_1))\\
     k_3 &= f(t^n + c_3 \Delta t, \bm{q}_N^n + \Delta t (a_{31} k_1 + a_{32} k_2))\\
     \vdots&\\
     k_i &= f\left(t^n + c_i \Delta t, \bm{q}_N^n + \Delta t \sum_{j=1}^s a_{ij} k_j \right)\\
  \end{aligned}
  $$

  and with

  $$
  f(t^n, \bm{q}_N^n) = - [\nabla \cdot \bm{F}(\bm{q}_N)]^n + [S(\bm{q}_N)]^n \, .
  $$

- Implicit time-stepping method

  This time stepping method which can be selected using the option `-implicit` is solved with Backward Differentiation Formula (BDF) method by default (similarly, any implicit time-stepping scheme available in PETSc can be chosen at runtime).
  The implicit formulation solves nonlinear systems for $\bm q_N$:

  $$
  \bm f(\bm q_N) \equiv \bm g(t^{n+1}, \bm{q}_N, \bm{\dot{q}}_N) = 0 \, ,
  $$ (eq-ts-implicit-ns)

  where the time derivative $\bm{\dot q}_N$ is defined by

  $$
  \bm{\dot{q}}_N(\bm q_N) = \alpha \bm q_N + \bm z_N
  $$

  in terms of $\bm z_N$ from prior state and $\alpha > 0$, both of which depend on the specific time integration scheme (backward difference formulas, generalized alpha, implicit Runge-Kutta, etc.).
  Each nonlinear system {eq}`eq-ts-implicit-ns` will correspond to a weak form, as explained below.
  In determining how difficult a given problem is to solve, we consider the Jacobian of {eq}`eq-ts-implicit-ns`,

  $$
  \frac{\partial \bm f}{\partial \bm q_N} = \frac{\partial \bm g}{\partial \bm q_N} + \alpha \frac{\partial \bm g}{\partial \bm{\dot q}_N}.
  $$

  The scalar "shift" $\alpha$ scales inversely with the time step $\Delta t$, so small time steps result in the Jacobian being dominated by the second term, which is a sort of "mass matrix", and typically well-conditioned independent of grid resolution with a simple preconditioner (such as Jacobi).
  In contrast, the first term dominates for large time steps, with a condition number that grows with the diameter of the domain and polynomial degree of the approximation space.
  Both terms are significant for time-accurate simulation and the setup costs of strong preconditioners must be balanced with the convergence rate of Krylov methods using weak preconditioners.

To obtain a finite element discretization, we first multiply the strong form {eq}`eq-vector-ns` by a test function $\bm v \in H^1(\Omega)$ and integrate,

$$
\int_{\Omega} \bm v \cdot \left(\frac{\partial \bm{q}_N}{\partial t} + \nabla \cdot \bm{F}(\bm{q}_N) - \bm{S}(\bm{q}_N) \right) \,dV = 0 \, , \; \forall \bm v \in \mathcal{V}_p\,,
$$

with $\mathcal{V}_p = \{ \bm v(\bm x) \in H^{1}(\Omega_e) \,|\, \bm v(\bm x_e(\bm X)) \in P_p(\bm{I}), e=1,\ldots,N_e \}$ a mapped space of polynomials containing at least polynomials of degree $p$ (with or without the higher mixed terms that appear in tensor product spaces).

Integrating by parts on the divergence term, we arrive at the weak form,

$$
\begin{aligned}
\int_{\Omega} \bm v \cdot \left( \frac{\partial \bm{q}_N}{\partial t} - \bm{S}(\bm{q}_N) \right)  \,dV
- \int_{\Omega} \nabla \bm v \!:\! \bm{F}(\bm{q}_N)\,dV & \\
+ \int_{\partial \Omega} \bm v \cdot \bm{F}(\bm q_N) \cdot \widehat{\bm{n}} \,dS
  &= 0 \, , \; \forall \bm v \in \mathcal{V}_p \,,
\end{aligned}
$$ (eq-weak-vector-ns)

where $\bm{F}(\bm q_N) \cdot \widehat{\bm{n}}$ is typically replaced with a boundary condition.

:::{note}
The notation $\nabla \bm v \!:\! \bm F$ represents contraction over both fields and spatial dimensions while a single dot represents contraction in just one, which should be clear from context, e.g., $\bm v \cdot \bm S$ contracts over fields while $\bm F \cdot \widehat{\bm n}$ contracts over spatial dimensions.
:::

We solve {eq}`eq-weak-vector-ns` using a Galerkin discretization (default) or a stabilized method, as is necessary for most real-world flows.

Galerkin methods produce oscillations for transport-dominated problems (any time the cell Péclet number is larger than 1), and those tend to blow up for nonlinear problems such as the Euler equations and (low-viscosity/poorly resolved) Navier-Stokes, in which case stabilization is necessary.
Our formulation follows {cite}`hughesetal2010`, which offers a comprehensive review of stabilization and shock-capturing methods for continuous finite element discretization of compressible flows.

- **SUPG** (streamline-upwind/Petrov-Galerkin)

  In this method, the weighted residual of the strong form {eq}`eq-vector-ns` is added to the Galerkin formulation {eq}`eq-weak-vector-ns`.
  The weak form for this method is given as

  $$
  \begin{aligned}
  \int_{\Omega} \bm v \cdot \left( \frac{\partial \bm{q}_N}{\partial t} - \bm{S}(\bm{q}_N) \right)  \,dV
  - \int_{\Omega} \nabla \bm v \!:\! \bm{F}(\bm{q}_N)\,dV & \\
  + \int_{\partial \Omega} \bm v \cdot \bm{F}(\bm{q}_N) \cdot \widehat{\bm{n}} \,dS & \\
  + \int_{\Omega} \bm{P}(\bm v)^T \, \left( \frac{\partial \bm{q}_N}{\partial t} \, + \,
  \nabla \cdot \bm{F} \, (\bm{q}_N) - \bm{S}(\bm{q}_N) \right) \,dV &= 0
  \, , \; \forall \bm v \in \mathcal{V}_p
  \end{aligned}
  $$ (eq-weak-vector-ns-supg)

  This stabilization technique can be selected using the option `-stab supg`.

- **SU** (streamline-upwind)

  This method is a simplified version of *SUPG* {eq}`eq-weak-vector-ns-supg` which is developed for debugging/comparison purposes. The weak form for this method is

  $$
  \begin{aligned}
  \int_{\Omega} \bm v \cdot \left( \frac{\partial \bm{q}_N}{\partial t} - \bm{S}(\bm{q}_N) \right)  \,dV
  - \int_{\Omega} \nabla \bm v \!:\! \bm{F}(\bm{q}_N)\,dV & \\
  + \int_{\partial \Omega} \bm v \cdot \bm{F}(\bm{q}_N) \cdot \widehat{\bm{n}} \,dS & \\
  + \int_{\Omega} \mathcal{P}(\bm v)^T \, \nabla \cdot \bm{F} \, (\bm{q}_N) \,dV
  & = 0 \, , \; \forall \bm v \in \mathcal{V}_p
  \end{aligned}
  $$ (eq-weak-vector-ns-su)

  This stabilization technique can be selected using the option `-stab su`.

In both {eq}`eq-weak-vector-ns-su` and {eq}`eq-weak-vector-ns-supg`, $\mathcal P$ is called the *perturbation to the test-function space*, since it modifies the original Galerkin method into *SUPG* or *SU* schemes.
It is defined as

$$
\mathcal P(\bm v) \equiv \left(\bm{\tau} \cdot \frac{\partial \bm{F}_{\text{adv}} (\bm{q}_N)}{\partial \bm{q}_N} \right)^T \, \nabla \bm v\,,
$$

where parameter $\bm{\tau} \in \mathbb R^{3\times 3}$ (spatial indices) or $\bm \tau \in \mathbb R^{5\times 5}$ (field indices) is an intrinsic time scale matrix.
This expression contains the flux Jacobian, which we express in variational notation by differentiating the advective flux $\bm F_{\text{adv}}$ of {eq}`eq-ns-flux`

$$
\begin{aligned}
\diff\bm F_{\text{adv}}(\diff\bm q; \bm q) &= \frac{\partial \bm F_{\text{adv}}}{\partial \bm q} \diff\bm q \\
&= \begin{pmatrix}
\diff\bm U \\
(\diff\bm U \otimes \bm U + \bm U \otimes \diff\bm U)/\rho - (\bm U \otimes \bm U)/\rho^2 \diff\rho + \diff P \bm I_3 \\
(E + P)\diff\bm U/\rho + (\diff E + \diff P)\bm U/\rho - (E + P) \bm U/\rho^2 \diff\rho
\end{pmatrix},
\end{aligned}
$$

where $\diff P$ is defined by differentiating {eq}`eq-state`.
In this notation, we may equivalently write the stabilization term as

$$
\mathcal P(\bm v)^T \bm r = \nabla \bm v \bm\tau \diff\bm F_{\text{adv}}(\bm r),
$$

where $\bm r$ is the strong form residual.
Note that both $\nabla \bm v$ and $\diff \bm F$ are $5\times 3$ matrices and that $\bm\tau$ can be defined with spatial indices, or field indices, leading to a stabilization term of $(\nabla \bm v)_{\alpha i} \tau_{ij} \diff \bm F_{\alpha j}$ for spatial or $(\nabla \bm v)_{\alpha i} \tau_{\alpha \beta} \diff \bm F_{\beta i}$ for field, where $\alpha,\beta$ are field indices and $i,j$ are spatial indices.

:::{dropdown} Stabilization scale $\bm\tau$
A velocity vector $\bm u$ can be pulled back to the reference element as $\bm u_{\bm X} = \nabla_{\bm x}\bm X \cdot \bm u$, with units of reference length (non-dimensional) per second.
To build intuition, consider a boundary layer element of dimension $(1, \epsilon)$, for which $\nabla_{\bm x} \bm X = \bigl(\begin{smallmatrix} 2 & \\ & 2/\epsilon \end{smallmatrix}\bigr)$.
So a small normal component of velocity will be amplified (by a factor of the aspect ratio $1/\epsilon$) in this transformation.
The ratio $\lVert \bm u \rVert / \lVert \bm u_{\bm X} \rVert$ is a covariant measure of (half) the element length in the direction of the velocity.
A contravariant measure of element length in the direction of a unit vector $\hat{\bm n}$ is given by $\lVert \bigl(\nabla_{\bm X} \bm x\bigr)^T \hat{\bm n} \rVert$.
While $\nabla_{\bm X} \bm x$ is readily computable, its inverse $\nabla_{\bm x} \bm X$ is needed directly in finite element methods and thus more convenient for our use.
If we consider a parallelogram, the covariant measure is larger than the contravariant measure for vectors pointing between acute corners and the opposite holds for vectors between oblique corners.

The cell Péclet number is classically defined by $\mathrm{Pe}_h = \lVert \bm u \rVert h / (2 \kappa)$ where $\kappa$ is the diffusivity (units of $m^2/s$).
This can be generalized to arbitrary grids by defining the local Péclet number

$$
\mathrm{Pe} = \frac{\lVert \bm u \rVert^2}{\lVert \bm u_{\bm X} \rVert \kappa}.
$$ (eq-peclet)

For scalar advection-diffusion, the stabilization is a scalar

$$
\tau = \frac{\xi(\mathrm{Pe})}{\lVert \bm u_{\bm X} \rVert},
$$ (eq-tau-advdiff)

where $\xi(\mathrm{Pe}) = \coth \mathrm{Pe} - 1/\mathrm{Pe}$ approaches 1 at large local Péclet number.
Note that $\tau$ has units of time and, in the transport-dominated limit, is proportional to element transit time in the direction of the propagating wave.
For advection-diffusion, $\bm F(q) = \bm u q$, and thus the perturbed test function is

$$
\mathcal P(v) = \tau \bm u \cdot \nabla v = \tau \bm u_{\bm X} \nabla_{\bm X} v.
$$ (eq-test-perturbation-advdiff)

See {cite}`hughesetal2010` equations 15-17 and 34-36 for further discussion of this formulation.

For the Navier-Stokes and Euler equations in primitive variables, {cite}`whiting2003hierarchical` defines a $5\times 5$ diagonal stabilization consisting of
1. continuity stabilization $\tau_c$
2. momentum stabilization $\tau_m$
3. energy stabilization $\tau_E$

However, since our equations are in conservative form, we follow {cite}`hughesetal2010` in defining a $3\times 3$ diagonal stabilization according to spatial criterion 2 (equation 27) as follows.

$$
\tau_{ii} = c_{\tau} \frac{2 \xi(\mathrm{Pe})}{(\lambda_{\max \text{abs}})_i \lVert \nabla_{x_i} \bm X \rVert}
$$ (eq-tau-conservative)

where $c_{\tau}$ is a multiplicative constant reported to be optimal at 0.5 for linear elements, $\hat{\bm n}_i$ is a unit vector in direction $i$, and $\nabla_{x_i} = \hat{\bm n}_i \cdot \nabla_{\bm x}$ is the derivative in direction $i$.
The flux Jacobian $\frac{\partial \bm F_{\text{adv}}}{\partial \bm q} \cdot \hat{\bm n}_i$ in each direction $i$ is a $5\times 5$ matrix with spectral radius $(\lambda_{\max \text{abs}})_i$ equal to the fastest wave speed.
The complete set of eigenvalues of the Euler flux Jacobian in direction $i$ are (e.g., {cite}`toro2009`)

$$
\Lambda_i = [u_i - a, u_i, u_i, u_i, u_i+a],
$$ (eq-eigval-advdiff)

where $u_i = \bm u \cdot \hat{\bm n}_i$ is the velocity component in direction $i$ and $a = \sqrt{\gamma P/\rho}$ is the sound speed for ideal gasses.
Note that the first and last eigenvalues represent nonlinear acoustic waves while the middle three are linearly degenerate, carrying a contact wave (temperature) and transverse components of momentum.
The fastest wave speed in direction $i$ is thus

$$
\lambda_{\max \text{abs}} \Bigl( \frac{\partial \bm F_{\text{adv}}}{\partial \bm q} \cdot \hat{\bm n}_i \Bigr) = |u_i| + a
$$ (eq-wavespeed)

Note that this wave speed is specific to ideal gases as $\gamma$ is an ideal gas parameter; other equations of state will yield a different acoustic wave speed.

:::

Currently, this demo provides three types of problems/physical models that can be selected at run time via the option `-problem`.
{ref}`problem-advection`, the problem of the transport of energy in a uniform vector velocity field, {ref}`problem-euler-vortex`, the exact solution to the Euler equations, and the so called {ref}`problem-density-current` problem.

(problem-advection)=

## Advection

A simplified version of system {eq}`eq-ns`, only accounting for the transport of total energy, is given by

$$
\frac{\partial E}{\partial t} + \nabla \cdot (\bm{u} E ) = 0 \, ,
$$ (eq-advection)

with $\bm{u}$ the vector velocity field. In this particular test case, a blob of total energy (defined by a characteristic radius $r_c$) is transported by two different wind types.

- **Rotation**

  In this case, a uniform circular velocity field transports the blob of total energy.
  We have solved {eq}`eq-advection` applying zero energy density $E$, and no-flux for $\bm{u}$ on the boundaries.

- **Translation**

  In this case, a background wind with a constant rectilinear velocity field, enters the domain and transports the blob of total energy out of the domain.

  For the inflow boundary conditions, a prescribed $E_{wind}$ is applied weakly on the inflow boundaries such that the weak form boundary integral in {eq}`eq-weak-vector-ns` is defined as

  $$
  \int_{\partial \Omega_{inflow}} \bm v \cdot \bm{F}(\bm q_N) \cdot \widehat{\bm{n}} \,dS = \int_{\partial \Omega_{inflow}} \bm v \, E_{wind} \, \bm u \cdot \widehat{\bm{n}} \,dS  \, ,
  $$

  For the outflow boundary conditions, we have used the current values of $E$, following {cite}`papanastasiou1992outflow` which extends the validity of the weak form of the governing equations to the outflow instead of replacing them with unknown essential or natural boundary conditions.
  The weak form boundary integral in {eq}`eq-weak-vector-ns` for outflow boundary conditions is defined as

  $$
  \int_{\partial \Omega_{outflow}} \bm v \cdot \bm{F}(\bm q_N) \cdot \widehat{\bm{n}} \,dS = \int_{\partial \Omega_{outflow}} \bm v \, E \, \bm u \cdot \widehat{\bm{n}} \,dS  \, ,
  $$

(problem-euler-vortex)=

## Isentropic Vortex

Three-dimensional Euler equations, which are simplified and nondimensionalized version of system {eq}`eq-ns` and account only for the convective fluxes, are given by

$$
\begin{aligned}
\frac{\partial \rho}{\partial t} + \nabla \cdot \bm{U} &= 0 \\
\frac{\partial \bm{U}}{\partial t} + \nabla \cdot \left( \frac{\bm{U} \otimes \bm{U}}{\rho} + P \bm{I}_3 \right) &= 0 \\
\frac{\partial E}{\partial t} + \nabla \cdot \left( \frac{(E + P)\bm{U}}{\rho} \right) &= 0 \, , \\
\end{aligned}
$$ (eq-euler)

Following the setup given in {cite}`zhang2011verification`, the mean flow for this problem is $\rho=1$, $P=1$, $T=P/\rho= 1$ (Specific Gas Constant, $R$, is 1), and $\bm{u}=(u_1,u_2,0)$ while the perturbation $\delta \bm{u}$, and $\delta T$ are defined as

$$
\begin{aligned} (\delta u_1, \, \delta u_2) &= \frac{\epsilon}{2 \pi} \, e^{0.5(1-r^2)} \, (-\bar{y}, \, \bar{x}) \, , \\ \delta T &= - \frac{(\gamma-1) \, \epsilon^2}{8 \, \gamma \, \pi^2} \, e^{1-r^2} \, , \\ \end{aligned}
$$

where $(\bar{x}, \, \bar{y}) = (x-x_c, \, y-y_c)$, $(x_c, \, y_c)$ represents the center of the domain, $r^2=\bar{x}^2 + \bar{y}^2$, and $\epsilon$ is the vortex strength ($\epsilon$ < 10).
There is no perturbation in the entropy $S=P/\rho^\gamma$ ($\delta S=0)$.

(problem-density-current)=

## Density Current

For this test problem (from {cite}`straka1993numerical`), we solve the full Navier-Stokes equations {eq}`eq-ns`, for which a cold air bubble (of radius $r_c$) drops by convection in a neutrally stratified atmosphere.
Its initial condition is defined in terms of the Exner pressure, $\pi(\bm{x},t)$, and potential temperature, $\theta(\bm{x},t)$, that relate to the state variables via

$$
\begin{aligned} \rho &= \frac{P_0}{( c_p - c_v)\theta(\bm{x},t)} \pi(\bm{x},t)^{\frac{c_v}{ c_p - c_v}} \, , \\ e &= c_v \theta(\bm{x},t) \pi(\bm{x},t) + \bm{u}\cdot \bm{u} /2 + g z \, , \end{aligned}
$$

where $P_0$ is the atmospheric pressure.
For this problem, we have used no-slip and non-penetration boundary conditions for $\bm{u}$, and no-flux for mass and energy densities.
## libCEED: Navier-Stokes Example

This page provides a description of the Navier-Stokes example for the libCEED library, based on PETSc.

The Navier-Stokes problem solves the compressible Navier-Stokes equations in three dimensions using an explicit time integration.
The state variables are mass density, momentum density, and energy density.

The main Navier-Stokes solver for libCEED is defined in [`navierstokes.c`](navierstokes.c) with different problem definitions according to the application of interest.

Build by using:

`make`

and run with:

```
./navierstokes -ceed [ceed] -problem [problem type] -degree [degree]
```

## Runtime options

% inclusion-fluids-marker

The Navier-Stokes mini-app is controlled via command-line options.
The following options are common among all problem types:

:::{list-table} Common Runtime Options
:header-rows: 1

* - Option
  - Description
  - Default value

* - `-ceed`
  - CEED resource specifier
  - `/cpu/self/opt/blocked`

* - `-test`
  - Run in test mode
  - `false`

* - `-compare_final_state_atol`
  - Test absolute tolerance
  - `1E-11`

* - `-compare_final_state_filename`
  - Test filename
  -

* - `-problem`
  - Problem to solve (`advection`, `advection2d`, `density_current`, or `euler_vortex`)
  - `density_current`

* - `-implicit`
  - Use implicit time integartor formulation
  -

* - `-degree`
  - Polynomial degree of tensor product basis (must be >= 1)
  - `1`

* - `-qextra`
  - Number of extra quadrature points
  - `2`

* - `-viz_refine`
  - Use regular refinement for visualization
  - `0`

* - `-output_freq`
  - Frequency of output, in number of steps
  - `10`

* - `-continue`
  - Continue from previous solution
  - `0`

* - `-output_dir`
  - Output directory
  - `.`

* - `-dm_plex_box_faces`
  - Number of faces in each linear direction
  - `3,3,3`

* - `-snes_view`
  - View PETSc `SNES` nonlinear solver configuration
  -

* - `-log_view`
  - View PETSc performance log
  -

* - `-help`
  - View comprehensive information about run-time options
  -
:::

For the 2D advection problem, the following additional command-line options are available:

:::{list-table} Advection2D Runtime Options
:header-rows: 1

* - Option
  - Description
  - Default value
  - Unit

* - `-lx`
  - Length scale in x direction
  - `8000`
  - `m`

* - `-ly`
  - Length scale in y direction
  - `8000`
  - `m`

* - `-rc`
  - Characteristic radius of thermal bubble
  - `1000`
  - `m`

* - `-units_meter`
  - 1 meter in scaled length units
  - `1E-2`
  -

* - `-units_second`
  - 1 second in scaled time units
  - `1E-2`
  -

* - `-units_kilogram`
  - 1 kilogram in scaled mass units
  - `1E-6`
  -

* - `-strong_form`
  - Strong (1) or weak/integrated by parts (0) residual
  - `0`
  -

* - `-stab`
  - Stabilization method (`none`, `su`, or `supg`)
  - `none`
  -

* - `-CtauS`
  - Scale coefficient for stabilization tau (nondimensional)
  - `0`
  -

* - `-wind_type`
  - Wind type in Advection (`rotation` or `translation`)
  - `rotation`
  -

* - `-wind_translation`
  - Constant wind vector when `-wind_type translation`
  - `1,0,0`
  -

* - `-E_wind`
  - Total energy of inflow wind when `-wind_type translation`
  - `1E6`
  - `J`
:::

An example of the `rotation` mode can be run with:

```
./navierstokes -problem advection2d -wind_type rotation -implicit -stab supg
```

and the `translation` mode with:

```
./navierstokes -problem advection2d -wind_type translation -wind_translation 1,-.5
```

For the 3D advection problem, the following additional command-line options are available:

:::{list-table} Advection3D Runtime Options
:header-rows: 1

* - Option
  - Description
  - Default value
  - Unit

* - `-lx`
  - Length scale in x direction
  - `8000`
  - `m`

* - `-ly`
  - Length scale in y direction
  - `8000`
  - `m`

* - `-lz`
  - Length scale in z direction
  - `4000`
  - `m`

* - `-rc`
  - Characteristic radius of thermal bubble
  - `1000`
  - `m`

* - `-units_meter`
  - 1 meter in scaled length units
  - `1E-2`
  -

* - `-units_second`
  - 1 second in scaled time units
  - `1E-2`
  -

* - `-units_kilogram`
  - 1 kilogram in scaled mass units
  - `1E-6`
  -

* - `-strong_form`
  - Strong (1) or weak/integrated by parts (0) residual
  - `0`
  -

* - `-stab`
  - Stabilization method (`none`, `su`, or `supg`)
  - `none`
  -

* - `-CtauS`
  - Scale coefficient for stabilization tau (nondimensional)
  - `0`
  -

* - `-wind_type`
  - Wind type in Advection (`rotation` or `translation`)
  - `rotation`
  -

* - `-wind_translation`
  - Constant wind vector when `-wind_type translation`
  - `1,0,0`
  -

* - `-E_wind`
  - Total energy of inflow wind when `-wind_type translation`
  - `1E6`
  - `J`

* - `-bubble_type`
  - `sphere` (3D) or `cylinder` (2D)
  - `shpere`
  -

* - `-bubble_continuity`
  - `smooth`, `back_sharp`, or `thick`
  - `smooth`
  -
:::

An example of the `rotation` mode can be run with:

```
./navierstokes -problem advection -wind_type rotation -implicit -stab supg
```

and the `translation` mode with:

```
./navierstokes -problem advection -wind_type translation -wind_translation .5,-1,0
```

For the Isentropic Vortex problem, the following additional command-line options are available:

:::{list-table} Isentropic Vortex Runtime Options
:header-rows: 1

* - Option
  - Description
  - Default value
  - Unit

* - `-lx`
  - Length scale in x direction
  - `1000`
  - `m`

* - `-ly`
  - Length scale in y direction
  - `1000`
  - `m`

* - `-lz`
  - Length scale in z direction
  - `1`
  - `m`

* - `-center`
  - Location of vortex center
  - `(lx,ly,lz)/2`
  - `(m,m,m)`

* - `-units_meter`
  - 1 meter in scaled length units
  - `1E-2`
  -

* - `-units_second`
  - 1 second in scaled time units
  - `1E-2`
  -

* - `-mean_velocity`
  - Background velocity vector
  - `(1,1,0)`
  -

* - `-vortex_strength`
  - Strength of vortex < 10
  - `5`
  -

* - `-c_tau`
  - Stabilization constant
  - `0.5`
  -
:::

This problem can be run with:

```
./navierstokes -problem euler_vortex -mean_velocity .5,-.8,0.
```

For the Density Current problem, the following additional command-line options are available:

:::{list-table} Euler Vortex Runtime Options
:header-rows: 1

* - Option
  - Description
  - Default value
  - Unit

* - `-lx`
  - Length scale in x direction
  - `8000`
  - `m`

* - `-ly`
  - Length scale in y direction
  - `8000`
  - `m`

* - `-lz`
  - Length scale in z direction
  - `4000`
  - `m`

* - `-center`
  - Location of bubble center
  - `(lx,ly,lz)/2`
  - `(m,m,m)`

* - `-dc_axis`
  - Axis of density current cylindrical anomaly, or `(0,0,0)` for spherically symmetric
  - `(0,0,0)`
  -

* - `-rc`
  - Characteristic radius of thermal bubble
  - `1000`
  - `m`

* - `-bc_wall`
  - Use wall boundary conditions on this list of faces
  - `-`
  -

* - `-bc_slip_x`
  - Use slip boundary conditions, for the x component, on this list of faces
  - `5,6`
  -

* - `-bc_slip_y`
  - Use slip boundary conditions, for the y component, on this list of faces
  - `3,4`
  -

* - `-bc_slip_z`
  - Use slip boundary conditions, for the z component, on this list of faces
  - `1,2`
  -

* - `-units_meter`
  - 1 meter in scaled length units
  - `1E-2`
  -

* - `-units_second`
  - 1 second in scaled time units
  - `1E-2`
  -

* - `-units_kilogram`
  - 1 kilogram in scaled mass units
  - `1E-6`
  -

* - `-units_Kelvin`
  - 1 Kelvin in scaled temperature units
  - `1`
  -

* - `-stab`
  - Stabilization method (`none`, `su`, or `supg`)
  - `none`
  -

* - `-c_tau`
  - Stabilization constant
  - `0.5`
  -

* - `-theta0`
  - Reference potential temperature
  - `300`
  - `K`

* - `-thetaC`
  - Perturbation of potential temperature
  - `-15`
  - `K`

* - `-P0`
  - Atmospheric pressure
  - `1E5`
  - `Pa`

* - `-N`
  - Brunt-Vaisala frequency
  - `0.01`
  - `1/s`

* - `-cv`
  - Heat capacity at constant volume
  - `717`
  - `J/(kg K)`

* - `-cp`
  - Heat capacity at constant pressure
  - `1004`
  - `J/(kg K)`

* - `-g`
  - Gravitational acceleration
  - `9.81`
  - `m/s^2`

* - `-lambda`
  - Stokes hypothesis second viscosity coefficient
  - `-2/3`
  -

* - `-mu`
  - Shear dynamic viscosity coefficient
  - `75`
  -  `Pa s`

* - `-k`
  - Thermal conductivity
  - `0.02638`
  - `W/(m K)`
:::

For the case of a square/cubic mesh, the list of face indices to be used with `-bc_wall` and/or `-bc_slip_x`, `-bc_slip_y`, and `-bc_slip_z` are:

* 2D:
  - faceMarkerBottom = 1
  - faceMarkerRight  = 2
  - faceMarkerTop    = 3
  - faceMarkerLeft   = 4
* 3D:
  - faceMarkerBottom = 1
  - faceMarkerTop    = 2
  - faceMarkerFront  = 3
  - faceMarkerBack   = 4
  - faceMarkerRight  = 5
  - faceMarkerLeft   = 6

This problem can be run with:

```
./navierstokes -problem density_current -dm_plex_box_faces 16,1,8 -degree 1 -lx 2000 -ly 125 -lz 1000 -rc 400. -bc_wall 1,2,5,6 -bc_slip_y 3,4 -viz_refine 2
```
# PETSc demos and BPs

(example-petsc-area)=

## Area

This example is located in the subdirectory {file}`examples/petsc`.
It demonstrates a simple usage of libCEED with PETSc to calculate the surface area of a closed surface.
The code uses higher level communication protocols for mesh handling in PETSc's DMPlex.
This example has the same mathematical formulation as {ref}`Ex1-Volume`, with the exception that the physical coordinates for this problem are $\bm{x}=(x,y,z)\in \mathbb{R}^3$, while the coordinates of the reference element are $\bm{X}=(X,Y) \equiv (X_0,X_1) \in \textrm{I} =[-1,1]^2$.

(example-petsc-area-cube)=

### Cube

This is one of the test cases of the computation of the {ref}`example-petsc-area` of a 2D manifold embedded in 3D.
This problem can be run with:

```
./area -problem cube
```

This example uses the following coordinate transformations for the computation of the geometric factors: from the physical coordinates on the cube, denoted by $\bar{\bm{x}}=(\bar{x},\bar{y},\bar{z})$, and physical coordinates on the discrete surface, denoted by $\bm{{x}}=(x,y)$, to $\bm{X}=(X,Y) \in \textrm{I}$ on the reference element, via the chain rule

$$
\frac{\partial \bm{x}}{\partial \bm{X}}_{(2\times2)} = \frac{\partial {\bm{x}}}{\partial \bar{\bm{x}}}_{(2\times3)} \frac{\partial \bar{\bm{x}}}{\partial \bm{X}}_{(3\times2)},
$$ (eq-coordinate-transforms-cube)

with Jacobian determinant given by

$$
\left| J \right| = \left\|col_1\left(\frac{\partial \bar{\bm{x}}}{\partial \bm{X}}\right)\right\| \left\|col_2 \left(\frac{\partial \bar{\bm{x}}}{\partial \bm{X}}\right) \right\|
$$ (eq-jacobian-cube)

We note that in equation {eq}`eq-coordinate-transforms-cube`, the right-most Jacobian matrix ${\partial\bar{\bm{x}}}/{\partial \bm{X}}_{(3\times2)}$ is provided by the library, while ${\partial{\bm{x}}}/{\partial \bar{ \bm{x}}}_{(2\times3)}$ is provided by the user as

$$
\left[ col_1\left(\frac{\partial\bar{\bm{x}}}{\partial \bm{X}}\right) / \left\| col_1\left(\frac{\partial\bar{\bm{x}}}{\partial \bm{X}}\right)\right\| , col_2\left(\frac{\partial\bar{\bm{x}}}{\partial \bm{X}}\right) / \left\| col_2\left(\frac{\partial\bar{\bm{x}}}{\partial \bm{X}}\right)\right\| \right]^T_{(2\times 3)}.
$$

(example-petsc-area-sphere)=

### Sphere

This problem computes the surface {ref}`example-petsc-area` of a tensor-product discrete sphere, obtained by projecting a cube inscribed in a sphere onto the surface of the sphere.
This discrete surface is sometimes referred to as a cubed-sphere (an example of such as a surface is given in figure {numref}`fig-cubed-sphere`).
This problem can be run with:

```
./area -problem sphere
```

(fig-cubed-sphere)=

:::{figure} ../../../../img/CubedSphere.svg
Example of a cubed-sphere, i.e., a tensor-product discrete sphere, obtained by
projecting a cube inscribed in a sphere onto the surface of the sphere.
:::

This example uses the following coordinate transformations for the computation of the geometric factors: from the physical coordinates on the sphere, denoted by $\overset{\circ}{\bm{x}}=(\overset{\circ}{x},\overset{\circ}{y},\overset{\circ}{z})$, and physical coordinates on the discrete surface, denoted by $\bm{{x}}=(x,y,z)$ (depicted, for simplicity, as coordinates on a circle and 1D linear element in figure {numref}`fig-sphere-coords`), to $\bm{X}=(X,Y) \in \textrm{I}$ on the reference element, via the chain rule

$$
\frac{\partial \overset{\circ}{\bm{x}}}{\partial \bm{X}}_{(3\times2)} = \frac{\partial \overset{\circ}{\bm{x}}}{\partial \bm{x}}_{(3\times3)} \frac{\partial\bm{x}}{\partial \bm{X}}_{(3\times2)} ,
$$ (eq-coordinate-transforms-sphere)

with Jacobian determinant given by

$$
\left| J \right| = \left| col_1\left(\frac{\partial \overset{\circ}{\bm{x}}}{\partial \bm{X}}\right) \times col_2 \left(\frac{\partial \overset{\circ}{\bm{x}}}{\partial \bm{X}}\right)\right| .
$$ (eq-jacobian-sphere)

(fig-sphere-coords)=

:::{figure} ../../../../img/SphereSketch.svg
Sketch of coordinates mapping between a 1D linear element and a circle.
In the case of a linear element the two nodes, $p_0$ and $p_1$, marked by red crosses, coincide with the endpoints of the element.
Two quadrature points, $q_0$ and $q_1$, marked by blue dots, with physical coordinates denoted by $\bm x(\bm X)$, are mapped to their corresponding radial projections on the circle, which have coordinates $\overset{\circ}{\bm{x}}(\bm x)$.
:::

We note that in equation {eq}`eq-coordinate-transforms-sphere`, the right-most Jacobian matrix ${\partial\bm{x}}/{\partial \bm{X}}_{(3\times2)}$ is provided by the library, while ${\partial \overset{\circ}{\bm{x}}}/{\partial \bm{x}}_{(3\times3)}$ is provided by the user with analytical derivatives.
In particular, for a sphere of radius 1, we have

$$
\overset{\circ}{\bm x}(\bm x) = \frac{1}{\lVert \bm x \rVert} \bm x_{(3\times 1)}
$$

and thus

$$
\frac{\partial \overset{\circ}{\bm{x}}}{\partial \bm{x}} = \frac{1}{\lVert \bm x \rVert} \bm I_{(3\times 3)} - \frac{1}{\lVert \bm x \rVert^3} (\bm x \bm x^T)_{(3\times 3)} .
$$

(example-petsc-bps)=

## Bakeoff problems and generalizations

The PETSc examples in this directory include a full suite of parallel {ref}`bakeoff problems <bps>` (BPs) using a "raw" parallel decomposition (see `bpsraw.c`) and using PETSc's `DMPlex` for unstructured grid management (see `bps.c`).
A generalization of these BPs to the surface of the cubed-sphere are available in `bpssphere.c`.

(example-petsc-bps-sphere)=

### Bakeoff problems on the cubed-sphere

For the $L^2$ projection problems, BP1-BP2, that use the mass operator, the coordinate transformations and the corresponding Jacobian determinant, equation {eq}`eq-jacobian-sphere`, are the same as in the {ref}`example-petsc-area-sphere` example.
For the Poisson's problem, BP3-BP6, on the cubed-sphere, in addition to equation {eq}`eq-jacobian-sphere`, the pseudo-inverse of $\partial \overset{\circ}{\bm{x}} / \partial \bm{X}$ is used to derive the contravariant metric tensor (please see figure {numref}`fig-sphere-coords` for a reference of the notation used).
We begin by expressing the Moore-Penrose (left) pseudo-inverse:

$$
\frac{\partial \bm{X}}{\partial \overset{\circ}{\bm{x}}}_{(2\times 3)} \equiv \left(\frac{\partial \overset{\circ}{\bm{x}}}{\partial \bm{X}}\right)_{(2\times 3)}^{+} =  \left(\frac{\partial \overset{\circ}{\bm{x}}}{\partial \bm{X}}_{(2\times3)}^T \frac{\partial\overset{\circ}{\bm{x}}}{\partial \bm{X}}_{(3\times2)} \right)^{-1} \frac{\partial \overset{\circ}{\bm{x}}}{\partial \bm{X}}_{(2\times3)}^T .
$$ (eq-dxcircdX-pseudo-inv)

This enables computation of gradients of an arbitrary function $u(\overset{\circ}{\bm x})$ in the embedding space as

$$
\frac{\partial u}{\partial \overset{\circ}{\bm x}}_{(1\times 3)} = \frac{\partial u}{\partial \bm X}_{(1\times 2)} \frac{\partial \bm X}{\partial \overset{\circ}{\bm x}}_{(2\times 3)}
$$

and thus the weak Laplacian may be expressed as

$$
\int_{\Omega} \frac{\partial v}{\partial \overset\circ{\bm x}} \left( \frac{\partial u}{\partial \overset\circ{\bm x}} \right)^T \, dS
    = \int_{\Omega} \frac{\partial v}{\partial \bm X} \underbrace{\frac{\partial \bm X}{\partial \overset\circ{\bm x}} \left( \frac{\partial \bm X}{\partial \overset\circ{\bm x}} \right)^T}_{\bm g_{(2\times 2)}}  \left(\frac{\partial u}{\partial \bm X} \right)^T \, dS
$$ (eq-weak-laplace-sphere)

where we have identified the $2\times 2$ contravariant metric tensor $\bm g$ (sometimes written $\bm g^{ij}$), and where now $\Omega$ represents the surface of the sphere, which is a two-dimensional closed surface embedded in the three-dimensional Euclidean space $\mathbb{R}^3$.
This expression can be simplified to avoid the explicit Moore-Penrose pseudo-inverse,

$$
\bm g = \left(\frac{\partial \overset{\circ}{\bm{x}}}{\partial \bm{X}}^T \frac{\partial\overset{\circ}{\bm{x}}}{\partial \bm{X}} \right)^{-1}_{(2\times 2)} \frac{\partial \overset{\circ}{\bm{x}}}{\partial \bm{X}}_{(2\times3)}^T \frac{\partial \overset{\circ}{\bm{x}}}{\partial \bm{X}}_{(3\times2)} \left(\frac{\partial \overset{\circ}{\bm{x}}}{\partial \bm{X}}^T \frac{\partial\overset{\circ}{\bm{x}}}{\partial \bm{X}} \right)^{-T}_{(2\times 2)} = \left(\frac{\partial \overset{\circ}{\bm{x}}}{\partial \bm{X}}^T \frac{\partial\overset{\circ}{\bm{x}}}{\partial \bm{X}} \right)^{-1}_{(2\times 2)}
$$

where we have dropped the transpose due to symmetry.
This allows us to simplify {eq}`eq-weak-laplace-sphere` as

$$
\int_{\Omega} \frac{\partial v}{\partial \overset\circ{\bm x}} \left( \frac{\partial u}{\partial \overset\circ{\bm x}} \right)^T \, dS     = \int_{\Omega} \frac{\partial v}{\partial \bm X} \underbrace{\left(\frac{\partial \overset{\circ}{\bm{x}}}{\partial \bm{X}}^T \frac{\partial\overset{\circ}{\bm{x}}}{\partial \bm{X}} \right)^{-1}}_{\bm g_{(2\times 2)}}  \left(\frac{\partial u}{\partial \bm X} \right)^T \, dS ,
$$

which is the form implemented in `qfunctions/bps/bp3sphere.h`.

(example-petsc-multigrid)=

## Multigrid

This example is located in the subdirectory {file}`examples/petsc`.
It investigates $p$-multigrid for the Poisson problem, equation {eq}`eq-variable-coeff-poisson`, using an unstructured high-order finite element discretization.
All of the operators associated with the geometric multigrid are implemented in libCEED.

$$
-\nabla\cdot \left( \kappa \left( x \right) \nabla x \right) = g \left( x \right)
$$ (eq-variable-coeff-poisson)

The Poisson operator can be specified with the decomposition given by the equation in figure {ref}`fig-operator-decomp`, and the restriction and prolongation operators given by interpolation basis operations, $\bm{B}$, and $\bm{B}^T$, respectively, act on the different grid levels with corresponding element restrictions, $\bm{G}$.
These three operations can be exploited by existing matrix-free multigrid software and smoothers.
Preconditioning based on the libCEED finite element operator decomposition is an ongoing area of research.
## libCEED + PETSc Examples

### CEED bakeoff problems with raw mesh management - bpsraw

This code solves the CEED bakeoff problems on a structured grid generated and referenced using only low-level communication primitives.

To build, run `make bpsraw`

To run, `./bpsraw -ceed [ceed-resource] -problem bp[1-6] -degree [degree]`

In addition to the common arguments, the following arguments may be set:

- `-local`             - Target number of locally owned DoFs per process

### CEED bakeoff problems with DMPlex - bps

This code solves the CEED bakeoff problems on a unstructured grid using DMPlex.
This example requires a PETSc version later than 3.11.3.

To build, run `make bps`

To run, `./bps -ceed [ceed-resource] -problem bp[1-6] -degree [degree]`

In addition to the common arguments, the following arguments may be set:

- `-mesh`              - Read mesh from file
- `-cells`             - Number of cells per dimension

#### Running a suite

Some run-time arguments can be passed lists, which allows a single `mpiexec` invocation to run many experiments.
For example

    mpiexec -n 64 ./bps -problem bp1,bp2,bp3,bp4 -degree 2,3,5,7                \
        -ceed /cpu/self/opt/serial,/cpu/self/xsmm/serial,/cpu/self/xsmm/blocked \
        -local_nodes 600,20000 | tee bps.log

which will sample from the `4*4*3=48` specified combinations, each of which will run a problem-size sweep of 600, 1200, 2400, 4800, 9600, 192000 FEM nodes per MPI rank. 
The resulting log file can be read by the Python plotting scripts in `benchmarks/`.

### CEED bakeoff problems with DMPlex and PCMG - multigrid

This code solves the CEED bakeoff problems on a unstructured grid using DMPlex with p-multigrid implemented in PCMG.
This example requires a PETSc version later than 3.11.3.

To build, run `make multigrid`

To run, `./multigrid -ceed [ceed-resource] -problem bp[1-6] -degree [degree]`

In addition to the common arguments, the following arguments may be set:

- `-mesh`              - Read mesh from file
- `-cells`             - Number of cells per dimension

### Command line arguments

The following arguments can be specified for all of the above examples:

- `-ceed`              - CEED resource specifier
- `-problem`           - CEED benchmark problem to solve
- `-degree`            - Polynomial degree of tensor product basis
- `-qextra`            - Number of extra quadrature points
- `-test`              - Testing mode (do not print unless error is large)
- `-benchmark`         - Benchmarking mode (prints benchmark statistics)

### libCEED example to compute surface area using DMPlex - area

This example uses the mass matrix to compute the surface area of a cube or a discrete cubed-sphere, defined via DMPlex.

To build, run `make area`

To run, `./area -problem cube -ceed [ceed-resource] -petscspace_degree [degree]`

or

`./area -problem sphere -ceed [ceed-resource] -petscspace_degree [degree]`

#### Command line arguments

The following arguments can be specified for the area example:

- `-ceed`              - CEED resource specifier
- `-problem`           - Problem to solve, either 'cube' or 'sphere'
- `-petscspace_degree` - Polynomial degree of tensor product basis
- `-qextra`            - Number of extra quadrature points
- `-test`              - Testing mode (do not print unless error is large)
- `-mesh`              - Read mesh from file

## libCEED + Nek5000 Examples

### Prerequisites

Nek5000 v18.0 or greater must be [installed](https://nek5000.mcs.anl.gov/getstarted/) to run
these examples.  It is assumed to exist at `../../../Nek5000` (a sibling to the
libCEED directory) or at a path defined in the environment variable `NEK5K_DIR`.
For example, you could set
```sh
    export NEK5K_DIR=/scratch/Nek5000
```
if that is where it is located.

The Nek5000 examples depend on the Nek5000 tools: `genbox`, `genmap`,
and `reatore2`. They can be built using
```sh
   ( cd $NEK5K_DIR/tools && ./maketools genbox genmap reatore2 )
```
See also the [Nek5000 documentation](https://nek5000.mcs.anl.gov/getstarted/).

### Building the Nek5000 examples

You can build the Nek5000 libCEED examples with the command `make bps`.

You can also build the Nek5000 libCEED examples by invoking `nek-examples.sh` script.
```sh
  ./nek-examples.sh -m
```

By default, the examples are built with MPI. To build the examples without MPI,
set the environment variable `MPI=0`.

Note: Nek5000 examples must be built sequentially. Due to the Nek5000 build
process, multiple examples cannot be built in parallel. At present, there is
only one Nek5000 example file to build, which handles both CEED BP 1 and
CEED BP 3.

### Running Nek5000 examples

You can run the Nek5000 libCEED examples by invoking `nek-examples.sh`
script. The syntax is:
```sh
  ./nek-examples.sh -c <ceed_backend> -e <example_name> \
                   -n <mpi_ranks> -b <box_geometry>
```
The different options that can be used for the script are listed below:
```
options:
   -h|-help     Print this usage information and exit
   -c|-ceed     Ceed backend to be used for the run (optional, default: /cpu/self)
   -e|-example  Example name (optional, default: bp1)
   -n|-np       Specify number of MPI ranks for the run (optional, default: 1)
   -t|-test     Run in test mode (not on by default)
   -b|-box      Box case in boxes sub-directory found along with this script (default: 2x2x2)
   -clean       clean the examples directory
   -m|-make     Make the examples
```
The only mandatory argument is `-b` or `-box` which sets the box geometry to be
used. This geometry should be found in `./boxes` directory.

For example, you can run bp1 as follows:
```sh
  ./run-nek-example.sh -ceed /cpu/self -e bp1 -n 4 -b 3
```
which is the same as running:
```sh
  ./run-nek-example.sh -b 3
```
(example-petsc-elasticity)=

# Solid mechanics mini-app

This example is located in the subdirectory {file}`examples/solids`.
It solves the steady-state static momentum balance equations using unstructured high-order finite/spectral element spatial discretizations.
As for the {ref}`example-petsc-navier-stokes` case, the solid mechanics elasticity example has been developed using PETSc, so that the pointwise physics (defined at quadrature points) is separated from the parallelization and meshing concerns.

In this mini-app, we consider three formulations used in solid mechanics applications: linear elasticity, Neo-Hookean hyperelasticity at small strain, and Neo-Hookean hyperelasticity at finite strain.
We provide the strong and weak forms of static balance of linear momentum in the small strain and finite strain regimes.
The stress-strain relationship (constitutive law) for each of the material models is provided.
Due to the nonlinearity of material models in Neo-Hookean hyperelasticity, the Newton linearization of the material models is provided.

:::{note}
Linear elasticity and small-strain hyperelasticity can both by obtained from the finite-strain hyperelastic formulation by linearization of geometric and constitutive nonlinearities.
The effect of these linearizations is sketched in the diagram below, where $\bm \sigma$ and $\bm \epsilon$ are stress and strain, respectively, in the small strain regime, while $\bm S$ and $\bm E$ are their finite-strain generalizations (second Piola-Kirchoff tensor and Green-Lagrange strain tensor, respectively) defined in the initial configuration, and $\mathsf C$ is a linearized constitutive model.

$$
\begin{CD}
  {\overbrace{\bm S(\bm E)}^{\text{Finite Strain Hyperelastic}}}
  @>{\text{constitutive}}>{\text{linearization}}>
  {\overbrace{\bm S = \mathsf C \bm E}^{\text{St. Venant-Kirchoff}}} \\
  @V{\text{geometric}}V{\begin{smallmatrix}\bm E \to \bm \epsilon \\ \bm S \to \bm \sigma \end{smallmatrix}}V
  @V{\begin{smallmatrix}\bm E \to \bm \epsilon \\ \bm S \to \bm \sigma \end{smallmatrix}}V{\text{geometric}}V \\
  {\underbrace{\bm \sigma(\bm \epsilon)}_\text{Small Strain Hyperelastic}}
  @>{\text{constitutive}}>\text{linearization}>
  {\underbrace{\bm \sigma = \mathsf C \bm \epsilon}_\text{Linear Elastic}}
\end{CD}
$$ (hyperelastic-cd)
:::

(running-elasticity)=

## Running the mini-app

```{include} README.md
:start-after: inclusion-solids-marker
```

(problem-linear-elasticity)=

## Linear Elasticity

The strong form of the static balance of linear momentum at small strain for the three-dimensional linear elasticity problem is given by {cite}`hughes2012finite`:

$$
\nabla \cdot \bm{\sigma} + \bm{g} = \bm{0}
$$ (lin-elas)

where $\bm{\sigma}$ and $\bm{g}$ are stress and forcing functions, respectively.
We multiply {eq}`lin-elas` by a test function $\bm v$ and integrate the divergence term by parts to arrive at the weak form: find $\bm u \in \mathcal V \subset H^1(\Omega)$ such that

$$
\int_{\Omega}{ \nabla \bm{v} \tcolon \bm{\sigma}} \, dV
- \int_{\partial \Omega}{\bm{v} \cdot \left(\bm{\sigma} \cdot \hat{\bm{n}}\right)} \, dS
- \int_{\Omega}{\bm{v} \cdot \bm{g}} \, dV
= 0, \quad \forall \bm v \in \mathcal V,
$$ (lin-elas-weak)

where $\bm{\sigma} \cdot \hat{\bm{n}}|_{\partial \Omega}$ is replaced by an applied force/traction boundary condition written in terms of the initial configuration.
When inhomogeneous Dirichlet boundary conditions are present, $\mathcal V$ is an affine space that satisfies those boundary conditions.

### Constitutive modeling

In their most general form, constitutive models define $\bm \sigma$ in terms of state variables.
In the model taken into consideration in the present mini-app, the state variables are constituted by the vector displacement field $\bm u$, and its gradient $\nabla \bm u$.
We begin by defining the symmetric (small/infintesimal) strain tensor as

$$
\bm{\epsilon} = \dfrac{1}{2}\left(\nabla \bm{u} + \nabla \bm{u}^T \right).
$$ (small-strain)

This constitutive model $\bm \sigma(\bm \epsilon)$ is a linear tensor-valued function of a tensor-valued input, but we will consider the more general nonlinear case in other models below.
In these cases, an arbitrary choice of such a function will generally not be invariant under orthogonal transformations and thus will not admissible as a physical model must not depend on the coordinate system chosen to express it.
In particular, given an orthogonal transformation $Q$, we desire

$$
Q \bm \sigma(\bm \epsilon) Q^T = \bm \sigma(Q \bm \epsilon Q^T),
$$ (elastic-invariance)

which means that we can change our reference frame before or after computing $\bm \sigma$, and get the same result either way.
Constitutive relations in which $\bm \sigma$ is uniquely determined by $\bm \epsilon$ while satisfying the invariance property {eq}`elastic-invariance` are known as Cauchy elastic materials.
Here, we define a strain energy density functional $\Phi(\bm \epsilon) \in \mathbb R$ and obtain the strain energy from its gradient,

$$
\bm \sigma(\bm \epsilon) = \frac{\partial \Phi}{\partial \bm \epsilon}.
$$ (strain-energy-grad)

:::{note}
The strain energy density functional cannot be an arbitrary function $\Phi(\bm \epsilon)$; it can only depend on *invariants*, scalar-valued functions $\gamma$ satisfying

$$
\gamma(\bm \epsilon) = \gamma(Q \bm \epsilon Q^T)
$$

for all orthogonal matrices $Q$.
:::

For the linear elasticity model, the strain energy density is given by

$$
\bm{\Phi} = \frac{\lambda}{2} (\operatorname{trace} \bm{\epsilon})^2 + \mu \bm{\epsilon} : \bm{\epsilon} .
$$

The constitutive law (stress-strain relationship) is therefore given by its gradient,

$$
\bm\sigma = \lambda (\operatorname{trace} \bm\epsilon) \bm I_3 + 2 \mu \bm\epsilon,
$$

where $\bm I_3$ is the $3 \times 3$ identity matrix, the colon represents a double contraction (over both indices of $\bm \epsilon$), and the Lamé parameters are given by

$$
\begin{aligned} \lambda &= \frac{E \nu}{(1 + \nu)(1 - 2 \nu)} \\ \mu &= \frac{E}{2(1 + \nu)} \end{aligned}.
$$

The constitutive law (stress-strain relationship) can also be written as

$$
\bm{\sigma} = \mathsf{C} \!:\! \bm{\epsilon}.
$$ (linear-stress-strain)

For notational convenience, we express the symmetric second order tensors $\bm \sigma$ and $\bm \epsilon$ as vectors of length 6 using the [Voigt notation](https://en.wikipedia.org/wiki/Voigt_notation).
Hence, the fourth order elasticity tensor $\mathsf C$ (also known as elastic moduli tensor or material stiffness tensor) can be represented as

$$
\mathsf C = \begin{pmatrix}
\lambda + 2\mu & \lambda & \lambda & & & \\
\lambda & \lambda + 2\mu & \lambda & & & \\
\lambda & \lambda & \lambda + 2\mu & & & \\
& & & \mu & & \\
& & & & \mu & \\
& & & & & \mu
\end{pmatrix}.
$$ (linear-elasticity-tensor)

Note that the incompressible limit $\nu \to \frac 1 2$ causes $\lambda \to \infty$, and thus $\mathsf C$ becomes singular.

(problem-hyper-small-strain)=

## Hyperelasticity at Small Strain

The strong and weak forms given above, in {eq}`lin-elas` and {eq}`lin-elas-weak`, are valid for Neo-Hookean hyperelasticity at small strain.
However, the strain energy density differs and is given by

$$
\bm{\Phi} = \lambda (1 + \operatorname{trace} \bm{\epsilon}) (\log(1 + \operatorname{trace} \bm\epsilon) - 1) + \mu \bm{\epsilon} : \bm{\epsilon} .
$$

As above, we have the corresponding constitutive law given by

$$
\bm{\sigma} = \lambda \log(1 + \operatorname{trace} \bm\epsilon) \bm{I}_3 + 2\mu \bm{\epsilon}
$$ (eq-neo-hookean-small-strain)

where $\bm{\epsilon}$ is defined as in {eq}`small-strain`.

### Newton linearization

Due to nonlinearity in the constitutive law, we require a Newton linearization of {eq}`eq-neo-hookean-small-strain`.
To derive the Newton linearization, we begin by expressing the derivative,

$$
\diff \bm{\sigma} = \dfrac{\partial \bm{\sigma}}{\partial \bm{\epsilon}} \tcolon \diff \bm{\epsilon}
$$

where

$$
\diff \bm{\epsilon} = \dfrac{1}{2}\left( \nabla \diff \bm{u} + \nabla \diff \bm{u}^T \right)
$$

and

$$
\diff \nabla \bm{u} = \nabla \diff \bm{u} .
$$

Therefore,

$$
\diff \bm{\sigma}  = \bar{\lambda} \cdot \operatorname{trace} \diff \bm{\epsilon} \cdot \bm{I}_3 + 2\mu \diff \bm{\epsilon}
$$ (derss)

where we have introduced the symbol

$$
\bar{\lambda} = \dfrac{\lambda}{1 + \epsilon_v }
$$

where volumetric strain is given by $\epsilon_v = \sum_i \epsilon_{ii}$.

Equation {eq}`derss` can be written in Voigt matrix notation as follows:

$$
\begin{pmatrix}
  \diff \sigma_{11} \\
  \diff \sigma_{22} \\
  \diff \sigma_{33} \\
  \diff \sigma_{23} \\
  \diff \sigma_{13} \\
  \diff \sigma_{12}
\end{pmatrix}  =
\begin{pmatrix}
  2 \mu +\bar{\lambda} & \bar{\lambda} & \bar{\lambda} & & & \\
  \bar{\lambda} & 2 \mu +\bar{\lambda} & \bar{\lambda} & & & \\
  \bar{\lambda} & \bar{\lambda} & 2 \mu +\bar{\lambda} & & & \\
  & & & \mu & & \\
  & & & & \mu & \\
  & & & & & \mu \\
\end{pmatrix}
\begin{pmatrix}
  \diff \epsilon_{11} \\
  \diff \epsilon_{22} \\
  \diff \epsilon_{33} \\
  2 \diff \epsilon_{23} \\
  2 \diff \epsilon_{13} \\
  2 \diff \epsilon_{12}
\end{pmatrix}.
$$ (mdss)

(problem-hyperelasticity-finite-strain)=

## Hyperelasticity at Finite Strain

In the *total Lagrangian* approach for the Neo-Hookean hyperelasticity problem, the discrete equations are formulated with respect to the initial configuration.
In this formulation, we solve for displacement $\bm u(\bm X)$ in the reference frame $\bm X$.
The notation for elasticity at finite strain is inspired by {cite}`holzapfel2000nonlinear` to distinguish between the current and initial configurations.
As explained in the {ref}`common-notation` section, we denote by capital letters the reference frame and by small letters the current one.

The strong form of the static balance of linear-momentum at *finite strain* (total Lagrangian) is given by:

$$
- \nabla_X \cdot \bm{P} - \rho_0 \bm{g} = \bm{0}
$$ (sblFinS)

where the $_X$ in $\nabla_X$ indicates that the gradient is calculated with respect to the initial configuration in the finite strain regime.
$\bm{P}$ and $\bm{g}$ are the *first Piola-Kirchhoff stress* tensor and the prescribed forcing function, respectively.
$\rho_0$ is known as the *initial* mass density.
The tensor $\bm P$ is not symmetric, living in the current configuration on the left and the initial configuration on the right.

$\bm{P}$ can be decomposed as

$$
\bm{P} = \bm{F} \, \bm{S},
$$ (1st2nd)

where $\bm S$ is the *second Piola-Kirchhoff stress* tensor, a symmetric tensor defined entirely in the initial configuration, and $\bm{F} = \bm I_3 + \nabla_X \bm u$ is the deformation gradient.
Different constitutive models can define $\bm S$.

### Constitutive modeling

For the constitutive modeling of hyperelasticity at finite strain, we begin by defining two symmetric tensors in the initial configuration, the right Cauchy-Green tensor

$$
\bm C = \bm F^T \bm F
$$

and the Green-Lagrange strain tensor

$$
\bm E = \frac 1 2 (\bm C - \bm I_3) = \frac 1 2 \Big( \nabla_X \bm u + (\nabla_X \bm u)^T + (\nabla_X \bm u)^T \nabla_X \bm u \Big),
$$ (eq-green-lagrange-strain)

the latter of which converges to the linear strain tensor $\bm \epsilon$ in the small-deformation limit.
The constitutive models considered, appropriate for large deformations, express $\bm S$ as a function of $\bm E$, similar to the linear case, shown in equation  {eq}`linear-stress-strain`, which  expresses the relationship between $\bm\sigma$ and $\bm\epsilon$.

Recall that the strain energy density functional can only depend upon invariants.
We will assume without loss of generality that $\bm E$ is diagonal and take its set of eigenvalues as the invariants.
It is clear that there can be only three invariants, and there are many alternate choices, such as $\operatorname{trace}(\bm E), \operatorname{trace}(\bm E^2), \lvert \bm E \rvert$, and combinations thereof.
It is common in the literature for invariants to be taken from $\bm C = \bm I_3 + 2 \bm E$ instead of $\bm E$.

For example, if we take the compressible Neo-Hookean model,

$$
\begin{aligned}
\Phi(\bm E) &= \frac{\lambda}{2}(\log J)^2 - \mu \log J + \frac \mu 2 (\operatorname{trace} \bm C - 3) \\
  &= \frac{\lambda}{2}(\log J)^2 - \mu \log J + \mu \operatorname{trace} \bm E,
\end{aligned}
$$ (neo-hookean-energy)

where $J = \lvert \bm F \rvert = \sqrt{\lvert \bm C \rvert}$ is the determinant of deformation (i.e., volume change) and $\lambda$ and $\mu$ are the Lamé parameters in the infinitesimal strain limit.

To evaluate {eq}`strain-energy-grad`, we make use of

$$
\frac{\partial J}{\partial \bm E} = \frac{\partial \sqrt{\lvert \bm C \rvert}}{\partial \bm E} = \lvert \bm C \rvert^{-1/2} \lvert \bm C \rvert \bm C^{-1} = J \bm C^{-1},
$$

where the factor of $\frac 1 2$ has been absorbed due to $\bm C = \bm I_3 + 2 \bm E.$
Carrying through the differentiation {eq}`strain-energy-grad` for the model {eq}`neo-hookean-energy`, we arrive at

$$
\bm S = \lambda \log J \bm C^{-1} + \mu (\bm I_3 - \bm C^{-1}).
$$ (neo-hookean-stress)

:::{tip}
An equivalent form of {eq}`neo-hookean-stress` is

$$
\bm S = \lambda \log J \bm C^{-1} + 2 \mu \bm C^{-1} \bm E,
$$ (neo-hookean-stress-stable)

which is more numerically stable for small $\bm E$, and thus preferred for computation.
Note that the product $\bm C^{-1} \bm E$ is also symmetric, and that $\bm E$ should be computed using {eq}`eq-green-lagrange-strain`.

Similarly, it is preferable to compute $\log J$ using `log1p`, especially in case of nearly incompressible materials.
To sketch this idea, suppose we have the $2\times 2$ non-symmetric matrix $\bm{F} = \left( \begin{smallmatrix} 1 + u_{0,0} & u_{0,1} \\ u_{1,0} & 1 + u_{1,1} \end{smallmatrix} \right)$.
Then we compute

$$
\log J = \mathtt{log1p}(u_{0,0} + u_{1,1} + u_{0,0} u_{1,1} - u_{0,1} u_{1,0}),
$$ (log1p)

which gives accurate results even in the limit when the entries $u_{i,j}$ are very small.
For example, if $u_{i,j} \sim 10^{-8}$, then naive computation of $\bm I_3 - \bm C^{-1}$ and $\log J$ will have a relative accuracy of order $10^{-8}$ in double precision and no correct digits in single precision.
When using the stable choices above, these quantities retain full $\varepsilon_{\text{machine}}$ relative accuracy.
:::

:::{dropdown} Mooney-Rivlin model
While the Neo-Hookean model depends on just two scalar invariants, $\mathbb I_1 = \trace \bm C = 3 + 2\trace \bm E$ and $J$, Mooney-Rivlin models depend on the additional invariant, $\mathbb I_2 = \frac 1 2 (\mathbb I_1^2 - \bm C \tcolon \bm C)$.
A coupled Mooney-Rivlin strain energy density (cf. Neo-Hookean {eq}`neo-hookean-energy`) is {cite}`holzapfel2000nonlinear`

$$
\Phi(\mathbb{I_1}, \mathbb{I_2}, J) = \frac{\lambda}{2}(\log J)^2 - (\mu_1 + 2\mu_2) \log J + \frac{\mu_1}{2}(\mathbb{I_1} - 3) + \frac{\mu_2}{2}(\mathbb{I_2} - 3).
$$ (mooney-rivlin-energy_coupled)

We differentiate $\Phi$ as in the Neo-Hookean case {eq}`neo-hookean-stress` to yield the second Piola-Kirchoff tensor,

$$
\begin{aligned}
\bm S &=  \lambda \log J \bm{C}^{-1} - (\mu_1 + 2\mu_2) \bm{C}^{-1} + \mu_1\bm I_3 + \mu_2(\mathbb{I_1} \bm I_3 - \bm C) \\
&= (\lambda \log J - \mu_1 - 2\mu_2) \bm C^{-1} + (\mu_1 + \mu_2 \mathbb I_1) \bm I_3 - \mu_2 \bm C,
\end{aligned}
$$ (mooney-rivlin-stress_coupled)

where we have used

$$
\begin{aligned}
\frac{\partial \mathbb{I_1}}{\partial \bm E} &= 2 \bm I_3, & \frac{\partial \mathbb{I_2}}{\partial \bm E} &= 2 \mathbb I_1 \bm I_3 - 2 \bm C, & \frac{\partial \log J}{\partial \bm E} &= \bm{C}^{-1}.
\end{aligned}
$$ (None)

This is a common model for vulcanized rubber, with a shear modulus (defined for the small-strain limit) of $\mu_1 + \mu_2$ that should be significantly smaller than the first Lamé parameter $\lambda$.
:::

:::{dropdown} Mooney-Rivlin strain energy comparison
We apply traction to a block and plot integrated strain energy $\Phi$ as a function of the loading paramater.

```{altair-plot}
:hide-code:

import altair as alt
import pandas as pd
def source_path(rel):
    import os
    return os.path.join(os.path.dirname(os.environ["DOCUTILSCONFIG"]), rel)

nh = pd.read_csv(source_path("examples/solids/tests-output/NH-strain.csv"))
nh["model"] = "Neo-Hookean"
nh["parameters"] = "E=2.8, nu=0.4"

mr = pd.read_csv(source_path("examples/solids/tests-output/MR-strain.csv"))
mr["model"] = "Mooney-Rivlin; Neo-Hookean equivalent"
mr["parameters"] = "mu_1=1, mu_2=0, nu=.4"

mr1 = pd.read_csv(source_path("examples/solids/tests-output/MR-strain1.csv"))
mr1["model"] = "Mooney-Rivlin"
mr1["parameters"] = "mu_1=0.5, mu_2=0.5, nu=.4"

df = pd.concat([nh, mr, mr1])
highlight = alt.selection_single(
   on = "mouseover",
   nearest = True,
   fields=["model", "parameters"],
)
base = alt.Chart(df).encode(
   alt.X("increment"),
   alt.Y("energy", scale=alt.Scale(type="sqrt")),
   alt.Color("model"),
   alt.Tooltip(("model", "parameters")),
   opacity=alt.condition(highlight, alt.value(1), alt.value(.5)),
   size=alt.condition(highlight, alt.value(2), alt.value(1)),
)
base.mark_point().add_selection(highlight) + base.mark_line()
```
:::

:::{note}
One can linearize {eq}`neo-hookean-stress` around $\bm E = 0$, for which $\bm C = \bm I_3 + 2 \bm E \to \bm I_3$ and $J \to 1 + \operatorname{trace} \bm E$, therefore {eq}`neo-hookean-stress` reduces to

$$
\bm S = \lambda (\trace \bm E) \bm I_3 + 2 \mu \bm E,
$$ (eq-st-venant-kirchoff)

which is the St. Venant-Kirchoff model (constitutive linearization without geometric linearization; see {eq}`hyperelastic-cd`).

This model can be used for geometrically nonlinear mechanics (e.g., snap-through of thin structures), but is inappropriate for large strain.

Alternatively, one can drop geometric nonlinearities, $\bm E \to \bm \epsilon$ and $\bm C \to \bm I_3$, while retaining the nonlinear dependence on $J \to 1 + \operatorname{trace} \bm \epsilon$, thereby yielding {eq}`eq-neo-hookean-small-strain` (see {eq}`hyperelastic-cd`).
:::

### Weak form

We multiply {eq}`sblFinS` by a test function $\bm v$ and integrate by parts to obtain the weak form for finite-strain hyperelasticity:
find $\bm u \in \mathcal V \subset H^1(\Omega_0)$ such that

$$
\int_{\Omega_0}{\nabla_X \bm{v} \tcolon \bm{P}} \, dV
 - \int_{\Omega_0}{\bm{v} \cdot \rho_0 \bm{g}} \, dV
 - \int_{\partial \Omega_0}{\bm{v} \cdot (\bm{P} \cdot \hat{\bm{N}})} \, dS
 = 0, \quad \forall \bm v \in \mathcal V,
$$ (hyperelastic-weak-form-initial)

where $\bm{P} \cdot \hat{\bm{N}}|_{\partial\Omega}$ is replaced by any prescribed force/traction boundary condition written in terms of the initial configuration.
This equation contains material/constitutive nonlinearities in defining $\bm S(\bm E)$, as well as geometric nonlinearities through $\bm P = \bm F\, \bm S$, $\bm E(\bm F)$, and the body force $\bm g$, which must be pulled back from the current configuration to the initial configuration.
Discretization of {eq}`hyperelastic-weak-form-initial` produces a finite-dimensional system of nonlinear algebraic equations, which we solve using Newton-Raphson methods.
One attractive feature of Galerkin discretization is that we can arrive at the same linear system by discretizing the Newton linearization of the continuous form; that is, discretization and differentiation (Newton linearization) commute.

### Newton linearization

To derive a Newton linearization of {eq}`hyperelastic-weak-form-initial`, we begin by expressing the derivative of {eq}`1st2nd` in incremental form,

$$
\diff \bm P = \frac{\partial \bm P}{\partial \bm F} \!:\! \diff \bm F = \diff \bm F\, \bm S + \bm F \underbrace{\frac{\partial \bm S}{\partial \bm E} \!:\! \diff \bm E}_{\diff \bm S}
$$ (eq-diff-P)

where

$$
\diff \bm E = \frac{\partial \bm E}{\partial \bm F} \!:\! \diff \bm F = \frac 1 2 \Big( \diff \bm F^T \bm F + \bm F^T \diff \bm F \Big)
$$

and $\diff\bm F = \nabla_X\diff\bm u$.
The quantity ${\partial \bm S} / {\partial \bm E}$ is known as the incremental elasticity tensor, and is analogous to the linear elasticity tensor $\mathsf C$ of {eq}`linear-elasticity-tensor`.
We now evaluate $\diff \bm S$ for the Neo-Hookean model {eq}`neo-hookean-stress`,

$$
\diff\bm S = \frac{\partial \bm S}{\partial \bm E} \!:\! \diff \bm E
= \lambda (\bm C^{-1} \!:\! \diff\bm E) \bm C^{-1}
  + 2 (\mu - \lambda \log J) \bm C^{-1} \diff\bm E \, \bm C^{-1},
$$ (eq-neo-hookean-incremental-stress)

where we have used

$$
\diff \bm C^{-1} = \frac{\partial \bm C^{-1}}{\partial \bm E} \!:\! \diff\bm E = -2 \bm C^{-1} \diff \bm E \, \bm C^{-1} .
$$

:::{note}
In the small-strain limit, $\bm C \to \bm I_3$ and $\log J \to 0$, thereby reducing {eq}`eq-neo-hookean-incremental-stress` to the St. Venant-Kirchoff model {eq}`eq-st-venant-kirchoff`.
:::

:::{dropdown} Newton linearization of Mooney-Rivlin
Similar to {eq}`eq-neo-hookean-incremental-stress`, we differentiate {eq}`mooney-rivlin-stress_coupled` using variational notation,

$$
\begin{aligned}
\diff\bm S &= \lambda (\bm C^{-1} \tcolon \diff\bm E) \bm C^{-1} \\
&\quad + 2(\mu_1 + 2\mu_2 - \lambda \log J) \bm C^{-1} \diff\bm E \bm C^{-1} \\
&\quad + 2 \mu_2 \Big[ \trace (\diff\bm E) \bm I_3 - \diff\bm E\Big] .
\end{aligned}
$$ (mooney-rivlin-dS-coupled)

Note that this agrees with {eq}`eq-neo-hookean-incremental-stress` if $\mu_1 = \mu, \mu_2 = 0$.
Moving from Neo-Hookean to Mooney-Rivlin modifies the second term and adds the third.
:::

:::{dropdown} Cancellation vs symmetry
Some cancellation is possible (at the expense of symmetry) if we substitute {eq}`eq-neo-hookean-incremental-stress` into {eq}`eq-diff-P`,

$$
\begin{aligned}
\diff \bm P &= \diff \bm F\, \bm S
  + \lambda (\bm C^{-1} : \diff \bm E) \bm F^{-T} + 2(\mu - \lambda \log J) \bm F^{-T} \diff\bm E \, \bm C^{-1} \\
&= \diff \bm F\, \bm S
  + \lambda (\bm F^{-T} : \diff \bm F) \bm F^{-T} + (\mu - \lambda \log J) \bm F^{-T} (\bm F^T \diff \bm F + \diff \bm F^T \bm F) \bm C^{-1} \\
&= \diff \bm F\, \bm S
  + \lambda (\bm F^{-T} : \diff \bm F) \bm F^{-T} + (\mu - \lambda \log J) \Big( \diff \bm F\, \bm C^{-1} + \bm F^{-T} \diff \bm F^T \bm F^{-T} \Big),
\end{aligned}
$$ (eq-diff-P-dF)

where we have exploited $\bm F \bm C^{-1} = \bm F^{-T}$ and

$$
\begin{aligned} \bm C^{-1} \!:\! \diff \bm E = \bm C_{IJ}^{-1} \diff \bm E_{IJ} &= \frac 1 2 \bm F_{Ik}^{-1} \bm F_{Jk}^{-1} (\bm F_{\ell I} \diff \bm F_{\ell J} + \diff \bm F_{\ell I} \bm F_{\ell J}) \\ &= \frac 1 2 \Big( \delta_{\ell k} \bm F_{Jk}^{-1} \diff \bm F_{\ell J} + \delta_{\ell k} \bm F_{Ik}^{-1} \diff \bm F_{\ell I} \Big) \\ &= \bm F_{Ik}^{-1} \diff \bm F_{kI} = \bm F^{-T} \!:\! \diff \bm F. \end{aligned}
$$

We prefer to compute with {eq}`eq-neo-hookean-incremental-stress` because {eq}`eq-diff-P-dF` is more expensive, requiring access to (non-symmetric) $\bm F^{-1}$ in addition to (symmetric) $\bm C^{-1} = \bm F^{-1} \bm F^{-T}$, having fewer symmetries to exploit in contractions, and being less numerically stable.
:::

:::{dropdown} $\diff\bm S$ in index notation
It is sometimes useful to express {eq}`eq-neo-hookean-incremental-stress` in index notation,

$$
\begin{aligned}
\diff\bm S_{IJ} &= \frac{\partial \bm S_{IJ}}{\partial \bm E_{KL}} \diff \bm E_{KL} \\
  &= \lambda (\bm C^{-1}_{KL} \diff\bm E_{KL}) \bm C^{-1}_{IJ} + 2 (\mu - \lambda \log J) \bm C^{-1}_{IK} \diff\bm E_{KL} \bm C^{-1}_{LJ} \\
  &= \underbrace{\Big( \lambda \bm C^{-1}_{IJ} \bm C^{-1}_{KL} + 2 (\mu - \lambda \log J) \bm C^{-1}_{IK} \bm C^{-1}_{JL} \Big)}_{\mathsf C_{IJKL}} \diff \bm E_{KL} \,,
\end{aligned}
$$ (eq-neo-hookean-incremental-stress-index)

where we have identified the effective elasticity tensor $\mathsf C = \mathsf C_{IJKL}$.
It is generally not desirable to store $\mathsf C$, but rather to use the earlier expressions so that only $3\times 3$ tensors (most of which are symmetric) must be manipulated.
That is, given the linearization point $\bm F$ and solution increment $\diff \bm F = \nabla_X (\diff \bm u)$ (which we are solving for in the Newton step), we compute $\diff \bm P$ via

1. recover $\bm C^{-1}$ and $\log J$ (either stored at quadrature points or recomputed),
2. proceed with $3\times 3$ matrix products as in {eq}`eq-neo-hookean-incremental-stress` or the second line of {eq}`eq-neo-hookean-incremental-stress-index` to compute $\diff \bm S$ while avoiding computation or storage of higher order tensors, and
3. conclude by {eq}`eq-diff-P`, where $\bm S$ is either stored or recomputed from its definition exactly as in the nonlinear residual evaluation.
:::

Note that the Newton linearization of {eq}`hyperelastic-weak-form-initial` may be written as a weak form for linear operators: find $\diff\bm u \in \mathcal V_0$ such that

$$
\int_{\Omega_0} \nabla_X \bm v \!:\! \diff\bm P dV = \text{rhs}, \quad \forall \bm v \in \mathcal V_0,
$$

where $\diff \bm P$ is defined by {eq}`eq-diff-P` and {eq}`eq-neo-hookean-incremental-stress`, and $\mathcal V_0$ is the homogeneous space corresponding to $\mathcal V$.

:::{note}
The decision of whether to recompute or store functions of the current state $\bm F$ depends on a roofline analysis {cite}`williams2009roofline,brown2010` of the computation and the cost of the constitutive model.
For low-order elements where flops tend to be in surplus relative to memory bandwidth, recomputation is likely to be preferable, where as the opposite may be true for high-order elements.
Similarly, analysis with a simple constitutive model may see better performance while storing little or nothing while an expensive model such as Arruda-Boyce {cite}`arruda1993largestretch`, which contains many special functions, may be faster when using more storage to avoid recomputation.
In the case where complete linearization is preferred, note the symmetry $\mathsf C_{IJKL} = \mathsf C_{KLIJ}$ evident in {eq}`eq-neo-hookean-incremental-stress-index`, thus $\mathsf C$ can be stored as a symmetric $6\times 6$ matrix, which has 21 unique entries.
Along with 6 entries for $\bm S$, this totals 27 entries of overhead compared to computing everything from $\bm F$.
This compares with 13 entries of overhead for direct storage of $\{ \bm S, \bm C^{-1}, \log J \}$, which is sufficient for the Neo-Hookean model to avoid all but matrix products.
:::

(problem-hyperelasticity-finite-strain-current-configuration)=

## Hyperelasticity in current configuration

In the preceeding discussion, all equations have been formulated in the initial configuration.
This may feel convenient in that the computational domain is clearly independent of the solution, but there are some advantages to defining the equations in the current configuration.

1. Body forces (like gravity), traction, and contact are more easily defined in the current configuration.
2. Mesh quality in the initial configuration can be very bad for large deformation.
3. The required storage and numerical representation can be smaller in the current configuration.

Most of the benefit in case 3 can be attained solely by moving the Jacobian representation to the current configuration {cite}`davydov2020matrix`, though residual evaluation may also be slightly faster in current configuration.
There are multiple commuting paths from the nonlinear weak form in initial configuration {eq}`hyperelastic-weak-form-initial` to the Jacobian weak form in current configuration {eq}`jacobian-weak-form-current`.
One may push forward to the current configuration and then linearize or linearize in initial configuration and then push forward, as summarized below.

$$
\begin{CD}
  {\overbrace{\nabla_X \bm{v} \tcolon \bm{FS}}^{\text{Initial Residual}}}
  @>{\text{push forward}}>{}>
  {\overbrace{\nabla_x \bm{v} \tcolon \bm{\tau}}^{\text{Current Residual}}} \\
  @V{\text{linearize}}V{\begin{smallmatrix} \diff\bm F = \nabla_X\diff\bm u \\ \diff\bm S(\diff\bm E) \end{smallmatrix}}V
  @V{\begin{smallmatrix} \diff\nabla_x\bm v = -\nabla_x\bm v \nabla_x \diff\bm u \\ \diff\bm\tau(\diff\bm\epsilon) \end{smallmatrix}}V{\text{linearize}}V \\
  {\underbrace{\nabla_X\bm{v}\tcolon \Big(\diff\bm{F}\bm{S} + \bm{F}\diff\bm{S}\Big)}_\text{Initial Jacobian}}
  @>{\text{push forward}}>{}>
  {\underbrace{\nabla_x\bm{v}\tcolon \Big(\diff\bm{\tau} -\bm{\tau}(\nabla_x \diff\bm{u})^T \Big)}_\text{Current Jacobian}}
\end{CD}
$$ (initial-current-linearize)

We will follow both paths for consistency and because both intermediate representations may be useful for implementation.

### Push forward, then linearize

The first term of {eq}`hyperelastic-weak-form-initial` can be rewritten in terms of the symmetric Kirchhoff stress tensor
$\bm{\tau}=J\bm{\sigma}=\bm{P}\bm{F}^T = \bm F \bm S \bm F^T$ as

$$
\nabla_X \bm{v} \tcolon \bm{P} = \nabla_X \bm{v} \tcolon \bm{\tau}\bm{F}^{-T} = \nabla_X \bm{v}\bm{F}^{-1} \tcolon \bm{\tau} = \nabla_x \bm{v} \tcolon \bm{\tau}
$$

therefore, the weak form in terms of $\bm{\tau}$ and $\nabla_x$ with integral over $\Omega_0$ is

$$
\int_{\Omega_0}{\nabla_x \bm{v} \tcolon \bm{\tau}} \, dV
 - \int_{\Omega_0}{\bm{v} \cdot \rho_0 \bm{g}} \, dV
 - \int_{\partial \Omega_0}{\bm{v}\cdot(\bm{P}\cdot\hat{\bm{N}})} \, dS
 = 0, \quad \forall \bm v \in \mathcal V.
$$ (hyperelastic-weak-form-current)

#### Linearize in current configuration

To derive a Newton linearization of {eq}`hyperelastic-weak-form-current`, first we define

$$
\nabla_x \diff \bm{u} = \nabla_X \diff \bm{u} \  \bm{F}^{-1} = \diff \bm{F} \bm{F}^{-1}
$$ (nabla_xdu)

and $\bm{\tau}$ for Neo-Hookean materials as the push forward of {eq}`neo-hookean-stress`

$$
\bm{\tau} = \bm{F}\bm{S}\bm{F}^T = \mu (\bm{b} - \bm I_3) + \lambda \log J \bm{I}_3,
$$ (tau-neo-hookean)

where $\bm{b} = \bm{F} \bm{F}^T$, is the left Cauchy-Green tensor.
Then by expanding the directional derivative of $\nabla_x \bm{v} \tcolon \bm{\tau}$, we arrive at

$$
\diff \ (\nabla_x \bm{v} \tcolon \bm{\tau}) = \diff \ (\nabla_x \bm{v})\tcolon \bm{\tau} + \nabla_x \bm{v} \tcolon \diff \bm{\tau} .
$$ (hyperelastic-linearization-current1)

The first term of {eq}`hyperelastic-linearization-current1` can be written as

$$
\begin{aligned} \diff \ (\nabla_x \bm{v})\tcolon \bm{\tau} &= \diff \ (\nabla_X \bm{v} \bm{F}^{-1})\tcolon \bm{\tau} = \Big(\underbrace{\nabla_X (\diff \bm{v})}_{0}\bm{F}^{-1} +  \nabla_X \bm{v}\diff \bm{F}^{-1}\Big)\tcolon \bm{\tau}\\   &= \Big(-\nabla_X \bm{v} \bm{F}^{-1}\diff\bm{F}\bm{F}^{-1}\Big)\tcolon \bm{\tau}=\Big(-\nabla_x \bm{v} \diff\bm{F}\bm{F}^{-1}\Big)\tcolon \bm{\tau}\\   &= \Big(-\nabla_x \bm{v} \nabla_x \diff\bm{u} \Big)\tcolon \bm{\tau}= -\nabla_x \bm{v}\tcolon\bm{\tau}(\nabla_x \diff\bm{u})^T \,, \end{aligned}
$$

where we have used $\diff \bm{F}^{-1}=-\bm{F}^{-1} \diff \bm{F} \bm{F}^{-1}$ and {eq}`nabla_xdu`.
Using this and {eq}`hyperelastic-linearization-current1` in {eq}`hyperelastic-weak-form-current` yields the weak form in the current configuration

$$
\int_{\Omega_0} \nabla_x \bm v \tcolon \Big(\diff\bm\tau - \bm\tau (\nabla_x \diff\bm u)^T \Big) = \text{rhs}.
$$ (jacobian-weak-form-current)

In the following, we will sometimes make use of the incremental strain tensor in the current configuration,

$$
\diff\bm\epsilon \equiv \frac{1}{2}\Big(\nabla_x \diff\bm{u} + (\nabla_x \diff\bm{u})^T   \Big) .
$$

:::{dropdown} Deriving $\diff\bm\tau$ for Neo-Hookean material
To derive a useful expression of $\diff\bm\tau$ for Neo-Hookean materials, we will use the representations

$$
\begin{aligned}
\diff \bm{b} &= \diff \bm{F} \bm{F}^T + \bm{F} \diff \bm{F}^T \\
&= \nabla_x \diff \bm{u} \ \bm{b} + \bm{b} \ (\nabla_x \diff \bm{u})^T \\
&= (\nabla_x \diff\bm u)(\bm b - \bm I_3) + (\bm b - \bm I_3) (\nabla_x \diff\bm u)^T + 2 \diff\bm\epsilon
\end{aligned}
$$

and

$$
\begin{aligned} \diff\ (\log J) &= \frac{\partial \log J}{\partial \bm{b}}\tcolon \diff \bm{b} = \frac{\partial J}{J\partial \bm{b}}\tcolon \diff \bm{b}=\frac{1}{2}\bm{b}^{-1}\tcolon \diff \bm{b} \\ &= \frac 1 2 \bm b^{-1} \tcolon \Big(\nabla_x \diff\bm u \ \bm b + \bm b (\nabla_x \diff\bm u)^T \Big) \\ &= \trace (\nabla_x \diff\bm u) \\ &= \trace \diff\bm\epsilon . \end{aligned}
$$

Substituting into {eq}`tau-neo-hookean` gives

$$
\begin{aligned}
\diff \bm{\tau} &= \mu \diff \bm{b} + \lambda \trace (\diff\bm\epsilon) \bm I_3 \\
&= \underbrace{2 \mu \diff\bm\epsilon + \lambda \trace (\diff\bm\epsilon) \bm I_3 - 2\lambda \log J \diff\bm\epsilon}_{\bm F \diff\bm S \bm F^T} \\
&\quad + (\nabla_x \diff\bm u)\underbrace{\Big( \mu (\bm b - \bm I_3) + \lambda \log J \bm I_3 \Big)}_{\bm\tau} \\
&\quad + \underbrace{\Big( \mu (\bm b - \bm I_3) + \lambda \log J \bm I_3 \Big)}_{\bm\tau}  (\nabla_x \diff\bm u)^T ,
\end{aligned}
$$ (dtau-neo-hookean)

where the final expression has been identified according to

$$
\diff\bm\tau = \diff\ (\bm F \bm S \bm F^T) = (\nabla_x \diff\bm u) \bm\tau + \bm F \diff\bm S \bm F^T + \bm\tau(\nabla_x \diff\bm u)^T.
$$
:::

Collecting terms, we may thus opt to use either of the two forms

$$
\begin{aligned}
\diff \bm{\tau} -\bm{\tau}(\nabla_x \diff\bm{u})^T &= (\nabla_x \diff\bm u)\bm\tau + \bm F \diff\bm S \bm F^T \\
&= (\nabla_x \diff\bm u)\bm\tau + \lambda \trace(\diff\bm\epsilon) \bm I_3 + 2(\mu - \lambda \log J) \diff\bm\epsilon,
\end{aligned}
$$ (cur_simp_Jac)

with the last line showing the especially compact representation available for Neo-Hookean materials.

### Linearize, then push forward

We can move the derivatives to the current configuration via

$$
\nabla_X \bm v \!:\! \diff\bm P = (\nabla_X \bm v) \bm F^{-1} \!:\! \diff \bm P \bm F^T = \nabla_x \bm v \!:\! \diff\bm P \bm F^T
$$

and expand

$$
\begin{aligned}
\diff\bm P \bm F^T &= \diff\bm F \bm S \bm F^T + \bm F \diff\bm S \bm F^T \\
&= \underbrace{\diff\bm F \bm F^{-1}}_{\nabla_x \diff\bm u} \underbrace{\bm F \bm S \bm F^T}_{\bm\tau} + \bm F \diff\bm S \bm F^T .
\end{aligned}
$$

:::{dropdown} Representation of $\bm F \diff\bm S \bm F^T$ for Neo-Hookean materials
Now we push {eq}`eq-neo-hookean-incremental-stress` forward via

$$
\begin{aligned}
\bm F \diff\bm S \bm F^T &= \lambda (\bm C^{-1} \!:\! \diff\bm E) \bm F \bm C^{-1} \bm F^T
  + 2 (\mu - \lambda \log J) \bm F \bm C^{-1} \diff\bm E \, \bm C^{-1} \bm F^T \\
    &= \lambda (\bm C^{-1} \!:\! \diff\bm E) \bm I_3 + 2 (\mu - \lambda \log J) \bm F^{-T} \diff\bm E \, \bm F^{-1} \\
    &= \lambda \operatorname{trace}(\nabla_x \diff\bm u) \bm I_3 + 2 (\mu - \lambda \log J) \diff\bm \epsilon
\end{aligned}
$$

where we have used

$$
\begin{aligned}
\bm C^{-1} \!:\! \diff\bm E &= \bm F^{-1} \bm F^{-T} \!:\! \bm F^T \diff\bm F \\
&= \operatorname{trace}(\bm F^{-1} \bm F^{-T} \bm F^T \diff \bm F) \\
&= \operatorname{trace}(\bm F^{-1} \diff\bm F) \\
&= \operatorname{trace}(\diff \bm F \bm F^{-1}) \\
&= \operatorname{trace}(\nabla_x \diff\bm u)
\end{aligned}
$$

and

$$
\begin{aligned}
\bm F^{-T} \diff\bm E \, \bm F^{-1} &= \frac 1 2 \bm F^{-T} (\bm F^T \diff\bm F + \diff\bm F^T \bm F) \bm F^{-1} \\
&= \frac 1 2 (\diff \bm F \bm F^{-1} + \bm F^{-T} \diff\bm F^T) \\
&= \frac 1 2 \Big(\nabla_x \diff\bm u + (\nabla_x\diff\bm u)^T \Big) \equiv \diff\bm\epsilon.
\end{aligned}
$$
:::

Collecting terms, the weak form of the Newton linearization for Neo-Hookean materials in the current configuration is

$$
\int_{\Omega_0} \nabla_x \bm v \!:\! \Big( (\nabla_x \diff\bm u) \bm\tau + \lambda \operatorname{trace}(\diff\bm\epsilon)\bm I_3 + 2(\mu - \lambda\log J)\diff \bm\epsilon \Big) dV = \text{rhs},
$$ (jacobian-weak-form-current2)

which equivalent to Algorithm 2 of {cite}`davydov2020matrix` and requires only derivatives with respect to the current configuration. Note that {eq}`cur_simp_Jac` and {eq}`jacobian-weak-form-current2` have recovered the same representation
using different algebraic manipulations.

:::{tip}
We define a second order *Green-Euler* strain tensor (cf. Green-Lagrange strain {eq}`eq-green-lagrange-strain`) as

$$
\bm e = \frac 1 2 \Big(\bm{b} - \bm{I}_3 \Big) = \frac 1 2 \Big( \nabla_X \bm{u} + (\nabla_X \bm{u})^T + \nabla_X \bm{u} \, (\nabla_X \bm{u})^T \Big).
$$ (green-euler-strain)

Then, the Kirchhoff stress tensor {eq}`tau-neo-hookean` can be written as

$$
\bm \tau = \lambda \log J \bm I_{3} + 2\mu \bm e,
$$ (tau-neo-hookean-stable)

which is more numerically stable for small strain, and thus preferred for computation. Note that the $\log J$ is computed via `log1p` {eq}`log1p`, as we discussed in the previous tip.
:::

### Jacobian representation

We have implemented four storage variants for the Jacobian in our finite strain hyperelasticity. In each case, some variables are computed during residual evaluation and used during Jacobian application.

:::{list-table} Four algorithms for Jacobian action in finite strain hyperelasticity problem
:header-rows: 1
:widths: auto

* - Option `-problem`
  - Static storage
  - Computed storage
  - \# scalars
  - Equations


* - `FSInitial-NH1`
  - $\nabla_{X} \hat X, \operatorname{det}\nabla_{\hat X} X$
  - $\nabla_X \bm u$
  - 19
  - {eq}`eq-diff-P` {eq}`eq-neo-hookean-incremental-stress`

* - `FSInitial-NH2`
  - $\nabla_{X} \hat X, \operatorname{det}\nabla_{\hat X} X$
  - $\nabla_X \bm u, \bm C^{-1}, \lambda \log J$
  - 26
  - {eq}`eq-diff-P` {eq}`eq-neo-hookean-incremental-stress`

* - `FSCurrent-NH1`
  - $\nabla_{X} \hat X, \operatorname{det}\nabla_{\hat X} X$
  - $\nabla_X \bm u$
  - 19
  - {eq}`jacobian-weak-form-current` {eq}`nabla_xdu`

* - `FSCurrent-NH2`
  - $\operatorname{det}\nabla_{\hat X} X$
  - $\nabla_x \hat X, \bm \tau, \lambda \log J$
  - 17
  - {eq}`jacobian-weak-form-current` {eq}`jacobian-weak-form-current2`
:::
# libCEED: Solid Mechanics Example

This page provides a description of the solid mechanics example for the
libCEED library, based on PETSc.

This code solves the steady-state static momentum balance equations using unstructured high-order finite/spectral element spatial discretizations.
In this mini-app, we consider three formulations used in solid mechanics applications: linear elasticity, Neo-Hookean hyperelasticity at small strain, and Neo-Hookean hyperelasticity at finite strain.
All three of these formulations are for compressible materials.

Build by using:

```
make
```

and run with:

```
./elasticity -mesh [.exo file] -degree [degree] -nu [nu] -E [E] [boundary options] -problem [problem type] -forcing [forcing] -ceed [ceed]
```

## Runtime options

% inclusion-solids-marker

The elasticity mini-app is controlled via command-line options, the following of which are mandatory.

:::{list-table} Mandatory Runtime Options
:header-rows: 1
:widths: 3 7

* - Option
  - Description
* - `-mesh [filename]`
  - Path to mesh file in any format supported by PETSc.
* - `-degree [int]`
  - Polynomial degree of the finite element basis
* - `-E [real]`
  - [Young's modulus](https://en.wikipedia.org/wiki/Young%27s_modulus), $E > 0$
* - `-nu [real]`
  - [Poisson's ratio](https://en.wikipedia.org/wiki/Poisson%27s_ratio), $\nu < 0.5$
* - `-bc_clamp [int list]`
  - List of face sets on which to displace by `-bc_clamp_[facenumber]_translate [x,y,z]`
    and/or `bc_clamp_[facenumber]_rotate [rx,ry,rz,c_0,c_1]`. Note: The default
    for a clamped face is zero displacement. All displacement is with respect to
    the initial configuration.
* - `-bc_traction [int list]`
  - List of face sets on which to set traction boundary conditions with the
    traction vector `-bc_traction_[facenumber] [tx,ty,tz]`
:::

:::{note}
This solver can use any mesh format that PETSc's `DMPlex` can read (Exodus, Gmsh, Med, etc.).
Our tests have primarily been using Exodus meshes created using [CUBIT]; sample meshes used for the example runs suggested here can be found in [this repository].
Note that many mesh formats require PETSc to be configured appropriately; e.g., `--download-exodusii` for Exodus support.
:::

Consider the specific example of the mesh seen below:

```{image} https://github.com/jeremylt/ceedSampleMeshes/raw/master/cylinderDiagram.png
```

With the sidesets defined in the figure, we provide here an example of a minimal set of command line options:

```
./elasticity -mesh [.exo file] -degree 4 -E 1e6 -nu 0.3 -bc_clamp 998,999 -bc_clamp_998_translate 0,-0.5,1
```

In this example, we set the left boundary, face set $999$, to zero displacement and the right boundary, face set $998$, to displace $0$ in the $x$ direction, $-0.5$ in the $y$, and $1$ in the $z$.

As an alternative to specifying a mesh with {code}`-mesh`, the user may use a DMPlex box mesh by specifying {code}`-dm_plex_box_faces [int list]`, {code}`-dm_plex_box_upper [real list]`, and {code}`-dm_plex_box_lower [real list]`.

As an alternative example exploiting {code}`-dm_plex_box_faces`, we consider a {code}`4 x 4 x 4` mesh where essential (Drichlet) boundary condition is placed on all sides. Sides 1 through 6 are rotated around $x$-axis:

```
./elasticity -problem FSInitial-NH1 -E 1 -nu 0.3 -num_steps 40 -snes_linesearch_type cp -dm_plex_box_faces 4,4,4 -bc_clamp 1,2,3,4,5,6 -bc_clamp_1_rotate 0,0,1,0,.3 -bc_clamp_2_rotate 0,0,1,0,.3 -bc_clamp_3_rotate 0,0,1,0,.3 -bc_clamp_4_rotate 0,0,1,0,.3 -bc_clamp_5_rotate 0,0,1,0,.3 -bc_clamp_6_rotate 0,0,1,0,.3
```

:::{note}
If the coordinates for a particular side of a mesh are zero along the axis of rotation, it may appear that particular side is clamped zero.
:::

On each boundary node, the rotation magnitude is computed: {code}`theta = (c_0 + c_1 * cx) * loadIncrement` where {code}`cx = kx * x + ky * y + kz * z`, with {code}`kx`, {code}`ky`, {code}`kz` are normalized values.

The command line options just shown are the minimum requirements to run the mini-app, but additional options may also be set as follows

:::{list-table} Additional Runtime Options
:header-rows: 1

* - Option
  - Description
  - Default value

* - `-ceed`
  - CEED resource specifier
  - `/cpu/self`

* - `-qextra`
  - Number of extra quadrature points
  - `0`

* - `-test`
  - Run in test mode
  -

* - `-problem`
  - Problem to solve (`Linear`, `SS-NH`, `FSInitial-NH1`, etc.)
  - `Linear`

* - `-forcing`
  -  Forcing term option (`none`, `constant`, or `mms`)
  - `none`

* - `-forcing_vec`
  -  Forcing vector
  - `0,-1,0`

* - `-multigrid`
  - Multigrid coarsening to use (`logarithmic`, `uniform` or `none`)
  - `logarithmic`

* - `-nu_smoother [real]`
  - Poisson's ratio for multigrid smoothers, $\nu < 0.5$
  -

* - `-num_steps`
  - Number of load increments for continuation method
  - `1` if `Linear` else `10`

* - `-view_soln`
  - Output solution at each load increment for viewing
  -

* - `-view_final_soln`
  - Output solution at final load increment for viewing
  -

* - `-snes_view`
  - View PETSc `SNES` nonlinear solver configuration
  -

* - `-log_view`
  - View PETSc performance log
  -

* - `-output_dir`
  - Output directory
  - `.`

* - `-help`
  - View comprehensive information about run-time options
  -
:::

To verify the convergence of the linear elasticity formulation on a given mesh with the method of manufactured solutions, run:

```
./elasticity -mesh [mesh] -degree [degree] -nu [nu] -E [E] -forcing mms
```

This option attempts to recover a known solution from an analytically computed forcing term.

### On algebraic solvers

This mini-app is configured to use the following Newton-Krylov-Multigrid method by default.

- Newton-type methods for the nonlinear solve, with the hyperelasticity models globalized using load increments.
- Preconditioned conjugate gradients to solve the symmetric positive definite linear systems arising at each Newton step.
- Preconditioning via $p$-version multigrid coarsening to linear elements, with algebraic multigrid (PETSc's `GAMG`) for the coarse solve.
  The default smoother uses degree 3 Chebyshev with Jacobi preconditioning.
  (Lower degree is often faster, albeit less robust; try {code}`-outer_mg_levels_ksp_max_it 2`, for example.)
  Application of the linear operators for all levels with degree $p > 1$ is performed matrix-free using analytic Newton linearization, while the lowest order $p = 1$ operators are assembled explicitly (using coloring at present).

Many related solvers can be implemented by composing PETSc command-line options.

### Nondimensionalization

Quantities such as the Young's modulus vary over many orders of magnitude, and thus can lead to poorly scaled equations.
One can nondimensionalize the model by choosing an alternate system of units, such that displacements and residuals are of reasonable scales.

:::{list-table} (Non)dimensionalization options
:header-rows: 1

* - Option
  - Description
  - Default value

* - :code:`-units_meter`
  - 1 meter in scaled length units
  - :code:`1`

* - :code:`-units_second`
  - 1 second in scaled time units
  - :code:`1`

* - :code:`-units_kilogram`
  - 1 kilogram in scaled mass units
  - :code:`1`
:::

For example, consider a problem involving metals subject to gravity.

:::{list-table} Characteristic units for metals
:header-rows: 1

* - Quantity
  - Typical value in SI units

* - Displacement, $\bm u$
  - $1 \,\mathrm{cm} = 10^{-2} \,\mathrm m$

* - Young's modulus, $E$
  - $10^{11} \,\mathrm{Pa} = 10^{11} \,\mathrm{kg}\, \mathrm{m}^{-1}\, \mathrm s^{-2}$

* - Body force (gravity) on volume, $\int \rho \bm g$
  - $5 \cdot 10^4 \,\mathrm{kg}\, \mathrm m^{-2} \, \mathrm s^{-2} \cdot (\text{volume} \, \mathrm m^3)$
:::

One can choose units of displacement independently (e.g., {code}`-units_meter 100` to measure displacement in centimeters), but $E$ and $\int \rho \bm g$ have the same dependence on mass and time, so cannot both be made of order 1.
This reflects the fact that both quantities are not equally significant for a given displacement size; the relative significance of gravity increases as the domain size grows.

### Diagnostic Quantities

Diagnostic quantities for viewing are provided when the command line options for visualization output, {code}`-view_soln` or {code}`-view_final_soln` are used.
The diagnostic quantities include displacement in the $x$ direction, displacement in the $y$ direction, displacement in the $z$ direction, pressure, $\operatorname{trace} \bm{E}$, $\operatorname{trace} \bm{E}^2$, $\lvert J \rvert$, and strain energy density.
The table below summarizes the formulations of each of these quantities for each problem type.

:::{list-table} Diagnostic quantities
   :header-rows: 1

   * - Quantity
     - Linear Elasticity
     - Hyperelasticity, Small Strain
     - Hyperelasticity, Finite Strain

   * - Pressure
     - $\lambda \operatorname{trace} \bm{\epsilon}$
     - $\lambda \log \operatorname{trace} \bm{\epsilon}$
     - $\lambda \log J$

   * - Volumetric Strain
     - $\operatorname{trace} \bm{\epsilon}$
     - $\operatorname{trace} \bm{\epsilon}$
     - $\operatorname{trace} \bm{E}$

   * - $\operatorname{trace} \bm{E}^2$
     - $\operatorname{trace} \bm{\epsilon}^2$
     - $\operatorname{trace} \bm{\epsilon}^2$
     - $\operatorname{trace} \bm{E}^2$

   * - $\lvert J \rvert$
     - $1 + \operatorname{trace} \bm{\epsilon}$
     - $1 + \operatorname{trace} \bm{\epsilon}$
     - $\lvert J \rvert$

   * - Strain Energy Density
     - $\frac{\lambda}{2} (\operatorname{trace} \bm{\epsilon})^2 + \mu \bm{\epsilon} : \bm{\epsilon}$
     - $\lambda (1 + \operatorname{trace} \bm{\epsilon}) (\log(1 + \operatorname{trace} \bm{\epsilon} ) - 1) + \mu \bm{\epsilon} : \bm{\epsilon}$
     - $\frac{\lambda}{2}(\log J)^2 + \mu \operatorname{trace} \bm{E} - \mu \log J$
:::

[cubit]: https://cubit.sandia.gov/
[this repository]: https://github.com/jeremylt/ceedSampleMeshes
# libCEED: Documentation

This page provides a brief description of the documentation for the libCEED library.

## Quick build

If you have Python and Doxygen installed, these two commands should build the documentation in `doc/sphinx/build/html/`.

```sh
pip install --user -r doc/sphinx/requirements.txt  # only needed once
make doc                                           # builds HTML site
```

## Sphinx

Sphinx is the tool used for libCEED's User Manual. Sphinx can produce documentation in different output formats: HTML, LaTeX (for printable PDF versions), ePub, Texinfo, manual pages, and plain text. Sphinx comes with a broad set of extensions for different features, for instance the automatic inclusion of documentation from docstrings and snippets of codes, support of todo items, highlighting of code, and math rendering.

To be able to contribute to libCEED's User Manual, Sphinx needs to be [installed](http://www.sphinx-doc.org/en/master/usage/installation.html) together with its desired extensions.

The Sphinx API documentation depends on Doxygen's XML output (via the `breathe` plugin).  Build these files in the `xml/` directory via:

```sh
make doxygen
```

If you are editing documentation, such as the reStructuredText files in `doc/sphinx/source`, you can rebuild incrementally via

```sh
make -C doc/sphinx html
```

which will HTML docs in the [doc/sphinx/build](./sphinx/build) directory. Use

```sh
make -C doc/sphinx latexpdf
```

to build PDF using the LaTeX toolchain (which must be installed).
This requires the `rsvg-convert` utility, which is likely available from your package manager under `librsvg` or `librsvg2-bin`.

For more Sphinx features, see

```sh
make -C doc/sphinx help
```

### Dependencies

Some of the extensions used require installation. They are distributed on [PyPI](https://pypi.org) and can be installed with `pip`. The extensions used for this project can be found in the [requiremenets file](./sphinx/requirements.txt) and can be readily installed by running:

```sh
pip install --user -r doc/sphinx/requirements.txt
```

from toplevel.
You can use `virtualenv` to keep these extensions isolated, and to check that all necessary extensions are included.

```sh
virtualenv VENV                              # create a virtual environment
. VENV/bin/active                            # activate the environment
pip install -r doc/sphinx/requirements.txt   # install dependencies inside VENV
make doc
```
---
title: 'libCEED: Fast algebra for high-order element-based discretizations'
tags:
  - high-performance computing
  - high-order methods
  - finite elements
  - spectral elements
  - matrix-free
authors:
  - name: Jed Brown
    orcid: 0000-0002-9945-0639
    affiliation: 1
  - name: Ahmad Abdelfattah
    orcid: 0000-0001-5054-4784
    affiliation: 3
  - name: Valeria Barra
    orcid: 0000-0003-1129-2056
    affiliation: 1
  - name: Natalie Beams
    orcid: 0000-0001-6060-4082
    affiliation: 3
  - name: Jean-Sylvain Camier
    orcid: 0000-0003-2421-1999
    affiliation: 2
  - name: Veselin Dobrev
    orcid: 0000-0003-1793-5622
    affiliation: 2
  - name: Yohann Dudouit
    orcid: 0000-0001-5831-561X
    affiliation: 2
  - name:  Leila Ghaffari
    orcid: 0000-0002-0965-214X
    affiliation: 1
  - name: Tzanio Kolev
    orcid: 0000-0002-2810-3090
    affiliation: 2
  - name: David Medina
    affiliation: 4
  - name: Will Pazner
    orcid: 0000-0003-4885-2934
    affiliation: 2
  - name: Thilina Ratnayaka
    orcid: 0000-0001-6102-6560
    affiliation: 5
  - name: Jeremy Thompson
    orcid: 0000-0003-2980-0899
    affiliation: 1
  - name: Stan Tomov
    orcid: 0000-0002-5937-7959
    affiliation: 3
affiliations:
 - name: University of Colorado at Boulder
   index: 1
 - name: Lawrence Livermore National Laboratory
   index: 2
 - name: University of Tennessee
   index: 3
 - name: Occalytics LLC
   index: 4
 - name: University of Illinois at Urbana-Champaign
   index: 5
date: 9 July 2021
bibliography: paper.bib
---

# Summary and statement of need

Finite element methods are widely used to solve partial differential equations (PDE) in science and engineering, but their standard implementation [@dealII92;@libMeshPaper;@LoggMardalWells2012] relies on assembling sparse matrices.
Sparse matrix multiplication and triangular operations perform a scalar multiply and add for each nonzero entry, just 2 floating point operations (flops) per scalar that must be loaded from memory [@williams2009roofline].
Modern hardware is capable of nearly 100 flops per scalar streamed from memory [@kruppcomparison] so sparse matrix operations cannot achieve more than about 2% utilization of arithmetic units.
Matrix assembly becomes even more problematic when the polynomial degree $p$ of the basis functions is increased, resulting in $O(p^d)$ storage and $O(p^{2d})$ compute per degree of freedom (DoF) in $d$ dimensions.
Methods pioneered by the spectral element community [@Orszag:1980; @deville2002highorder] exploit problem structure to reduce costs to $O(1)$ storage and $O(p)$ compute per DoF, with very high utilization of modern CPUs and GPUs.
Unfortunately, high-quality implementations have been relegated to applications and intrusive frameworks that are often difficult to extend to new problems or incorporate into legacy applications, especially when strong preconditioners are required.

`libCEED`, the Code for Efficient Extensible Discretization [@libceed-user-manual], is a lightweight library that provides a purely algebraic interface for linear and nonlinear operators and preconditioners with element-based discretizations.
`libCEED` provides portable performance via run-time selection of implementations optimized for CPUs and GPUs, including support for just-in-time (JIT) compilation.
It is designed for convenient use in new and legacy software, and offers interfaces in C99 [@C99-lang], Fortran77 [@Fortran77-lang], Python [@Python-lang], Julia [@Julia-lang], and Rust [@Rust-lang].
Users and library developers can integrate `libCEED` at a low level into existing applications in place of existing matrix-vector products without significant refactoring of their own discretization infrastructure.
Alternatively, users can utilize integrated `libCEED` support in MFEM [@MFEMlibrary; @mfem-paper].

In addition to supporting applications and discretization libraries, `libCEED` provides a platform for performance engineering and co-design, as well as an algebraic interface for solvers research like adaptive $p$-multigrid, much like how sparse matrix libraries enable development and deployment of algebraic multigrid solvers.

# Concepts and interface

Consider finite element discretization of a problem based on a weak form with one weak derivative: find $u$ such that

$$ v^T F(u) := \int_\Omega v \cdot f_0(u, \nabla u) + \nabla v \!:\! f_1(u, \nabla u) = 0 \quad \forall v, $$

where the functions $f_0$ and $f_1$ define the physics and possible stabilization of the problem [@Brown:2010] and the functions $u$ and $v$ live in a suitable space.
Integrals in the weak form are evaluated by summing over elements $e$,

$$ F(u) = \sum_e \mathcal E_e^T B_e^T W_e f(B_e \mathcal E_e u), $$

where $\mathcal E_e$ restricts to element $e$, $B_e$ evaluates solution values and derivatives to quadrature points, $f$ acts independently at quadrature points, and $W_e$ is a (diagonal) weighting at quadrature points.
By grouping the operations $W_e$ and $f$ into a point-block diagonal $D$ and stacking the restrictions $\mathcal E_e$ and basis actions $B_e$ for each element, we can express the global residual in operator notation (\autoref{fig:decomposition}), where $\mathcal P$ is an optional external operator, such as the parallel restriction in MPI-based [@gropp2014using] solvers.
Inhomogeneous Neumann, Robin, and nonlinear boundary conditions can be added in a similar fashion by adding terms integrated over boundary faces while Dirichlet boundary conditions can be added by setting the target values prior to applying the operator representing the weak form.
Similar face integral terms can also be used to represent discontinuous Galerkin formulations.

![`libCEED` uses a logical decomposition to define element-based discretizations, with optimized implementations of the action and preconditioning ingredients. \label{fig:decomposition}](img/libCEED-2-trim.pdf)

`libCEED`'s native C interface is object-oriented, providing data types for each logical object in the decomposition.

Symbol        `libCEED` type             Description
------        ------------             -----------
$D$           `CeedQFunction`          User-defined action at quadrature points
$B$           `CeedBasis`              Basis evaluation to quadrature (dense/structured)
$\mathcal E$  `CeedElemRestriction`    Restriction to each element (sparse/boolean)
$A$           `CeedOperator`           Linear or nonlinear operator acting on L-vectors

`libCEED` implementations ("backends") are free to reorder and fuse computational steps (including eliding memory to store intermediate representations) so long as the mathematical properties of the operator $A$ are preserved.
A `CeedOperator` is composed of one or more operators defined as in \autoref{fig:decomposition}, and acts on a `CeedVector`, which typically encapsulates zero-copy access to host or device memory provided by the caller.
The element restriction $\mathcal E$ requires mesh topology and a numbering of DoFs, and may be a no-op when data is already composed by element (such as with discontinuous Galerkin methods).
The discrete basis $B$ is the purely algebraic expression of a finite element basis (shape functions) and quadrature; it often possesses structure that is exploited to speed up its action.
Some constructors are provided for arbitrary polynomial degree $H^1$ Lagrange bases with a tensor-product representation due to the computational efficiency of computing solution values and derivatives at quadrature points via tensor contractions.
However, the user can define a `CeedBasis` for arbitrary element topology including tetrahedra, prisms, and other realizations of abstract polytopes, by providing quadrature weights and the matrices used to compute solution values and derivatives at quadrature points from the DoFs on the element.

The physics (weak form) is expressed through `CeedQFunction`, which can either be defined by the user or selected from a gallery distributed with `libCEED`.
These pointwise functions do not depend on element resolution, topology, or basis degree (see \autoref{fig:schematic}), in contrast to systems like FEniCS where UFL forms specify basis degree at compile time.
This isolation is valuable for $hp$-refinement and adaptivity (where $h$ commonly denotes the average element size and $p$ the polynomial degree of the basis functions; see @babuska1994hpfem) and $p$-multigrid solvers; mixed-degree, mixed-topology, and $h$-nonconforming finite element methods are readily expressed by composition.
Additionally, a single source implementation (in vanilla C or C++) for the `CeedQFunction`s can be used on CPUs or GPUs (transparently using the @NVRTCwebsite, HIPRTC, or OCCA [@OCCAwebsite] run-time compilation features).

`libCEED` provides computation of the true operator diagonal for preconditioning with Jacobi and Chebyshev as well as direct assembly of sparse matrices (e.g., for coarse operators in multigrid) and construction of $p$-multigrid prolongation and restriction operators.
Preconditioning matrix-free operators is an active area of research; support for domain decomposition methods and inexact subdomain solvers based on the fast diagonalization method [@lottes2005hms] are in active development.

![A schematic of element restriction and basis applicator operators for
elements with different topology. This sketch shows the independence of Q-functions
(in this case representing a Laplacian) on element resolution, topology, and basis degree.\label{fig:schematic}](img/QFunctionSketch.pdf)

# High-level languages

`libCEED` provides high-level interfaces in Python, Julia, and Rust, each of which is maintained and tested as part of the main repository, but distributed through each language's respective package manager.

The Python interface uses CFFI, the C Foreign Function Interface [@python-cffi]. CFFI allows reuse of most C declarations and requires only a minimal adaptation of some of them. The C and Python APIs are mapped in a nearly 1:1 correspondence. For instance, a `CeedVector` object is exposed as `libceed.Vector` in Python, and supports no-copy host and GPU device interperability with Python arrays from the NumPy [@NumPy] or Numba [@Numba] packages. The interested reader can find more details on `libCEED`'s Python interface in @libceed-paper-proc-scipy-2020.

The Julia interface, referred to as `LibCEED.jl`, provides both a low-level interface, which is generated automatically from `libCEED`'s C header files, and a high-level interface. The high-level interface takes advantage of Julia's metaprogramming and just-in-time compilation capabilities to enable concise definition of Q-functions that work on both CPUs and GPUs, along with their composition into operators as in \autoref{fig:decomposition}.

The Rust interface also wraps automatically-generated bindings from the `libCEED` C header files, offering increased safety due to Rust ownership and borrow checking, and more convenient definition of Q-functions (e.g., via closures).

# Backends

\autoref{fig:libCEEDBackends} shows a subset of the backend implementations (backends) available in `libCEED`.
GPU implementations are available via pure @CUDAwebsite and pure @HIPwebsite, as well as the OCCA [@OCCAwebsite] and MAGMA [@MAGMAwebsite] libraries. CPU implementations are available via pure C and AVX intrinsics as well as the LIBXSMM library [@LIBXSMM]. `libCEED` provides a dynamic interface such that users only need to write a single source (no need for templates/generics) and can select the desired specialized implementation at run time. Moreover, each process or thread can instantiate an arbitrary number of backends on an arbitrary number of devices.

![`libCEED` provides the algebraic core for element-based discretizations, with specialized implementations
(backends) for heterogeneous architectures.\label{fig:libCEEDBackends}](img/libCEEDBackends.png)

# Performance benchmarks

The Exascale Computing Project (ECP) co-design Center for Efficient Exascale Discretization [@CEEDwebsite] has defined a suite of Benchmark Problems (BPs) to test and compare the performance of high-order finite element implementations [@Fischer2020scalability; @CEED-ECP-paper]. \autoref{fig:bp3} compares the performance of `libCEED` solving BP3 (CG iteration on a 3D Poisson problem) or CPU and GPU systems of similar (purchase/operating and energy) cost. These tests use PETSc [@PETScUserManual] for unstructured mesh management and parallel solvers with GPU-aware communication [@zhang2021petscsf]; a similar implementation with comparable performance is available through MFEM.

![Performance for BP3 using the \texttt{xsmm/blocked} backend on a 2-socket AMD EPYC 7452 (32-core, 2.35GHz) and the \texttt{cuda/gen} backend on LLNL's Lassen system with NVIDIA V100 GPUs. Each curve represents fixing the basis degree $p$ and varying the number of elements. The CPU enables faster solution of smaller problem sizes (as in strong scaling) while the GPU is more efficient for applications that can afford to wait for larger sizes. Note that the CPU exhibits a performance drop when the working set becomes too large for L3 cache (128 MB/socket) while no such drop exists for the GPU. (This experiment was run with release candidates of PETSc 3.14 and libCEED 0.7 using gcc-10 on EPYC and clang-10/CUDA-10 on Lassen.) \label{fig:bp3}](img/bp3-2020.pdf)

# Demo applications and integration

To highlight the ease of library reuse for solver composition and leverage `libCEED`'s full capability for real-world applications, `libCEED` comes with a suite of application examples, including problems of interest to the fluid dynamics and solid mechanics communities.
The fluid dynamics example solves the 2D and 3D compressible Navier-Stokes equations using SU/SUPG stabilization and implicit, explicit, or IMEX time integration; \autoref{fig:NSvortices} shows vortices arising in the "density current" [@straka1993numerical] when a cold bubble of air reaches the ground.
The solid mechanics example solves static linear elasticity and hyperelasticity with load continuation and Newton-Krylov solvers with $p$-multigrid preconditioners; \autoref{fig:Solids} shows a twisted Neo-Hookean beam. Both of these examples have been developed using PETSc, where `libCEED` provides the matrix-free operator and preconditioner ingredient evaluation and PETSc provides the unstructured mesh management and parallel solvers.

![Vortices develop as a cold air bubble drops to the ground.\label{fig:NSvortices}](img/Vortices.png)

![Strain energy density in a twisted Neo-Hookean beam.\label{fig:Solids}](img/SolidTwistExample.jpeg)

`libCEED` also includes additional examples with PETSc, MFEM, and Nek5000 [@Nekwebsite].

If MFEM is built with `libCEED` support, existing MFEM users can pass `-d ceed-cuda:/gpu/cuda/gen` to use a `libCEED` CUDA backend, and similarly for other backends.
The `libCEED` implementations, accessed in this way, currently provide MFEM users with the fastest operator action on CPUs and GPUs (CUDA and HIP/ROCm) without writing any `libCEED` Q-functions.

# Acknowledgements

This research is supported by the Exascale Computing Project (17-SC-20-SC), a collaborative effort of two U.S. Department of Energy organizations (Office of Science and the National Nuclear Security Administration) responsible for the planning and preparation of a capable exascale ecosystem, including software, applications, hardware, advanced system engineering and early testbed platforms, in support of the nations exascale computing imperative. We thank Lawrence Livermore National Laboratory for access to the Lassen and Corona machines.

# References
# Welcome to libCEED's User Manual!

```{toctree}
:caption: Contents
:maxdepth: 4

intro
gettingstarted
libCEEDapi
examples/index
ffi
api/index
precision
libCEEDdev
Contributing <CONTRIBUTING>
Code of Conduct <CODE_OF_CONDUCT>
releasenotes
```

# Indices and tables

- {ref}`genindex`
- {ref}`search`

```{bibliography}
```
# Interface Concepts

This page provides a brief description of the theoretical foundations and the
practical implementation of the libCEED library.

(theoretical-framework)=

## Theoretical Framework

In finite element formulations, the weak form of a Partial Differential Equation
(PDE) is evaluated on a subdomain $\Omega_e$ (element) and the local results
are composed into a larger system of equations that models the entire problem on
the global domain $\Omega$. In particular, when high-order finite elements or
spectral elements are used, the resulting sparse matrix representation of the global
operator is computationally expensive, with respect to both the memory transfer and
floating point operations needed for its evaluation. libCEED provides an interface
for matrix-free operator description that enables efficient evaluation on a variety
of computational device types (selectable at run time). We present here the notation
and the mathematical formulation adopted in libCEED.

We start by considering the discrete residual $F(u)=0$ formulation
in weak form. We first define the $L^2$ inner product between real-valued functions

$$
\langle v, u \rangle = \int_\Omega v u d \bm{x},
$$

where $\bm{x} \in \mathbb{R}^d \supset \Omega$.

We want to find $u$ in a suitable space $V_D$,
such that

$$
\langle  \bm v,  \bm f(u) \rangle = \int_\Omega  \bm v \cdot  \bm f_0 (u, \nabla u) + \nabla \bm v :  \bm f_1 (u, \nabla u) = 0
$$ (residual)

for all $\bm v$ in the corresponding homogeneous space $V_0$, where $\bm f_0$
and $\bm f_1$ contain all possible sources in the problem. We notice here that
$\bm f_0$ represents all terms in {eq}`residual` which multiply the (possibly vector-valued) test
function $\bm v$ and $\bm f_1$ all terms which multiply its gradient $\nabla \bm v$.
For an n-component problems in $d$ dimensions, $\bm f_0 \in \mathbb{R}^n$ and
$\bm f_1 \in \mathbb{R}^{nd}$.

:::{note}
The notation $\nabla \bm v \!:\! \bm f_1$ represents contraction over both
fields and spatial dimensions while a single dot represents contraction in just one,
which should be clear from context, e.g., $\bm v \cdot \bm f_0$ contracts only over
fields.
:::

:::{note}
In the code, the function that represents the weak form at quadrature
points is called the {ref}`CeedQFunction`. In the {ref}`Examples` provided with the
library (in the {file}`examples/` directory), we store the term $\bm f_0$ directly
into `v`, and the term $\bm f_1$ directly into `dv` (which stands for
$\nabla \bm v$). If equation {eq}`residual` only presents a term of the
type $\bm f_0$, the {ref}`CeedQFunction` will only have one output argument,
namely `v`. If equation {eq}`residual` also presents a term of the type
$\bm f_1$, then the {ref}`CeedQFunction` will have two output arguments, namely,
`v` and `dv`.
:::

## Finite Element Operator Decomposition

Finite element operators are typically defined through weak formulations of
partial differential equations that involve integration over a computational
mesh. The required integrals are computed by splitting them as a sum over the
mesh elements, mapping each element to a simple *reference* element (e.g. the
unit square) and applying a quadrature rule in reference space.

This sequence of operations highlights an inherent hierarchical structure
present in all finite element operators where the evaluation starts on *global
(trial) degrees of freedom (dofs) or nodes on the whole mesh*, restricts to
*dofs on subdomains* (groups of elements), then moves to independent
*dofs on each element*, transitions to independent *quadrature points* in
reference space, performs the integration, and then goes back in reverse order
to global (test) degrees of freedom on the whole mesh.

This is illustrated below for the simple case of symmetric linear operator on
third order ($Q_3$) scalar continuous ($H^1$) elements, where we use
the notions **T-vector**, **L-vector**, **E-vector** and **Q-vector** to represent
the sets corresponding to the (true) degrees of freedom on the global mesh, the split
local degrees of freedom on the subdomains, the split degrees of freedom on the
mesh elements, and the values at quadrature points, respectively.

We refer to the operators that connect the different types of vectors as:

- Subdomain restriction $\bm{P}$
- Element restriction $\bm{G}$
- Basis (Dofs-to-Qpts) evaluator $\bm{B}$
- Operator at quadrature points $\bm{D}$

More generally, when the test and trial space differ, they get their own
versions of $\bm{P}$, $\bm{G}$ and $\bm{B}$.

(fig-operator-decomp)=

:::{figure} ../../img/libCEED.png
Operator Decomposition
:::

Note that in the case of adaptive mesh refinement (AMR), the restrictions
$\bm{P}$ and $\bm{G}$ will involve not just extracting sub-vectors,
but evaluating values at constrained degrees of freedom through the AMR interpolation.
There can also be several levels of subdomains ($\bm P_1$, $\bm P_2$,
etc.), and it may be convenient to split $\bm{D}$ as the product of several
operators ($\bm D_1$, $\bm D_2$, etc.).

### Terminology and Notation

Vector representation/storage categories:

- True degrees of freedom/unknowns, **T-vector**:

  > - each unknown $i$ has exactly one copy, on exactly one processor, $rank(i)$
  > - this is a non-overlapping vector decomposition
  > - usually includes any essential (fixed) dofs.
  >
  > ```{image} ../../img/T-vector.svg
  > ```

- Local (w.r.t. processors) degrees of freedom/unknowns, **L-vector**:

  > - each unknown $i$ has exactly one copy on each processor that owns an
  >   element containing $i$
  > - this is an overlapping vector decomposition with overlaps only across
  >   different processors---there is no duplication of unknowns on a single
  >   processor
  > - the shared dofs/unknowns are the overlapping dofs, i.e. the ones that have
  >   more than one copy, on different processors.
  >
  > ```{image} ../../img/L-vector.svg
  > ```

- Per element decomposition, **E-vector**:

  > - each unknown $i$ has as many copies as the number of elements that contain
  >   $i$
  > - usually, the copies of the unknowns are grouped by the element they belong
  >   to.
  >
  > ```{image} ../../img/E-vector.svg
  > ```

- In the case of AMR with hanging nodes (giving rise to hanging dofs):

  > - the **L-vector** is enhanced with the hanging/dependent dofs
  > - the additional hanging/dependent dofs are duplicated when they are shared
  >   by multiple processors
  > - this way, an **E-vector** can be derived from an **L-vector** without any
  >   communications and without additional computations to derive the dependent
  >   dofs
  > - in other words, an entry in an **E-vector** is obtained by copying an entry
  >   from the corresponding **L-vector**, optionally switching the sign of the
  >   entry (for $H(\mathrm{div})$---and $H(\mathrm{curl})$-conforming spaces).
  >
  > ```{image} ../../img/L-vector-AMR.svg
  > ```

- In the case of variable order spaces:

  > - the dependent dofs (usually on the higher-order side of a face/edge) can
  >   be treated just like the hanging/dependent dofs case.

- Quadrature point vector, **Q-vector**:

  > - this is similar to **E-vector** where instead of dofs, the vector represents
  >   values at quadrature points, grouped by element.

- In many cases it is useful to distinguish two types of vectors:

  > - **X-vector**, or **primal X-vector**, and **X'-vector**, or **dual X-vector**
  > - here X can be any of the T, L, E, or Q categories
  > - for example, the mass matrix operator maps a **T-vector** to a **T'-vector**
  > - the solutions vector is a **T-vector**, and the RHS vector is a **T'-vector**
  > - using the parallel prolongation operator, one can map the solution
  >   **T-vector** to a solution **L-vector**, etc.

Operator representation/storage/action categories:

- Full true-dof parallel assembly, **TA**, or **A**:

  > - ParCSR or similar format
  > - the T in TA indicates that the data format represents an operator from a
  >   **T-vector** to a **T'-vector**.

- Full local assembly, **LA**:

  > - CSR matrix on each rank
  > - the parallel prolongation operator, $\bm{P}$, (and its transpose) should use
  >   optimized matrix-free action
  > - note that $\bm{P}$ is the operator mapping T-vectors to L-vectors.

- Element matrix assembly, **EA**:

  > - each element matrix is stored as a dense matrix
  > - optimized element and parallel prolongation operators
  > - note that the element prolongation operator is the mapping from an
  >   **L-vector** to an **E-vector**.

- Quadrature-point/partial assembly, **QA** or **PA**:

  > - precompute and store $w\det(J)$ at all quadrature points in all mesh elements
  > - the stored data can be viewed as a **Q-vector**.

- Unassembled option,  **UA** or **U**:

  > - no assembly step
  > - the action uses directly the mesh node coordinates, and assumes specific
  >   form of the coefficient, e.g. constant, piecewise-constant, or given as a
  >   **Q-vector** (Q-coefficient).

### Partial Assembly

Since the global operator $\bm{A}$ is just a series of variational restrictions
with $\bm{B}$, $\bm{G}$ and $\bm{P}$, starting from its
point-wise kernel $\bm{D}$, a "matvec" with $\bm{A}$ can be
performed by evaluating and storing some of the innermost variational restriction
matrices, and applying the rest of the operators "on-the-fly". For example, one can
compute and store a global matrix on **T-vector** level. Alternatively, one can compute
and store only the subdomain (**L-vector**) or element (**E-vector**) matrices and
perform the action of $\bm{A}$ using matvecs with $\bm{P}$ or
$\bm{P}$ and $\bm{G}$. While these options are natural for
low-order discretizations, they are not a good fit for high-order methods due to
the amount of FLOPs needed for their evaluation, as well as the memory transfer
needed for a matvec.

Our focus in libCEED, instead, is on **partial assembly**, where we compute and
store only $\bm{D}$ (or portions of it) and evaluate the actions of
$\bm{P}$, $\bm{G}$ and $\bm{B}$ on-the-fly.
Critically for performance, we take advantage of the tensor-product structure of the
degrees of freedom and quadrature points on *quad* and *hex* elements to perform the
action of $\bm{B}$ without storing it as a matrix.

Implemented properly, the partial assembly algorithm requires optimal amount of
memory transfers (with respect to the polynomial order) and near-optimal FLOPs
for operator evaluation. It consists of an operator *setup* phase, that
evaluates and stores $\bm{D}$ and an operator *apply* (evaluation) phase that
computes the action of $\bm{A}$ on an input vector. When desired, the setup
phase may be done as a side-effect of evaluating a different operator, such as a
nonlinear residual. The relative costs of the setup and apply phases are
different depending on the physics being expressed and the representation of
$\bm{D}$.

### Parallel Decomposition

After the application of each of the first three transition operators,
$\bm{P}$, $\bm{G}$ and $\bm{B}$, the operator evaluation
is decoupled  on their ranges, so $\bm{P}$, $\bm{G}$ and
$\bm{B}$ allow us to "zoom-in" to subdomain, element and quadrature point
level, ignoring the coupling at higher levels.

Thus, a natural mapping of $\bm{A}$ on a parallel computer is to split the
**T-vector** over MPI ranks (a non-overlapping decomposition, as is typically
used for sparse matrices), and then split the rest of the vector types over
computational devices (CPUs, GPUs, etc.) as indicated by the shaded regions in
the diagram above.

One of the advantages of the decomposition perspective in these settings is that
the operators $\bm{P}$, $\bm{G}$, $\bm{B}$ and
$\bm{D}$ clearly separate the MPI parallelism
in the operator ($\bm{P}$) from the unstructured mesh topology
($\bm{G}$), the choice of the finite element space/basis ($\bm{B}$)
and the geometry and point-wise physics $\bm{D}$. These components also
naturally fall in different classes of numerical algorithms -- parallel (multi-device)
linear algebra for $\bm{P}$, sparse (on-device) linear algebra for
$\bm{G}$, dense/structured linear algebra (tensor contractions) for
$\bm{B}$ and parallel point-wise evaluations for $\bm{D}$.

Currently in libCEED, it is assumed that the host application manages the global
**T-vectors** and the required communications among devices (which are generally
on different compute nodes) with **P**. Our API is thus focused on the
**L-vector** level, where the logical devices, which in the library are
represented by the {ref}`Ceed` object, are independent. Each MPI rank can use one or
more {ref}`Ceed`s, and each {ref}`Ceed`, in turn, can represent one or more physical
devices, as long as libCEED backends support such configurations. The idea is
that every MPI rank can use any logical device it is assigned at runtime. For
example, on a node with 2 CPU sockets and 4 GPUs, one may decide to use 6 MPI
ranks (each using a single {ref}`Ceed` object): 2 ranks using 1 CPU socket each, and
4 using 1 GPU each. Another choice could be to run 1 MPI rank on the whole node
and use 5 {ref}`Ceed` objects: 1 managing all CPU cores on the 2 sockets and 4
managing 1 GPU each. The communications among the devices, e.g. required for
applying the action of $\bm{P}$, are currently out of scope of libCEED. The
interface is non-blocking for all operations involving more than O(1) data,
allowing operations performed on a coprocessor or worker threads to overlap with
operations on the host.

## API Description

The libCEED API takes an algebraic approach, where the user essentially
describes in the *frontend* the operators **G**, **B** and **D** and the library
provides *backend* implementations and coordinates their action to the original
operator on **L-vector** level (i.e. independently on each device / MPI task).

One of the advantages of this purely algebraic description is that it already
includes all the finite element information, so the backends can operate on
linear algebra level without explicit finite element code. The frontend
description is general enough to support a wide variety of finite element
algorithms, as well as some other types algorithms such as spectral finite
differences. The separation of the front- and backends enables applications to
easily switch/try different backends. It also enables backend developers to
impact many applications from a single implementation.

Our long-term vision is to include a variety of backend implementations in
libCEED, ranging from reference kernels to highly optimized kernels targeting
specific devices (e.g. GPUs) or specific polynomial orders. A simple reference
backend implementation is provided in the file
[ceed-ref.c](https://github.com/CEED/libCEED/blob/main/backends/ref/ceed-ref.c).

On the frontend, the mapping between the decomposition concepts and the code
implementation is as follows:

- **L-**, **E-** and **Q-vector** are represented as variables of type {ref}`CeedVector`.
  (A backend may choose to operate incrementally without forming explicit **E-** or
  **Q-vectors**.)
- $\bm{G}$ is represented as variable of type {ref}`CeedElemRestriction`.
- $\bm{B}$ is represented as variable of type {ref}`CeedBasis`.
- the action of $\bm{D}$ is represented as variable of type {ref}`CeedQFunction`.
- the overall operator $\bm{G}^T \bm{B}^T \bm{D} \bm{B} \bm{G}$
  is represented as variable of type
  {ref}`CeedOperator` and its action is accessible through {c:func}`CeedOperatorApply()`.

To clarify these concepts and illustrate how they are combined in the API,
consider the implementation of the action of a simple 1D mass matrix
(cf. [tests/t500-operator.c](https://github.com/CEED/libCEED/blob/main/tests/t500-operator.c)).

```{literalinclude} ../../../tests/t500-operator.c
:language: c
:linenos: true
```

The constructor

```{literalinclude} ../../../tests/t500-operator.c
:end-at: CeedInit
:language: c
:start-at: CeedInit
```

creates a logical device `ceed` on the specified *resource*, which could also be
a coprocessor such as `"/nvidia/0"`. There can be any number of such devices,
including multiple logical devices driving the same resource (though performance
may suffer in case of oversubscription). The resource is used to locate a
suitable backend which will have discretion over the implementations of all
objects created with this logical device.

The `setup` routine above computes and stores $\bm{D}$, in this case a
scalar value in each quadrature point, while `mass` uses these saved values to perform
the action of $\bm{D}$. These functions are turned into the {ref}`CeedQFunction`
variables `qf_setup` and `qf_mass` in the {c:func}`CeedQFunctionCreateInterior()` calls:

```{literalinclude} ../../../tests/t500-operator.c
:end-before: //! [QFunction Create]
:language: c
:start-after: //! [QFunction Create]
```

A {ref}`CeedQFunction` performs independent operations at each quadrature point and
the interface is intended to facilitate vectorization.  The second argument is
an expected vector length. If greater than 1, the caller must ensure that the
number of quadrature points `Q` is divisible by the vector length. This is
often satisfied automatically due to the element size or by batching elements
together to facilitate vectorization in other stages, and can always be ensured
by padding.

In addition to the function pointers (`setup` and `mass`), {ref}`CeedQFunction`
constructors take a string representation specifying where the source for the
implementation is found. This is used by backends that support Just-In-Time
(JIT) compilation (i.e., CUDA and OCCA) to compile for coprocessors.
For full support across all backends, these {ref}`CeedQFunction` source files must only contain constructs mutually supported by C99, C++11, and CUDA.
For example, explicit type casting of void pointers and explicit use of compatible arguments for {code}`math` library functions is required, and variable-length array (VLA) syntax for array reshaping is only available via libCEED's {code}`CEED_Q_VLA` macro.

Different input and output fields are added individually, specifying the field
name, size of the field, and evaluation mode.

The size of the field is provided by a combination of the number of components
the effect of any basis evaluations.

The evaluation mode (see {ref}`CeedBasis-Typedefs and Enumerations`) `CEED_EVAL_INTERP`
for both input and output fields indicates that the mass operator only contains terms of
the form

$$
\int_\Omega v \cdot f_0 (u, \nabla u)
$$

where $v$ are test functions (see the {ref}`theoretical-framework`).
More general operators, such as those of the form

$$
\int_\Omega v \cdot f_0 (u, \nabla u) + \nabla v : f_1 (u, \nabla u)
$$

can be expressed.

For fields with derivatives, such as with the basis evaluation mode
(see {ref}`CeedBasis-Typedefs and Enumerations`) `CEED_EVAL_GRAD`, the size of the
field needs to reflect both the number of components and the geometric dimension.
A 3-dimensional gradient on four components would therefore mean the field has a size of
12\.

The $\bm{B}$ operators for the mesh nodes, `basis_x`, and the unknown field,
`basis_u`, are defined in the calls to the function {c:func}`CeedBasisCreateTensorH1Lagrange()`.
In this example, both the mesh and the unknown field use $H^1$ Lagrange finite
elements of order 1 and 4 respectively (the `P` argument represents the number of 1D
degrees of freedom on each element). Both basis operators use the same integration rule,
which is Gauss-Legendre with 8 points (the `Q` argument).

```{literalinclude} ../../../tests/t500-operator.c
:end-before: //! [Basis Create]
:language: c
:start-after: //! [Basis Create]
```

Other elements with this structure can be specified in terms of the `Q×P`
matrices that evaluate values and gradients at quadrature points in one
dimension using {c:func}`CeedBasisCreateTensorH1()`. Elements that do not have tensor
product structure, such as symmetric elements on simplices, will be created
using different constructors.

The $\bm{G}$ operators for the mesh nodes, `elem_restr_x`, and the unknown field,
`elem_restr_u`, are specified in the {c:func}`CeedElemRestrictionCreate()`. Both of these
specify directly the dof indices for each element in the `ind_x` and `ind_u`
arrays:

```{literalinclude} ../../../tests/t500-operator.c
:end-before: //! [ElemRestr Create]
:language: c
:start-after: //! [ElemRestr Create]
```

```{literalinclude} ../../../tests/t500-operator.c
:end-before: //! [ElemRestrU Create]
:language: c
:start-after: //! [ElemRestrU Create]
```

If the user has arrays available on a device, they can be provided using
`CEED_MEM_DEVICE`. This technique is used to provide no-copy interfaces in all
contexts that involve problem-sized data.

For discontinuous Galerkin and for applications such as Nek5000 that only
explicitly store **E-vectors** (inter-element continuity has been subsumed by
the parallel restriction $\bm{P}$), the element restriction $\bm{G}$
is the identity and {c:func}`CeedElemRestrictionCreateStrided()` is used instead.
We plan to support other structured representations of $\bm{G}$ which will
be added according to demand.
There are two common approaches for supporting non-conforming elements: applying the node constraints via $\bm P$ so that the **L-vector** can be processed uniformly and applying the constraints via $\bm G$ so that the **E-vector** is uniform.
The former can be done with the existing interface while the latter will require a generalization to element restriction that would define field values at constrained nodes as linear combinations of the values at primary nodes.

These operations, $\bm{P}$, $\bm{B}$, and $\bm{D}$,
are combined with a {ref}`CeedOperator`. As with {ref}`CeedQFunction`s, operator fields are added
separately with a matching field name, basis ($\bm{B}$), element restriction
($\bm{G}$), and **L-vector**. The flag
`CEED_VECTOR_ACTIVE` indicates that the vector corresponding to that field will
be provided to the operator when {c:func}`CeedOperatorApply()` is called. Otherwise the
input/output will be read from/written to the specified **L-vector**.

With partial assembly, we first perform a setup stage where $\bm{D}$ is evaluated
and stored. This is accomplished by the operator `op_setup` and its application
to `X`, the nodes of the mesh (these are needed to compute Jacobians at
quadrature points). Note that the corresponding {c:func}`CeedOperatorApply()` has no basis
evaluation on the output, as the quadrature data is not needed at the dofs:

```{literalinclude} ../../../tests/t500-operator.c
:end-before: //! [Setup Create]
:language: c
:start-after: //! [Setup Create]
```

```{literalinclude} ../../../tests/t500-operator.c
:end-before: //! [Setup Set]
:language: c
:start-after: //! [Setup Set]
```

```{literalinclude} ../../../tests/t500-operator.c
:end-before: //! [Setup Apply]
:language: c
:start-after: //! [Setup Apply]
```

The action of the operator is then represented by operator `op_mass` and its
{c:func}`CeedOperatorApply()` to the input **L-vector** `U` with output in `V`:

```{literalinclude} ../../../tests/t500-operator.c
:end-before: //! [Operator Create]
:language: c
:start-after: //! [Operator Create]
```

```{literalinclude} ../../../tests/t500-operator.c
:end-before: //! [Operator Set]
:language: c
:start-after: //! [Operator Set]
```

```{literalinclude} ../../../tests/t500-operator.c
:end-before: //! [Operator Apply]
:language: c
:start-after: //! [Operator Apply]
```

A number of function calls in the interface, such as {c:func}`CeedOperatorApply()`, are
intended to support asynchronous execution via their last argument,
`CeedRequest*`. The specific (pointer) value used in the above example,
`CEED_REQUEST_IMMEDIATE`, is used to express the request (from the user) for the
operation to complete before returning from the function call, i.e. to make sure
that the result of the operation is available in the output parameters
immediately after the call. For a true asynchronous call, one needs to provide
the address of a user defined variable. Such a variable can be used later to
explicitly wait for the completion of the operation.

## Gallery of QFunctions

LibCEED provides a gallery of built-in {ref}`CeedQFunction`s in the {file}`gallery/` directory.
The available QFunctions are the ones associated with the mass, the Laplacian, and
the identity operators. To illustrate how the user can declare a {ref}`CeedQFunction`
via the gallery of available QFunctions, consider the selection of the
{ref}`CeedQFunction` associated with a simple 1D mass matrix
(cf. [tests/t410-qfunction.c](https://github.com/CEED/libCEED/blob/main/tests/t410-qfunction.c)).

```{literalinclude} ../../../tests/t410-qfunction.c
:language: c
:linenos: true
```

## Interface Principles and Evolution

LibCEED is intended to be extensible via backends that are packaged with the
library and packaged separately (possibly as a binary containing proprietary
code). Backends are registered by calling

```{literalinclude} ../../../backends/ref/ceed-ref.c
:end-before: //! [Register]
:language: c
:start-after: //! [Register]
```

typically in a library initializer or "constructor" that runs automatically.
`CeedInit` uses this prefix to find an appropriate backend for the resource.

Source (API) and binary (ABI) stability are important to libCEED. Prior to
reaching version 1.0, libCEED does not implement strict [semantic versioning](https://semver.org) across the entire interface. However, user code,
including libraries of {ref}`CeedQFunction`s, should be source and binary
compatible moving from 0.x.y to any later release 0.x.z. We have less experience
with external packaging of backends and do not presently guarantee source or
binary stability, but we intend to define stability guarantees for libCEED 1.0.
We'd love to talk with you if you're interested in packaging backends
externally, and will work with you on a practical stability policy.
# Getting Started

```{include} ./README.md
:start-after: gettingstarted-inclusion-marker
```
# Developer Notes

## Style Guide

Please check your code for style issues by running

`make style`

In addition to those automatically enforced style rules, libCEED tends to follow the following code style conventions:

- Variable names: `snake_case`
- Strut members: `snake_case`
- Function and method names: `PascalCase` or language specific style
- Type names: `PascalCase` or language specific style
- Constant names: `CAPS_SNAKE_CASE` or language specific style

Also, documentation files should have one sentence per line to help make git diffs clearer and less disruptive.

## Clang-tidy

Please check your code for common issues by running

`make tidy`

which uses the `clang-tidy` utility included in recent releases of Clang.  This
tool is much slower than actual compilation (`make -j8` parallelism helps).  To
run on a single file, use

`make interface/ceed.c.tidy`

for example.  All issues reported by `make tidy` should be fixed.

## Include-What-You-Use

Header inclusion for source files should follow the principal of 'include what you use' rather than relying upon transitive `#include` to define all symbols.

Every symbol that is used in the source file `foo.c` should be defined in `foo.c`, `foo.h`, or in a header file `#include`d in one of these two locations.
Please check your code by running the tool [`include-what-you-use`](https://include-what-you-use.org/) to see recommendations for changes to your source.
Most issues reported by `include-what-you-use` should be fixed; however this rule is flexible to account for differences in header file organization in external libraries.
If you have `include-what-you-use` installed in a sibling directory to libCEED or set the environment variable `IWYU_CC`, then you can use the makefile target `make iwyu`.

Header files should be listed in alphabetical order, with installed headers preceding local headers and `ceed` headers being listed first.
The `ceed-f64.h` and `ceed-f32.h` headers should only be included in `ceed.h`.

```c
#include <ceed.h>
#include <ceed/backend.h>
#include <stdbool.h>
#include <string.h>
#include "ceed-avx.h"
```

## Shape

Backends often manipulate tensors of dimension greater than 2.  It is
awkward to pass fully-specified multi-dimensional arrays using C99 and
certain operations will flatten/reshape the tensors for computational
convenience.  We frequently use comments to document shapes using a
lexicographic ordering.  For example, the comment

```c
// u has shape [dim, num_comp, Q, num_elem]
```

means that it can be traversed as

```c
for (d=0; d<dim; d++)
  for (c=0; c<num_comp; c++)
    for (q=0; q<Q; q++)
      for (e=0; e<num_elem; e++)
        u[((d*num_comp + c)*Q + q)*num_elem + e] = ...
```

This ordering is sometimes referred to as row-major or C-style.  Note
that flattening such as

```c
// u has shape [dim, num_comp, Q*num_elem]
```

and

```c
// u has shape [dim*num_comp, Q, num_elem]
```

are purely implicit -- one just indexes the same array using the
appropriate convention.

## `restrict` Semantics

QFunction arguments can be assumed to have `restrict` semantics. That is, each input and output array must reside in distinct memory without overlap.

## CeedVector Array Access Semantics

Backend implementations are expected to separately track 'owned' and 'borrowed' memory locations.
Backends are responsible for freeing 'owned' memory; 'borrowed' memory is set by the user and backends only have read/write access to 'borrowed' memory.
For any given precision and memory type, a backend should only have 'owned' or 'borrowed' memory, not both.

Backends are responsible for tracking which memory locations contain valid data.
If the user calls {c:func}`CeedVectorTakeArray` on the only memory location that contains valid data, then the {ref}`CeedVector` is left in an *invalid state*.
To repair an *invalid state*, the user must set valid data by calling {c:func}`CeedVectorSetValue`, {c:func}`CeedVectorSetArray`, or {c:func}`CeedVectorGetArrayWrite`.

Some checks for consistency and data validity with {ref}`CeedVector` array access are performed at the interface level.
All backends may assume that array access will conform to these guidelines:

- Borrowed memory

  - {ref}`CeedVector` access to borrowed memory is set with {c:func}`CeedVectorSetArray` with `copy_mode = CEED_USE_POINTER` and revoked with {c:func}`CeedVectorTakeArray`.
    The user must first call {c:func}`CeedVectorSetArray` with `copy_mode = CEED_USE_POINTER` for the appropriate precision and memory type before calling {c:func}`CeedVectorTakeArray`.
  - {c:func}`CeedVectorTakeArray` cannot be called on a vector in a *invalid state*.

- Owned memory

  - Owned memory can be allocated by calling {c:func}`CeedVectorSetValue` or by calling {c:func}`CeedVectorSetArray` with `copy_mode = CEED_COPY_VALUES`.
  - Owned memory can be set by calling {c:func}`CeedVectorSetArray` with `copy_mode = CEED_OWN_POINTER`.
  - Owned memory can also be allocated by calling {c:func}`CeedVectorGetArrayWrite`.
    The user is responsible for manually setting the contents of the array in this case.

- Data validity

  - Internal syncronization and user calls to {c:func}`CeedVectorSync` cannot be made on a vector in an *invalid state*.
  - Calls to {c:func}`CeedVectorGetArray` and {c:func}`CeedVectorGetArrayRead` cannot be made on a vector in an *invalid state*.
  - Calls to {c:func}`CeedVectorSetArray` and {c:func}`CeedVectorSetValue` can be made on a vector in an *invalid state*.
  - Calls to {c:func}`CeedVectorGetArrayWrite` can be made on a vector in an *invalid* state.
    Data syncronization is not required for the memory location returned by {c:func}`CeedVectorGetArrayWrite`.
    The caller should assume that all data at the memory location returned by {c:func}`CeedVectorGetArrayWrite` is *invalid*.

## Internal Layouts

Ceed backends are free to use any **E-vector** and **Q-vector** data layout, to include never fully forming these vectors, so long as the backend passes the `t5**` series tests and all examples.
There are several common layouts for **L-vectors**, **E-vectors**, and **Q-vectors**, detailed below:

- **L-vector** layouts

  - **L-vectors** described by a {ref}`CeedElemRestriction` have a layout described by the `offsets` array and `comp_stride` parameter.
    Data for node `i`, component `j`, element `k` can be found in the **L-vector** at index `offsets[i + k*elem_size] + j*comp_stride`.
  - **L-vectors** described by a strided {ref}`CeedElemRestriction` have a layout described by the `strides` array.
    Data for node `i`, component `j`, element `k` can be found in the **L-vector** at index `i*strides[0] + j*strides[1] + k*strides[2]`.

- **E-vector** layouts

  - If possible, backends should use {c:func}`CeedElemRestrictionSetELayout()` to use the `t2**` tests.
    If the backend uses a strided **E-vector** layout, then the data for node `i`, component `j`, element `k` in the **E-vector** is given by `i*layout[0] + j*layout[1] + k*layout[2]`.
  - Backends may choose to use a non-strided **E-vector** layout; however, the `t2**` tests will not function correctly in this case and the tests will need to be whitelisted for the backend to pass the test suite.

- **Q-vector** layouts

  - When the size of a {ref}`CeedQFunction` field is greater than `1`, data for quadrature point `i` component `j` can be found in the **Q-vector** at index `i + Q*j`.
    Backends are free to provide the quadrature points in any order.
  - When the {ref}`CeedQFunction` field has `emode` `CEED_EVAL_GRAD`, data for quadrature point `i`, component `j`, derivative `k` can be found in the **Q-vector** at index `i + Q*j + Q*size*k`.
  - Note that backend developers must take special care to ensure that the data in the **Q-vectors** for a field with `emode` `CEED_EVAL_NONE` is properly ordered when the backend uses different layouts for **E-vectors** and **Q-vectors**.

## Backend Inheritance

There are three mechanisms by which a Ceed backend can inherit implementation from another Ceed backend.
These options are set in the backend initialization routine.

1. Delegation - Developers may use {c:func}`CeedSetDelegate()` to set a backend that will provide the implementation of any unimplemented Ceed objects.
2. Object delegation  - Developers may use {c:func}`CeedSetObjectDelegate()` to set a backend that will provide the implementation of a specific unimplemented Ceed object.
   Object delegation has higher precedence than delegation.
3. Operator fallback - Developers may use {c:func}`CeedSetOperatorFallbackResource()` to set a {ref}`Ceed` resource that will provide the implementation of unimplemented {ref}`CeedOperator` methods.
   A fallback {ref}`Ceed` with this resource will only be instantiated if a method is called that is not implemented by the parent {ref}`Ceed`.
   In order to use the fallback mechanism, the parent {ref}`Ceed` and fallback resource must use compatible **E-vector** and **Q-vector** layouts.
# Introduction

Historically, conventional high-order finite element methods were rarely used for
industrial problems because the Jacobian rapidly loses sparsity as the order is
increased, leading to unaffordable solve times and memory requirements
{cite}`brown2010`. This effect typically limited the order of accuracy to at most
quadratic, especially because quadratic finite element formulations are computationally advantageous in terms of
floating point operations (FLOPS) per degree of freedom (DOF)---see
{numref}`fig-assembledVsmatrix-free`---, despite the fast convergence and favorable
stability properties offered by higher order discretizations. Nowadays, high-order
numerical methods, such as the spectral element method (SEM)---a special case of
nodal p-Finite Element Method (FEM) which can reuse the interpolation nodes for
quadrature---are employed, especially with (nearly) affine elements, because
linear constant coefficient problems can be very efficiently solved using the
fast diagonalization method combined with a multilevel coarse solve. In
{numref}`fig-assembledVsmatrix-free` we analyze and compare the theoretical costs,
of different configurations: assembling the sparse matrix representing the action
of the operator (labeled as *assembled*), non assembling the matrix and storing
only the metric terms needed as an operator setup-phase (labeled as *tensor-qstore*)
and non assembling  the matrix and computing the metric terms on the fly and storing
a compact representation of the linearization at quadrature points (labeled as
*tensor*). In the right panel, we show the cost in terms of FLOPS/DOF. This metric for
computational efficiency made sense historically, when the performance was mostly
limited by processors' clockspeed. A more relevant performance plot for current
state-of-the-art high-performance machines (for which the bottleneck of performance is
mostly in the memory bandwith) is shown in the left panel of
{numref}`fig-assembledVsmatrix-free`, where the memory bandwith is measured in terms of
bytes/DOF. We can see that high-order methods, implemented properly with only partial
assembly, require optimal amount of memory transfers (with respect to the
polynomial order) and near-optimal FLOPs for operator evaluation. Thus, high-order
methods in matrix-free representation not only possess favorable properties, such as
higher accuracy and faster convergence to solution, but also manifest an efficiency gain
compared to their corresponding assembled representations.

(fig-assembledvsmatrix-free)=

:::{figure} ../../img/TensorVsAssembly.png
Comparison of memory transfer and floating point operations per
degree of freedom for different representations of a linear operator for a PDE in
3D with $b$ components and variable coefficients arising due to Newton
linearization of a material nonlinearity. The representation labeled as *tensor*
computes metric terms on the fly and stores a compact representation of the
linearization at quadrature points. The representation labeled as *tensor-qstore*
pulls the metric terms into the stored representation. The *assembled* representation
uses a (block) CSR format.
:::

Furthermore, software packages that provide high-performance implementations have often
been special-purpose and intrusive. libCEED {cite}`libceed-joss-paper` is a new library that offers a purely
algebraic interface for matrix-free operator representation and supports run-time
selection of implementations tuned for a variety of computational device types,
including CPUs and GPUs. libCEED's purely algebraic interface can unobtrusively be
integrated in new and legacy software to provide performance portable interfaces.
While libCEED's focus is on high-order finite elements, the approach is algebraic
and thus applicable to other discretizations in factored form. libCEED's role, as
a lightweight portable library that allows a wide variety of applications to share
highly optimized discretization kernels, is illustrated in
{numref}`fig-libCEED-backends`, where a non-exhaustive list of specialized
implementations (backends) is provided. libCEED provides a low-level Application
Programming Interface (API) for user codes so that applications with their own
discretization infrastructure (e.g., those in [PETSc](https://www.mcs.anl.gov/petsc/),
[MFEM](https://mfem.org/) and [Nek5000](https://nek5000.mcs.anl.gov/)) can evaluate
and use the core operations provided by libCEED. GPU implementations are available via
pure [CUDA](https://developer.nvidia.com/about-cuda) as well as the
[OCCA](http://github.com/libocca/occa) and [MAGMA](https://bitbucket.org/icl/magma)
libraries. CPU implementations are available via pure C and AVX intrinsics as well as
the [LIBXSMM](http://github.com/hfp/libxsmm) library. libCEED provides a unified
interface, so that users only need to write a single source code and can select the
desired specialized implementation at run time. Moreover, each process or thread can
instantiate an arbitrary number of backends.

(fig-libceed-backends)=

:::{figure} ../../img/libCEEDBackends.png
The role of libCEED as a lightweight, portable library which provides a low-level
API for efficient, specialized implementations. libCEED allows different applications
to share highly optimized discretization kernels.
:::
# libCEED Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, religion, or sexual identity
and orientation.

We pledge to act and interact in ways that contribute to an open, welcoming,
diverse, inclusive, and healthy community.

## Our Standards

Examples of behavior that contributes to a positive environment for our
community include:

* Demonstrating empathy and kindness toward other people
* Being respectful of differing opinions, viewpoints, and experiences
* Giving and gracefully accepting constructive feedback
* Accepting responsibility and apologizing to those affected by our mistakes,
  and learning from the experience
* Focusing on what is best not just for us as individuals, but for the
  overall community

Examples of unacceptable behavior include:

* The use of sexualized language or imagery, and sexual attention or
  advances of any kind
* Trolling, insulting or derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or email
  address, without their explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

## Enforcement Responsibilities

Community leaders are responsible for clarifying and enforcing our standards of
acceptable behavior and will take appropriate and fair corrective action in
response to any behavior that they deem inappropriate, threatening, offensive,
or harmful.

Community leaders have the right and responsibility to remove, edit, or reject
comments, commits, code, wiki edits, issues, and other contributions that are
not aligned to this Code of Conduct, and will communicate reasons for moderation
decisions when appropriate.

## Scope

This Code of Conduct applies within all community spaces, and also applies when
an individual is officially representing the community in public spaces.
Examples of representing our community include using an official e-mail address,
posting via an official social media account, or acting as an appointed
representative at an online or offline event.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported to the community leaders responsible for enforcement at
jed@jedbrown.org, valeria.barra@colorado.edu, or tzanio@llnl.gov.
All complaints will be reviewed and investigated promptly and fairly.

All community leaders are obligated to respect the privacy and security of the
reporter of any incident.

## Enforcement Guidelines

Community leaders will follow these Community Impact Guidelines in determining
the consequences for any action they deem in violation of this Code of Conduct:

### 1. Correction

**Community Impact**: Use of inappropriate language or other behavior deemed
unprofessional or unwelcome in the community.

**Consequence**: A private, written warning from community leaders, providing
clarity around the nature of the violation and an explanation of why the
behavior was inappropriate. A public apology may be requested.

### 2. Warning

**Community Impact**: A violation through a single incident or series
of actions.

**Consequence**: A warning with consequences for continued behavior. No
interaction with the people involved, including unsolicited interaction with
those enforcing the Code of Conduct, for a specified period of time. This
includes avoiding interactions in community spaces as well as external channels
like social media. Violating these terms may lead to a temporary or
permanent ban.

### 3. Temporary Ban

**Community Impact**: A serious violation of community standards, including
sustained inappropriate behavior.

**Consequence**: A temporary ban from any sort of interaction or public
communication with the community for a specified period of time. No public or
private interaction with the people involved, including unsolicited interaction
with those enforcing the Code of Conduct, is allowed during this period.
Violating these terms may lead to a permanent ban.

### 4. Permanent Ban

**Community Impact**: Demonstrating a pattern of violation of community
standards, including sustained inappropriate behavior,  harassment of an
individual, or aggression toward or disparagement of classes of individuals.

**Consequence**: A permanent ban from any sort of public interaction within
the community.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 2.0, available at
https://www.contributor-covenant.org/version/2/0/code_of_conduct.html.

Community Impact Guidelines were inspired by [Mozilla's code of conduct
enforcement ladder](https://github.com/mozilla/diversity).

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see the FAQ at
https://www.contributor-covenant.org/faq. Translations are available at
https://www.contributor-covenant.org/translations.
# libCEED: How to Contribute

Contributions to libCEED are encouraged.
<!---
Please use a pull request to the appropriate branch ('stable' for
backward-compatible bug fixes for the last stable release, main' for
new features and everything else).
-->
Please make your commits well-organized and
[atomic](https://en.wikipedia.org/wiki/Atomic_commit#Atomic_commit_convention),
using `git rebase --interactive` as needed.  Check that tests
(including "examples") pass using `make prove-all`.  If adding a new
feature, please add or extend a test so that your new feature is
tested.

In typical development, every commit should compile, be covered by the
test suite, and pass all tests.  This improves the efficiency of
reviewing and facilitates use of
[`git bisect`](https://git-scm.com/docs/git-bisect).

Open an issue or RFC (request for comments) pull request to discuss
any significant changes before investing time.  It is useful to create
a WIP (work in progress) pull request for any long-running development
so that others can be aware of your work and help to avoid creating
merge conflicts.

Write commit messages for a reviewer of your pull request and for a
future developer (maybe you) that bisects and finds that a bug was
introduced in your commit.  The assumptions that are clear in your
mind while committing are likely not in the mind of whomever (possibly
you) needs to understand it in the future.

Give credit where credit is due using tags such as `Reported-by:
Helpful User <helpful@example.com>` or
[`Co-authored-by: Snippet Mentor <code.by@comment.com>`](https://help.github.com/en/github/committing-changes-to-your-project/creating-a-commit-with-multiple-authors#creating-co-authored-commits-on-the-command-line).
Please use a real name and email for your author information (`git
config user.name` and `user.email`).  If your author information or
email becomes inconsistent (look at `git shortlog -se`), please edit
`.mailmap` to obtain your preferred name and email address.

When contributors make a major contribution and support it, their names
are included in the automatically generated user-manual documentation.

Please avoid "merging from upstream" (like merging 'main' into your
feature branch) unless there is a specific reason to do so, in which
case you should explain why in the merge commit.
[Rationale](https://lwn.net/Articles/328436/) from
[Junio](https://gitster.livejournal.com/42247.html) and
[Linus](http://yarchive.net/comp/linux/git_merges_from_upstream.html).

You can use `make style` to help conform to coding conventions of the
project, but try to avoid mixing whitespace or formatting changes with
content changes (see atomicity above).

By submitting a pull request, you are affirming the following.

## [Developer's Certificate of Origin 1.1](https://developercertificate.org/)

By making a contribution to this project, I certify that:

(a) The contribution was created in whole or in part by me and I
    have the right to submit it under the open source license
    indicated in the file; or

(b) The contribution is based upon previous work that, to the best
    of my knowledge, is covered under an appropriate open source
    license and I have the right under that license to submit that
    work with modifications, whether created in whole or in part
    by me, under the same open source license (unless I am
    permitted to submit under a different license), as indicated
    in the file; or

(c) The contribution was provided directly to me by some other
    person who certified (a), (b) or (c) and I have not modified
    it.

(d) I understand and agree that this project and the contribution
    are public and that a record of the contribution (including all
    personal information I submit with it, including my sign-off) is
    maintained indefinitely and may be redistributed consistent with
    this project or the open source license(s) involved.

## Authorship

libCEED contains components authored by many individuals.  It is
important that contributors receive appropriate recognition through
informal and academically-recognized credit systems such as
publications.  Status as a named author on the users manual and
libCEED software publications will be granted for those who

1. make significant contributions to libCEED (in implementation,
  documentation, conceptualization, review, etc.) and
2. maintain and support those contributions.

Maintainers will do their best to notice when contributions reach this
level and add your name to `AUTHORS`, but please email or create an
issue if you believe your contributions have met these criteria and
haven't yet been acknowledged.

Authors of publications about libCEED as a whole, including
DOI-bearing archives, shall offer co-authorship to all individuals
listed in the `AUTHORS` file.  Authors of publications claiming
specific libCEED contributions shall evaluate those listed in
`AUTHORS` and offer co-authorship to those who made significant
intellectual contributions to the work.

Note that there is no co-authorship expectation for those publishing
about use of libCEED (versus creation of new features in libCEED), but
see the [citing section](https://libceed.readthedocs.io/en/latest/gettingstarted/#how-to-cite)
and use your judgment regarding significance of support/advice you may have
received in developing your use case and interpreting results.
# Floating Point Precision

Currently, libCEED supports two options for {code}`CeedScalar` : double and single.  The default is to use 
double precision.  Users wishing to set {code}`CeedScalar` to single precision should edit `include/ceed/ceed.h` and change

```{literalinclude} ../../../include/ceed/ceed.h
:end-at: "#include \"ceed-f64.h\""
:language: c
:start-at: "#include \"ceed-f64.h\""
```

to include {code}`ceed-f32.h` instead, then recompile the library.
Tests can be run using `make test FC=` because the Fortran tests do not support single precision at this time.

## Language-specific notes

 - **C**: {code}`CEED_SCALAR_TYPE` will be defined to match one of the values of the {code}`CeedScalarType` {code}`enum`, and can be used 
       for compile-time checking of {code}`CeedScalar`'s type; see, e.g., {code}`tests/t314-basis.c`.

 - **Fortran**: There is no definition of {code}`CeedScalar` available in the Fortran header.  The user is responsible for ensuring
            that data used in Fortran code is of the correct type ({code}`real*8` or {code}`real*4`) for libCEED's current configuration.

 - **Julia**: After compiling the single precision version of libCEED, instruct LibCEED.jl to use this library with the {code}`set_libceed_path!`
              function and restart the Julia session. LibCEED.jl will configure itself to use the appropriate type for {code}`CeedScalar`. 

 - **Python**: Make sure to replace the {code}`ceed-f64.h` inclusion rather than commenting it out, to guarantee that the Python
           bindings will pick the correct precision.
           The {c:func}`scalar_type()` function has been added to the {code}`Ceed` class for convenience.  It returns a string 
           corresponding to a numpy datatype matching that of {code}`CeedScalar`.

 - **Rust**: The {code}`Scalar` type corresponds to {code}`CeedScalar`.

**This is work in progress!**  The ability to use single precision is an initial step in ongoing development of mixed-precision support in libCEED.
A current GitHub [issue](https://github.com/CEED/libCEED/issues/778) contains discussions related to this development.
# Julia, Python, and Rust Interfaces

libCEED provides high-level interfaces using the Julia, Python, and Rust
programming languages.

More information about the Julia interface can be found at the [LibCEED.jl
documentation](http://ceed.exascaleproject.org/libCEED-julia-docs/dev/).

Usage of the Python interface is illustrated through a sequence of [Jupyter
Notebook tutorials](https://github.com/CEED/libCEED/tree/main/examples/python). More
information on the Python interface is available in the [SciPy paper](https://doi.org/10.25080/Majora-342d178e-00c).

More information about the Rust interface can be found at the [Rust interface
documentation](https://docs.rs/libceed).
# Changes/Release Notes

On this page we provide a summary of the main API changes, new features and examples
for each release of libCEED.

(main)=

## Current `main` branch

### Interface changes

- Update {c:func}`CeedQFunctionGetFields` and {c:func}`CeedOperatorGetFields` to include number of fields.
- Promote to the public API: QFunction and Operator field objects, `CeedQFunctionField` and `CeedOperatorField`, and associated getters, {c:func}`CeedQFunctionGetFields`; {c:func}`CeedQFunctionFieldGetName`; {c:func}`CeedQFunctionFieldGetSize`; {c:func}`CeedQFunctionFieldGetEvalMode`; {c:func}`CeedOperatorGetFields`; {c:func}`CeedOperatorFieldGetElemRestriction`; {c:func}`CeedOperatorFieldGetBasis`; and {c:func}`CeedOperatorFieldGetVector`.
- Clarify and document conditions where `CeedQFunction` and `CeedOperator` become immutable and no further fields or suboperators can be added.
- Add {c:func}`CeedOperatorLinearAssembleQFunctionBuildOrUpdate` to reduce object creation overhead in assembly of CeedOperator preconditioning ingredients.
- Promote {c:func}`CeedOperatorCheckReady`to the public API to facilitate interactive interfaces.
- Warning added when compiling OCCA backend to alert users that this backend is experimental.
- `ceed-backend.h`, `ceed-hash.h`, and `ceed-khash.h` removed. Users should use `ceed/backend.h`, `ceed/hash.h`, and `ceed/khash.h`.
- Added {c:func}`CeedQFunctionGetKernelName`; refactored {c:func}`CeedQFunctionGetSourcePath` to exclude function kernel name.
- Clarify documentation for {c:func}`CeedVectorTakeArray`; this function will error if {c:func}`CeedVectorSetArray` with `copy_mode == CEED_USE_POINTER` was not previously called for the corresponding `CeedMemType`.
- Added {c:func}`CeedVectorGetArrayWrite` that allows access to uninitalized arrays; require initalized data for {c:func}`CeedVectorGetArray`.
- Added {c:func}`CeedQFunctionContextRegisterDouble` and {c:func}`CeedQFunctionContextRegisterInt32` with {c:func}`CeedQFunctionContextSetDouble` and {c:func}`CeedQFunctionContextSetInt32` to facilitate easy updating of {c:struct}`CeedQFunctionContext` data by user defined field names.
- Added {c:func}`CeedQFunctionContextGetFieldDescriptions` to retreive user defined descriptions of fields that are registered with `CeedQFunctionContextRegister*`.

### New features

- `CeedScalar` can now be set as `float` or `double` at compile time.
- Added JiT utilities in `ceed/jit-tools.h` to reduce duplicated code in GPU backends.
- Added support for JiT of QFunctions with `#include "relative/path/local-file.h"` statements for additional local files. Note that files included with `""` are searched relative to the current file first, then by compiler paths (as with `<>` includes). To use this feature, one should adhere to relative paths only, not compiler flags like `-I`, which the JiT will not be aware of.
- Remove need to guard library headers in QFunction source for code generation backends.
- `CeedDebugEnv()` macro created to provide debugging outputs when Ceed context is not present.
- Added {c:func}`CeedStringAllocCopy` to reduce repeated code for copying strings internally.
- Added {c:func}`CeedPathConcatenate` to facilitate loading kernel source files with a path relative to the current file.

### Maintainability

- Refactored preconditioner support internally to facilitate future development and improve GPU completeness/test coverage.
- `Include-what-you-use` makefile target added as `make iwyu`.
- Create backend constant `CEED_FIELD_MAX` to reduce magic numbers in codebase.
- Put GPU JiTed kernel source code into separate files.

(v0-9)=

## v0.9 (Jul 6, 2021)

### Interface changes

- Minor modification in error handling macro to silence pedantic warnings when compiling with Clang, but no functional impact.

### New features

- Add {c:func}`CeedVectorAXPY` and {c:func}`CeedVectorPointwiseMult` as a convenience for stand-alone testing and internal use.
- Add `CEED_QFUNCTION_HELPER` macro to properly annotate QFunction helper functions for code generation backends.
- Add `CeedPragmaOptimizeOff` macro for code that is sensitive to floating point errors from fast math optimizations.
- Rust support: split `libceed-sys` crate out of `libceed` and [publish both on crates.io](https://crates.io/crates/libceed).

### Performance improvements

### Examples

- Solid mechanics mini-app updated to explore the performance impacts of various formulations in the initial and current configurations.
- Fluid mechanics example adds GPU support and improves modularity.

### Deprecated backends

- The `/cpu/self/tmpl` and `/cpu/self/tmpl/sub` backends have been removed. These backends were intially added to test the backend inheritance mechanism, but this mechanism is now widely used and tested in multiple backends.

(v0-8)=

## v0.8 (Mar 31, 2021)

### Interface changes

- Error handling improved to include enumerated error codes for C interface return values.
- Installed headers that will follow semantic versioning were moved to {code}`include/ceed` directory. These headers have been renamed from {code}`ceed-*.h` to {code}`ceed/*.h`. Placeholder headers with the old naming schema are currently provided, but these headers will be removed in the libCEED v0.9 release.

### New features

- Julia and Rust interfaces added, providing a nearly 1-1 correspondence with the C interface, plus some convenience features.
- Static libraries can be built with `make STATIC=1` and the pkg-config file is installed accordingly.
- Add {c:func}`CeedOperatorLinearAssembleSymbolic` and {c:func}`CeedOperatorLinearAssemble` to support full assembly of libCEED operators.

### Performance improvements

- New HIP MAGMA backends for hipMAGMA library users: `/gpu/hip/magma` and `/gpu/hip/magma/det`.
- New HIP backends for improved tensor basis performance: `/gpu/hip/shared` and `/gpu/hip/gen`.

### Examples

- {ref}`example-petsc-elasticity` example updated with traction boundary conditions and improved Dirichlet boundary conditions.
- {ref}`example-petsc-elasticity` example updated with Neo-Hookean hyperelasticity in current configuration as well as improved Neo-Hookean hyperelasticity exploring storage vs computation tradeoffs.
- {ref}`example-petsc-navier-stokes` example updated with isentropic traveling vortex test case, an analytical solution to the Euler equations that is useful for testing boundary conditions, discretization stability, and order of accuracy.
- {ref}`example-petsc-navier-stokes` example updated with support for performing convergence study and plotting order of convergence by polynomial degree.

(v0-7)=

## v0.7 (Sep 29, 2020)

### Interface changes

- Replace limited {code}`CeedInterlaceMode` with more flexible component stride {code}`compstride` in {code}`CeedElemRestriction` constructors.
  As a result, the {code}`indices` parameter has been replaced with {code}`offsets` and the {code}`nnodes` parameter has been replaced with {code}`lsize`.
  These changes improve support for mixed finite element methods.
- Replace various uses of {code}`Ceed*Get*Status` with {code}`Ceed*Is*` in the backend API to match common nomenclature.
- Replace {code}`CeedOperatorAssembleLinearDiagonal` with {c:func}`CeedOperatorLinearAssembleDiagonal` for clarity.
- Linear Operators can be assembled as point-block diagonal matrices with {c:func}`CeedOperatorLinearAssemblePointBlockDiagonal`, provided in row-major form in a {code}`ncomp` by {code}`ncomp` block per node.
- Diagonal assemble interface changed to accept a {ref}`CeedVector` instead of a pointer to a {ref}`CeedVector` to reduce memory movement when interfacing with calling code.
- Added {c:func}`CeedOperatorLinearAssembleAddDiagonal` and {c:func}`CeedOperatorLinearAssembleAddPointBlockDiagonal` for improved future integration with codes such as MFEM that compose the action of {ref}`CeedOperator`s external to libCEED.
- Added {c:func}`CeedVectorTakeAray` to sync and remove libCEED read/write access to an allocated array and pass ownership of the array to the caller.
  This function is recommended over {c:func}`CeedVectorSyncArray` when the {code}`CeedVector` has an array owned by the caller that was set by {c:func}`CeedVectorSetArray`.
- Added {code}`CeedQFunctionContext` object to manage user QFunction context data and reduce copies between device and host memory.
- Added {c:func}`CeedOperatorMultigridLevelCreate`, {c:func}`CeedOperatorMultigridLevelCreateTensorH1`, and {c:func}`CeedOperatorMultigridLevelCreateH1` to facilitate creation of multigrid prolongation, restriction, and coarse grid operators using a common quadrature space.

### New features

- New HIP backend: `/gpu/hip/ref`.
- CeedQFunction support for user `CUfunction`s in some backends

### Performance improvements

- OCCA backend rebuilt to facilitate future performance enhancements.
- Petsc BPs suite improved to reduce noise due to multiple calls to {code}`mpiexec`.

### Examples

- {ref}`example-petsc-elasticity` example updated with strain energy computation and more flexible boundary conditions.

### Deprecated backends

- The `/gpu/cuda/reg` backend has been removed, with its core features moved into `/gpu/cuda/ref` and `/gpu/cuda/shared`.

(v0-6)=

## v0.6 (Mar 29, 2020)

libCEED v0.6 contains numerous new features and examples, as well as expanded
documentation in [this new website](https://libceed.readthedocs.io).

### New features

- New Python interface using [CFFI](https://cffi.readthedocs.io/) provides a nearly
  1-1 correspondence with the C interface, plus some convenience features.  For instance,
  data stored in the {cpp:type}`CeedVector` structure are available without copy as
  {py:class}`numpy.ndarray`.  Short tutorials are provided in
  [Binder](https://mybinder.org/v2/gh/CEED/libCEED/main?urlpath=lab/tree/examples/tutorials/).
- Linear QFunctions can be assembled as block-diagonal matrices (per quadrature point,
  {c:func}`CeedOperatorAssembleLinearQFunction`) or to evaluate the diagonal
  ({c:func}`CeedOperatorAssembleLinearDiagonal`).  These operations are useful for
  preconditioning ingredients and are used in the libCEED's multigrid examples.
- The inverse of separable operators can be obtained using
  {c:func}`CeedOperatorCreateFDMElementInverse` and applied with
  {c:func}`CeedOperatorApply`.  This is a useful preconditioning ingredient,
  especially for Laplacians and related operators.
- New functions: {c:func}`CeedVectorNorm`, {c:func}`CeedOperatorApplyAdd`,
  {c:func}`CeedQFunctionView`, {c:func}`CeedOperatorView`.
- Make public accessors for various attributes to facilitate writing composable code.
- New backend: `/cpu/self/memcheck/serial`.
- QFunctions using variable-length array (VLA) pointer constructs can be used with CUDA
  backends.  (Single source is coming soon for OCCA backends.)
- Fix some missing edge cases in CUDA backend.

### Performance Improvements

- MAGMA backend performance optimization and non-tensor bases.
- No-copy optimization in {c:func}`CeedOperatorApply`.

### Interface changes

- Replace {code}`CeedElemRestrictionCreateIdentity` and
  {code}`CeedElemRestrictionCreateBlocked` with more flexible
  {c:func}`CeedElemRestrictionCreateStrided` and
  {c:func}`CeedElemRestrictionCreateBlockedStrided`.
- Add arguments to {c:func}`CeedQFunctionCreateIdentity`.
- Replace ambiguous uses of {cpp:enum}`CeedTransposeMode` for L-vector identification
  with {cpp:enum}`CeedInterlaceMode`.  This is now an attribute of the
  {cpp:type}`CeedElemRestriction` (see {c:func}`CeedElemRestrictionCreate`) and no
  longer passed as `lmode` arguments to {c:func}`CeedOperatorSetField` and
  {c:func}`CeedElemRestrictionApply`.

### Examples

libCEED-0.6 contains greatly expanded examples with {ref}`new documentation <Examples>`.
Notable additions include:

- Standalone {ref}`ex2-surface` ({file}`examples/ceed/ex2-surface`): compute the area of
  a domain in 1, 2, and 3 dimensions by applying a Laplacian.

- PETSc {ref}`example-petsc-area` ({file}`examples/petsc/area.c`): computes surface area
  of domains (like the cube and sphere) by direct integration on a surface mesh;
  demonstrates geometric dimension different from topological dimension.

- PETSc {ref}`example-petsc-bps`:

  - {file}`examples/petsc/bpsraw.c` (formerly `bps.c`): transparent CUDA support.
  - {file}`examples/petsc/bps.c` (formerly `bpsdmplex.c`): performance improvements
    and transparent CUDA support.
  - {ref}`example-petsc-bps-sphere` ({file}`examples/petsc/bpssphere.c`):
    generalizations of all CEED BPs to the surface of the sphere; demonstrates geometric
    dimension different from topological dimension.

- {ref}`example-petsc-multigrid` ({file}`examples/petsc/multigrid.c`): new p-multigrid
  solver with algebraic multigrid coarse solve.

- {ref}`example-petsc-navier-stokes` ({file}`examples/fluids/navierstokes.c`; formerly
  `examples/navier-stokes`): unstructured grid support (using PETSc's `DMPlex`),
  implicit time integration, SU/SUPG stabilization, free-slip boundary conditions, and
  quasi-2D computational domain support.

- {ref}`example-petsc-elasticity` ({file}`examples/solids/elasticity.c`): new solver for
  linear elasticity, small-strain hyperelasticity, and globalized finite-strain
  hyperelasticity using p-multigrid with algebraic multigrid coarse solve.

(v0-5)=

## v0.5 (Sep 18, 2019)

For this release, several improvements were made. Two new CUDA backends were added to
the family of backends, of which, the new `cuda-gen` backend achieves state-of-the-art
performance using single-source {ref}`CeedQFunction`. From this release, users
can define Q-Functions in a single source code independently of the targeted backend
with the aid of a new macro `CEED QFUNCTION` to support JIT (Just-In-Time) and CPU
compilation of the user provided {ref}`CeedQFunction` code. To allow a unified
declaration, the {ref}`CeedQFunction` API has undergone a slight change:
the `QFunctionField` parameter `ncomp` has been changed to `size`. This change
requires setting the previous value of `ncomp` to `ncomp*dim` when adding a
`QFunctionField` with eval mode `CEED EVAL GRAD`.

Additionally, new CPU backends
were included in this release, such as the `/cpu/self/opt/*` backends (which are
written in pure C and use partial **E-vectors** to improve performance) and the
`/cpu/self/ref/memcheck` backend (which relies upon the
[Valgrind](http://valgrind.org/) Memcheck tool to help verify that user
{ref}`CeedQFunction` have no undefined values).
This release also included various performance improvements, bug fixes, new examples,
and improved tests. Among these improvements, vectorized instructions for
{ref}`CeedQFunction` code compiled for CPU were enhanced by using `CeedPragmaSIMD`
instead of `CeedPragmaOMP`, implementation of a {ref}`CeedQFunction` gallery and
identity Q-Functions were introduced, and the PETSc benchmark problems were expanded
to include unstructured meshes handling were. For this expansion, the prior version of
the PETSc BPs, which only included data associated with structured geometries, were
renamed `bpsraw`, and the new version of the BPs, which can handle data associated
with any unstructured geometry, were called `bps`. Additionally, other benchmark
problems, namely BP2 and BP4 (the vector-valued versions of BP1 and BP3, respectively),
and BP5 and BP6 (the collocated versions---for which the quadrature points are the same
as the Gauss Lobatto nodes---of BP3 and BP4 respectively) were added to the PETSc
examples. Furthermoew, another standalone libCEED example, called `ex2`, which
computes the surface area of a given mesh was added to this release.

Backends available in this release:

| CEED resource (`-ceed`)  | Backend                                             |
|--------------------------|-----------------------------------------------------|
| `/cpu/self/ref/serial`   | Serial reference implementation                     |
| `/cpu/self/ref/blocked`  | Blocked reference implementation                    |
| `/cpu/self/ref/memcheck` | Memcheck backend, undefined value checks            |
| `/cpu/self/opt/serial`   | Serial optimized C implementation                   |
| `/cpu/self/opt/blocked`  | Blocked optimized C implementation                  |
| `/cpu/self/avx/serial`   | Serial AVX implementation                           |
| `/cpu/self/avx/blocked`  | Blocked AVX implementation                          |
| `/cpu/self/xsmm/serial`  | Serial LIBXSMM implementation                       |
| `/cpu/self/xsmm/blocked` | Blocked LIBXSMM implementation                      |
| `/cpu/occa`              | Serial OCCA kernels                                 |
| `/gpu/occa`              | CUDA OCCA kernels                                   |
| `/omp/occa`              | OpenMP OCCA kernels                                 |
| `/ocl/occa`              | OpenCL OCCA kernels                                 |
| `/gpu/cuda/ref`          | Reference pure CUDA kernels                         |
| `/gpu/cuda/reg`          | Pure CUDA kernels using one thread per element      |
| `/gpu/cuda/shared`       | Optimized pure CUDA kernels using shared memory     |
| `/gpu/cuda/gen`          | Optimized pure CUDA kernels using code generation   |
| `/gpu/magma`             | CUDA MAGMA kernels                                  |

Examples available in this release:

:::{list-table}
:header-rows: 1
:widths: auto
* - User code
  - Example
* - `ceed`
  - * ex1 (volume)
    * ex2 (surface)
* - `mfem`
  - * BP1 (scalar mass operator)
    * BP3 (scalar Laplace operator)
* - `petsc`
  - * BP1 (scalar mass operator)
    * BP2 (vector mass operator)
    * BP3 (scalar Laplace operator)
    * BP4 (vector Laplace operator)
    * BP5 (collocated scalar Laplace operator)
    * BP6 (collocated vector Laplace operator)
    * Navier-Stokes
* - `nek5000`
  - * BP1 (scalar mass operator)
    * BP3 (scalar Laplace operator)
:::

(v0-4)=

## v0.4 (Apr 1, 2019)

libCEED v0.4 was made again publicly available in the second full CEED software
distribution, release CEED 2.0. This release contained notable features, such as
four new CPU backends, two new GPU backends, CPU backend optimizations, initial
support for operator composition, performance benchmarking, and a Navier-Stokes demo.
The new CPU backends in this release came in two families. The `/cpu/self/*/serial`
backends process one element at a time and are intended for meshes with a smaller number
of high order elements. The `/cpu/self/*/blocked` backends process blocked batches of
eight interlaced elements and are intended for meshes with higher numbers of elements.
The `/cpu/self/avx/*` backends rely upon AVX instructions to provide vectorized CPU
performance. The `/cpu/self/xsmm/*` backends rely upon the
[LIBXSMM](http://github.com/hfp/libxsmm) package to provide vectorized CPU
performance. The `/gpu/cuda/*` backends provide GPU performance strictly using CUDA.
The `/gpu/cuda/ref` backend is a reference CUDA backend, providing reasonable
performance for most problem configurations. The `/gpu/cuda/reg` backend uses a simple
parallelization approach, where each thread treats a finite element. Using just in time
compilation, provided by nvrtc (NVidia Runtime Compiler), and runtime parameters, this
backend unroll loops and map memory address to registers. The `/gpu/cuda/reg` backend
achieve good peak performance for 1D, 2D, and low order 3D problems, but performance
deteriorates very quickly when threads run out of registers.

A new explicit time-stepping Navier-Stokes solver was added to the family of libCEED
examples in the `examples/petsc` directory (see {ref}`example-petsc-navier-stokes`).
This example solves the time-dependent Navier-Stokes equations of compressible gas
dynamics in a static Eulerian three-dimensional frame, using structured high-order
finite/spectral element spatial discretizations and explicit high-order time-stepping
(available in PETSc). Moreover, the Navier-Stokes example was developed using PETSc,
so that the pointwise physics (defined at quadrature points) is separated from the
parallelization and meshing concerns.

Backends available in this release:

| CEED resource (`-ceed`)  | Backend                                             |
|--------------------------|-----------------------------------------------------|
| `/cpu/self/ref/serial`   | Serial reference implementation                     |
| `/cpu/self/ref/blocked`  | Blocked reference implementation                    |
| `/cpu/self/tmpl`         | Backend template, defaults to `/cpu/self/blocked`   |
| `/cpu/self/avx/serial`   | Serial AVX implementation                           |
| `/cpu/self/avx/blocked`  | Blocked AVX implementation                          |
| `/cpu/self/xsmm/serial`  | Serial LIBXSMM implementation                       |
| `/cpu/self/xsmm/blocked` | Blocked LIBXSMM implementation                      |
| `/cpu/occa`              | Serial OCCA kernels                                 |
| `/gpu/occa`              | CUDA OCCA kernels                                   |
| `/omp/occa`              | OpenMP OCCA kernels                                 |
| `/ocl/occa`              | OpenCL OCCA kernels                                 |
| `/gpu/cuda/ref`          | Reference pure CUDA kernels                         |
| `/gpu/cuda/reg`          | Pure CUDA kernels using one thread per element      |
| `/gpu/magma`             | CUDA MAGMA kernels                                  |

Examples available in this release:

:::{list-table}
:header-rows: 1
:widths: auto
* - User code
  - Example
* - `ceed`
  - * ex1 (volume)
* - `mfem`
  - * BP1 (scalar mass operator)
    * BP3 (scalar Laplace operator)
* - `petsc`
  - * BP1 (scalar mass operator)
    * BP3 (scalar Laplace operator)
    * Navier-Stokes
* - `nek5000`
  - * BP1 (scalar mass operator)
    * BP3 (scalar Laplace operator)
:::

(v0-3)=

## v0.3 (Sep 30, 2018)

Notable features in this release include active/passive field interface, support for
non-tensor bases, backend optimization, and improved Fortran interface. This release
also focused on providing improved continuous integration, and many new tests with code
coverage reports of about 90%. This release also provided a significant change to the
public interface: a {ref}`CeedQFunction` can take any number of named input and output
arguments while {ref}`CeedOperator` connects them to the actual data, which may be
supplied explicitly to `CeedOperatorApply()` (active) or separately via
`CeedOperatorSetField()` (passive). This interface change enables reusable libraries
of CeedQFunctions and composition of block solvers constructed using
{ref}`CeedOperator`. A concept of blocked restriction was added to this release and
used in an optimized CPU backend. Although this is typically not visible to the user,
it enables effective use of arbitrary-length SIMD while maintaining cache locality.
This CPU backend also implements an algebraic factorization of tensor product gradients
to perform fewer operations than standard application of interpolation and
differentiation from nodes to quadrature points. This algebraic formulation
automatically supports non-polynomial and non-interpolatory bases, thus is more general
than the more common derivation in terms of Lagrange polynomials on the quadrature points.

Backends available in this release:

| CEED resource (`-ceed`) | Backend                                             |
|-------------------------|-----------------------------------------------------|
| `/cpu/self/blocked`     | Blocked reference implementation                    |
| `/cpu/self/ref`         | Serial reference implementation                     |
| `/cpu/self/tmpl`        | Backend template, defaults to `/cpu/self/blocked`   |
| `/cpu/occa`             | Serial OCCA kernels                                 |
| `/gpu/occa`             | CUDA OCCA kernels                                   |
| `/omp/occa`             | OpenMP OCCA kernels                                 |
| `/ocl/occa`             | OpenCL OCCA kernels                                 |
| `/gpu/magma`            | CUDA MAGMA kernels                                  |

Examples available in this release:

:::{list-table}
:header-rows: 1
:widths: auto
* - User code
  - Example
* - `ceed`
  - * ex1 (volume)
* - `mfem`
  - * BP1 (scalar mass operator)
    * BP3 (scalar Laplace operator)
* - `petsc`
  - * BP1 (scalar mass operator)
    * BP3 (scalar Laplace operator)
* - `nek5000`
  - * BP1 (scalar mass operator)
    * BP3 (scalar Laplace operator)
:::

(v0-21)=

## v0.21 (Sep 30, 2018)

A MAGMA backend (which relies upon the
[MAGMA](https://bitbucket.org/icl/magma) package) was integrated in libCEED for this
release. This initial integration set up the framework of using MAGMA and provided the
libCEED functionality through MAGMA kernels as one of libCEED’s computational backends.
As any other backend, the MAGMA backend provides extended basic data structures for
{ref}`CeedVector`, {ref}`CeedElemRestriction`, and {ref}`CeedOperator`, and implements
the fundamental CEED building blocks to work with the new data structures.
In general, the MAGMA-specific data structures keep the libCEED pointers to CPU data
but also add corresponding device (e.g., GPU) pointers to the data. Coherency is handled
internally, and thus seamlessly to the user, through the functions/methods that are
provided to support them.

Backends available in this release:

| CEED resource (`-ceed`) | Backend                         |
|-------------------------|---------------------------------|
| `/cpu/self`             | Serial reference implementation |
| `/cpu/occa`             | Serial OCCA kernels             |
| `/gpu/occa`             | CUDA OCCA kernels               |
| `/omp/occa`             | OpenMP OCCA kernels             |
| `/ocl/occa`             | OpenCL OCCA kernels             |
| `/gpu/magma`            | CUDA MAGMA kernels              |

Examples available in this release:

:::{list-table}
:header-rows: 1
:widths: auto
* - User code
  - Example
* - `ceed`
  - * ex1 (volume)
* - `mfem`
  - * BP1 (scalar mass operator)
    * BP3 (scalar Laplace operator)
* - `petsc`
  - * BP1 (scalar mass operator)
* - `nek5000`
  - * BP1 (scalar mass operator)
:::

(v0-2)=

## v0.2 (Mar 30, 2018)

libCEED was made publicly available the first full CEED software distribution, release
CEED 1.0. The distribution was made available using the Spack package manager to provide
a common, easy-to-use build environment, where the user can build the CEED distribution
with all dependencies. This release included a new Fortran interface for the library.
This release also contained major improvements in the OCCA backend (including a new
`/ocl/occa` backend) and new examples. The standalone libCEED example was modified to
compute the volume volume of a given mesh (in 1D, 2D, or 3D) and placed in an
`examples/ceed` subfolder. A new `mfem` example to perform BP3 (with the application
of the Laplace operator) was also added to this release.

Backends available in this release:

| CEED resource (`-ceed`) | Backend                         |
|-------------------------|---------------------------------|
| `/cpu/self`             | Serial reference implementation |
| `/cpu/occa`             | Serial OCCA kernels             |
| `/gpu/occa`             | CUDA OCCA kernels               |
| `/omp/occa`             | OpenMP OCCA kernels             |
| `/ocl/occa`             | OpenCL OCCA kernels             |

Examples available in this release:

:::{list-table}
:header-rows: 1
:widths: auto
* - User code
  - Example
* - `ceed`
  - * ex1 (volume)
* - `mfem`
  - * BP1 (scalar mass operator)
    * BP3 (scalar Laplace operator)
* - `petsc`
  - * BP1 (scalar mass operator)
* - `nek5000`
  - * BP1 (scalar mass operator)
:::

(v0-1)=

## v0.1 (Jan 3, 2018)

Initial low-level API of the CEED project. The low-level API provides a set of Finite
Elements kernels and components for writing new low-level kernels. Examples include:
vector and sparse linear algebra, element matrix assembly over a batch of elements,
partial assembly and action for efficient high-order operators like mass, diffusion,
advection, etc. The main goal of the low-level API is to establish the basis for the
high-level API. Also, identifying such low-level kernels and providing a reference
implementation for them serves as the basis for specialized backend implementations.
This release contained several backends: `/cpu/self`, and backends which rely upon the
[OCCA](http://github.com/libocca/occa) package, such as `/cpu/occa`,
`/gpu/occa`, and `/omp/occa`.
It also included several examples, in the `examples` folder:
A standalone code that shows the usage of libCEED (with no external
dependencies) to apply the Laplace operator, `ex1`; an `mfem` example to perform BP1
(with the application of the mass operator); and a `petsc` example to perform BP1
(with the application of the mass operator).

Backends available in this release:

| CEED resource (`-ceed`) | Backend                         |
|-------------------------|---------------------------------|
| `/cpu/self`             | Serial reference implementation |
| `/cpu/occa`             | Serial OCCA kernels             |
| `/gpu/occa`             | CUDA OCCA kernels               |
| `/omp/occa`             | OpenMP OCCA kernels             |

Examples available in this release:

| User code             | Example                           |
|-----------------------|-----------------------------------|
| `ceed`                | ex1 (scalar Laplace operator)     |
| `mfem`                | BP1 (scalar mass operator)        |
| `petsc`               | BP1 (scalar mass operator)        |
```
# libceed-sys: unsafe bindings to libCEED

This is the documentation for the low level (unsafe) Rust bindings to the libCEED C
interface. See the [libCEED user manual](https://libceed.readthedocs.io) for usage
information. Note that most Rust users will prefer the higher level (safe) Rust
interface in the [`libceed` crate](https://docs.rs/libceed).

libCEED is a low-level API for for the efficient high-order discretization methods
developed by the ECP co-design Center for Efficient Exascale Discretizations (CEED).
While our focus is on high-order finite elements, the approach is mostly algebraic
and thus applicable to other discretizations in factored form.

## Usage

To use low level libCEED bindings in a Rust package, the following `Cargo.toml`
can be used.
```toml
[dependencies]
libceed-sys = "0.9.0"
```

For a development version of the libCEED Rust bindings, use the following `Cargo.toml`.
```toml
[dependencies]
libceed-sys = { git = "https://github.com/CEED/libCEED", branch = "main" }
```

Supported features:
* `static` (default): link to static libceed.a
* `system`: use libceed from a system directory (otherwise, install from source)

## Development

To develop libCEED, use `cargo build` in the `rust/libceed-sys` directory to
install a local copy and build the bindings.

If you need custom flags for the C project, we recommend using `make -C c-src
configure` to cache arguments in `c-src/config.mk`. If that file exists during
`cargo build` then edits will prompt recompilation of the bindings.

### Shared libraries
If one is developing libCEED C source and testing multiple language bindings at
once, a few seconds can be cut out of the edit/compile/test loop by disabling
the `static` feature and using

```bash
export LD_LIBRARY_PATH=$CEED_DIR/lib
export PKG_CONFIG_PATH=$CEED_DIR/lib/pkgconfig
```

#### Without system
If you disable the `static` feature and are not using a system version from a
standard path/somewhere that can be found by pkg-config, then you'll need to set
`LD_LIBRARY_PATH` to the appropriate target directory for doctests to be able to
find it. This might look like

```bash
export LD_LIBRARY_PATH=$CEED_DIR/target/debug/build/libceed-sys-d1ea22c6e1ad3f23/out/lib
```

where the precise hash value is printed during `cargo build --verbose` or you
can find it with `find target -name libceed.so`. This mode of development is
more fragile than the default (which uses static libraries).

Note that the `LD_LIBRARY_PATH` workarounds will become unnecessary if [this
issue](https://github.com/rust-lang/cargo/issues/1592) is resolved -- it's
currently closed, but the problem still exists.

## License: BSD-2-Clause

## Contributing

The `libceed-sys` crate is developed within the [libCEED
repository](https://github.com/CEED/libCEED). See the [contributing
guidelines](https://libceed.readthedocs.io/en/latest/CONTRIBUTING/) for details.
# libceed: efficient, extensible discretization

[![GitHub Actions](https://github.com/CEED/libCEED/actions/workflows/rust-test-with-style.yml/badge.svg)](https://github.com/CEED/libCEED/actions/workflows/rust-test-with-style.yml)
[![Documentation](https://docs.rs/libceed/badge.svg)](https://docs.rs/libceed)

This crate provides an interface to [libCEED](https://libceed.readthedocs.io), which is a performance-portable library for extensible element-based discretization for partial differential equations and related computational problems.
The formulation is algebraic and intended to be lightweight and easy to incorporate in higher level abstractions.
See the [libCEED user manual](https://libceed.readthedocs.io) for details on [interface concepts](https://libceed.readthedocs.io/en/latest/libCEEDapi/) and extensive examples.

![libCEED operator decomposition](https://libceed.readthedocs.io/en/latest/_images/libCEED.png)

## Usage

To call libCEED from a Rust package, the following `Cargo.toml` can be used.
```toml
[dependencies]
libceed = "0.9.0"
```

For a development version of the libCEED Rust bindings, use the following `Cargo.toml`.
```toml
[dependencies]
libceed = { git = "https://github.com/CEED/libCEED", branch = "main" }
```

```rust
extern crate libceed;

fn main() -> libceed::Result<()> {
    let ceed = libceed::Ceed::init("/cpu/self/ref");
    let xc = ceed.vector_from_slice(&[0., 0.5, 1.0])?;
    let xs = xc.view()?;
    assert_eq!(xs[..], [0., 0.5, 1.0]);
    Ok(())
}
```

This crate provides modules for each object, but they are usually created from the `Ceed` object as with the vector above.
The resource string passed to `Ceed::init` is used to identify the "backend", which includes algorithmic strategies and hardware such as NVIDIA and AMD GPUs.
See the [libCEED documentation](https://libceed.readthedocs.io/en/latest/gettingstarted/#backends) for more information on available backends.

## Examples

Examples of libCEED can be found in the [libCEED repository](https://github.com/CEED/libCEED) under the `examples/rust` directory.

## Documentation

This crate uses `katexit` to render equations in the documentation.
To build the [documentation](https://docs.rs/libceed) locally with `katexit` enabled, use

```bash
cargo doc --features=katexit
```

## License: BSD-2-Clause

## Contributing

The `libceed` crate is developed within the [libCEED repository](https://github.com/CEED/libCEED).
See the [contributing guidelines](https://libceed.readthedocs.io/en/latest/CONTRIBUTING/) for details.
