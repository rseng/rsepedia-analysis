# Roadmap of LibSWIFFT

The list below details, in no particular order, a number of directions for future development of LibSWIFFT. We welcome [contributions](CONTRIBUTING.md) in each of these as well as suggestions for additional directions.

- Support for higher-modulus variations of SWIFFT.
- Improved build modularity: multi-versioned build installation, packaging for availability via `find_package` in cmake.
- Build-support for additional platforms, operating systems and toolchains.
- Improved test coverage: numerical edge cases.
- Automatic runtime selection of the best SWIFFT implementation for the native platform.
- Support for parallel processing using OpenMP.
- SWIFFT based hashing for long inputs.
- GPU kernels for SWIFFT functions.

## Out of Scope for LibSWIFFT

The list below gives a few examples of relevant directions that are out of scope for LibSWIFFT and should go into separate projects.

- Wrappers of LibSWIFFT in other programming languages.
- Integrations of LibSWIFFT with established cryptography libraries.
- Implementations of proofs of knowledge or digital signatures based on LibSWIFFT.
# Release Checklist

To prepare a pull request for release:

- All code changes should already be in their own commits.
- Decide the version of this release.
- Run the tests-executable to ensure all tests pass.
- Run performance tests (see section below) on appropriate machines to ensure performance is no worse than that of the previous version.
- Run coverage tests (see section below) to ensure code is sufficiently covered.
- Run doxygen (see section below) to test it correctly processes the source code.
- Change the version in `CMakeLists.txt`.
- Repeat the above tests.
- If any test or doxygen issues are found, abort the release until they are fixed.
- Commit with a comment like `version 1.0.0` using the decided version.
- Tag with the decided version preceded by `v`, e.g., `v1.0.0`.
- Create a pull-request.

To prepare a pull request for development, use tag with suffix `+dev`, e.g., `v1.0.0+dev`.

## Running Performance Tests

Performance tests are executed using a release build:

    mkdir -p build/release
    cd build/release
    cmake -DCMAKE_BUILD_TYPE=Release ../..
    make
    ./test/swifft_catch "[swifftperf]"

For an OpenMP release build, add the option `-DSWIFFT_ENABLE_OPENMP=On` to the above `cmake` command line.

## Running Coverage Tests

Coverage tests are executed using a debug build:

    mkdir -p build/debug
    cd build/debug
    cmake -DCMAKE_BUILD_TYPE=Debug ../..
    make
    ./test/swifft_catch
    for file in $(find . -name "*.gcda"); do echo "=== $file ==="; gcov $file; done | less
    for file in $(find . -name "*swifft*.gcov"); do echo "=== $file ==="; cat $file; done | less

This procedure is specific to GCC builds and does not cover performance testing code.

## Generating Documentation

To generate documentation:

- Ensure `sphinx` documentation generator is installed and working
- Ensure `doxygen` is installed and working
- Go to the `docs` directory
- Run `make html`

The doxygen configuration was set up using these steps:

- Run `doxygen -g doxygen.conf`
- Edit `doxygen.conf`:
  - Set `PROJECT_NAME` to `LibSWIFFT`
  - Set `PROJECT_NUMBER` to the decided version
  - Set `RECURSIVE` to `YES`.
  - Set `EXTRACT_ALL` to `YES`
  - Set `EXCLUDE` to `build`

The doxygen generation can be ran using `doxygen Doxyfile`.
# Intellectual Property

The original public version (v1.0.0) of LibSWIFFT is copyright (C) 2020 Yaron Gvili and [Gvili Tech Ltd](http://gvili.tech). Later versions may include contributions with separate copyrights.

The LibSWIFFT [logo](assets/LibSWIFFT-logo.png) and [icon](assets/LibSWIFFT-icon.png) are trademarks of Gvili Tech Ltd.

LibSWIFFT is open source software. For licensing terms, see [LICENSE.txt](LICENSE.txt).
![LibSWIFFT](assets/LibSWIFFT-logo.png)

# LibSWIFFT - A fast C/C++ library for the SWIFFT secure homomorphic hash function

- [Official Repository](https://github.com/gvilitechltd/LibSWIFFT)
- [JOSS paper](https://doi.org/10.21105/joss.03040) [![DOI](https://joss.theoj.org/papers/10.21105/joss.03040/status.svg)](https://doi.org/10.21105/joss.03040)
- [Documentation](http://libswifft.readthedocs.io/en/latest/)
- [![Documentation Status](https://readthedocs.org/projects/libswifft/badge/?version=latest)](http://libswifft.readthedocs.io/en/latest/)
- [![Packaging Status](https://github.com/gvilitechltd/libswifft/actions/workflows/new-code.yml/badge.svg?branch=main)](https://github.com/gvilitechltd/libswifft/actions/workflows/new-code.yml/)
- **Quick start on Linux**: Ensure `docker` is installed and runnable, clone this repository and go to its root directory, and run the command `docker build . && docker run --rm -it $(docker build -q .)` to build and test LibSWIFFT.
- **Quick performance comparison on Linux**: Ensure `docker` is installed and runnable, clone this repository and go to its root directory, and run the command `docker build . -f Dockerfile.compare-to-K2SN-MSS && docker run --rm -it $(docker build -q . -f Dockerfile.compare-to-K2SN-MSS)` to build and compare performance with K2SN-MSS.

## Introducing LibSWIFFT

LibSWIFFT is a production-ready C/C++ library providing SWIFFT, one of the fastest available secure hash functions, which is also collision-resistant as well as facilitates zero-knowledge proofs of knowledge of a preimage (ZKPoKP) and post-quantum digital signatures. It is based on academic work from 2007 described further below.

**Why now, in early 2021?** In late 2017, NIST has started a [process for standardizing post-quantum cryptography](https://csrc.nist.gov/Projects/post-quantum-cryptography/post-quantum-cryptography-standardization), suggesting that it believes it may not be too long before a practical quantum-computer that threatens critical security standards (including Internet ones) based on classical cryptography will become a reality. About two years later, Google announced it had achieved [quantum supremacy](https://en.wikipedia.org/wiki/Quantum_supremacy), by completing in 200 seconds a task they claimed would have taken a classical supercomputer about 10,000 years to complete. Though IBM, maker of the most powerful supercomputer at the time, disputed this claim and asserted the supercomputer would take only about 2.5 days for the task, it is clear quantum computing technology is advancing quickly. Consequently, post-quantum cryptography is becoming more relevant today and perhaps even urgent to develop.

**Why another implementation of SWIFFT?** LibSWIFFT is a reliable building block for fast and scalable cryptographic protocols. It is simple to use and maintain, has clean APIs, is well documented and tested, and is at least as fast as other implementations of SWIFFT and often faster. Other implementations of SWIFFT are:

- The [8-bit](https://github.com/anon1985/Swifft-avx2-8) and [16-bit](https://github.com/anon1985/K2SN-MSS/tree/master/swifft16) AVX2 implementations for [K2SN-MSS](https://eprint.iacr.org/2019/442.pdf). Both are arguably not as easy to use nor as well documented and tested as LibSWIFFT. The former is slower and uses less memory while the latter is about as fast as LibSWIFFT for AVX2 yet does not support AVX512.
- [The original implementation](https://github.com/micciancio/SWIFFT) written in 2007. It is minimal non-production code. [The AVX2 implementation for K2SN-MSS](https://eprint.iacr.org/2019/442.pdf) is reported to be 25% faster.

An invocation of the tests-executable of LibSWIFFT running single-threaded using AVX2 on an Intel Skylake microarchitecture (Intel(R) Core(TM) i7-10875H CPU @ 2.30GHz):

    $ ./test/swifft_catch "swifft takes at most 2000 cycles per call"
    Filters: swifft takes at most 2000 cycles per call
    running 1*10000000 rounds: cycles/rounds=1097.94 cycles/byte=4.28882 Giga-cycles/sec=2.30399 MB/sec=512.322 cycles/rdtsc=16

demonstrates that LibSWIFFT is quite fast on short inputs (here, 256 bytes), often used in practical zero-knowledge proofs and post-quantum digital signatures. This is more than an order of magnitude faster than the [originally reported](https://www.alonrosen.net/PAPERS/lattices/swifft.pdf) 40MB/sec on a 3.2 GHz Intel Pentium 4. This is also faster than [K2SN-MSS's binary 16-bit SWIFFT function implementation](https://github.com/gvilitechltd/K2SN-MSS/tree/swifftperf) (for an input of 128 bytes), which is the fastest one in the K2SN-MSS implementation, for the same executaion settings, i.e. running single-threaded using AVX2 on an Intel Skylake microarchitecture (Intel(R) Core(TM) i7-10875H CPU @ 2.30GHz):

    $ ./tester
    1000000 SWIFFT16 rounds: cycles/round=737.363098 cycles/byte=5.760649

It also compares well with modern hash functions:

- [Blake3](https://github.com/BLAKE3-team/BLAKE3) - cryptographic hash function achieving about [3-to-4 cycles/byte using AVX512 on short inputs](https://github.com/BLAKE3-team/BLAKE3-specs/blob/master/blake3.pdf) and are non-homomorphic nor facilitating proofs of knowledge of a preimage.
- [Seahash](https://docs.rs/seahash/4.0.1/seahash/index.html) - a hash function achieving ~0.24 cycles/byte but is non-cryptographic.

## On SWIFFT

[SWIFFT](https://en.wikipedia.org/wiki/SWIFFT) is a family of homomorphic hash functions provided as part of a candidate submission to the [NIST hash function competition](https://en.wikipedia.org/wiki/NIST_hash_function_competition) in 2008 by Vadim Lyubashevsky, Daniele Micciancio, Chris Peikert, and Alon Rosen. The family has an interesting set of characteristics:

1. **Provably collision-resistant**: it is provably computationally hard to find a binary-valued `x` given its image under `f` that is chosen randomly.
2. **Universal**: for any `x_1,x_2` in the domain the probability that `f(x_1) = f(x_2)` over the choice of `f` in the family is the inverse of the size of the range.
3. **Regular**: if `x` is distributed uniformly then `f(x)` is distributed uniformly for any `f` in the family.
4. **Randomness extractor**: even if `x` is not distributed uniformly, `f(x)` is distributed uniformly over the choice of `f` in the family.
5. **Homomorphic**: each `f` in the family is homomorphic over the domain.
6. **Facilitates ZKPoKPs**: the resulting hash value can be readily used in ZKPoKPs.
7. **Facilitates post-quantum digital signatures**: its homomorphism property enables short post-quantum hash-based digital signature schemes, such as [K2SN-MSS](https://eprint.iacr.org/2019/442.pdf).
8. **Constant-time**: any `f` in the family is free of data-dependent branching and therefore facilitates avoidance of timing side-channel attacks.

Nevertheless, the family is not pseudorandom, due to the homomorphism property, nor claimed to behave like a random oracle.

## On LibSWIFFT

LibSWIFFT was implemented with reference to the [SWIFFTX submission to NIST](https://csrc.nist.gov/projects/hash-functions/sha-3-project) and provides the same SWIFFT hash function that is part of the submission. High speed is achieved using various code optimization techniques, including SIMD instructions that are very natural for the implementation of the SWIFFT function. Compared to the SWIFFT code in the submission, LibSWIFFT adds the following:

1. Automatic library initialization using build-time generation of internal tables.
2. Convenient APIs, including for homomorphic operations and parallel variations based on OpenMP, for computing SWIFFT on short inputs.
3. Support for input vectors of either binary-valued (in {0,1}) or trinary-valued (in {-1,0,1}) elements.
4. Bug fixes with respect to the reference submission, in particular related to the homomorphism property.
5. Performance improvements compared to the reference submission.
6. Support for newer CPU instruction sets: AVX, AVX2, and AVX512.
7. Over 30 test-cases providing excellent coverage of the APIs and the mathematical properties of SWIFFT.

Formally, LibSWIFFT provides a single hash function that maps from an input domain `Z_2^{2048}` (taking 256B) to an output domain `Z_{257}^{64}` (taking 128B, at 2B per element) and then to a compact domain `Z_{256}^{64}` (taking 64B). The computation of the first map is done over `Z_{257}`. The homomorphism property applies to the input and output domains, but not to the compact domain, and is revealed when the binary-valued input domain is naturally embedded in `Z_{257}^{2048}`. Generally, it is computationally hard to find a binary-valued pre-image given an output computed as the sum of `N` outputs corresponding to known binary-valued pre-images. On the other hand, it is easy to find a small-valued pre-image (over `Z_{257}^{2048}`) when `N` is small, since it is simply the sum of the known pre-images due to the homomorphism property.

## Using LibSWIFFT

LibSWIFFT is intended to be used by cryptography researchers and by software developers knowledgeable in cryptography programming. LibSWIFFT is most useful in use cases that require provable-security and speed on short inputs. It may also be interesting in use cases that take advantage of its uncommon homomorphism property. Future versions of LibSWIFFT may target a larger audience.

The main LibSwifft C API is documented in `include/libswifft/swifft.h`. The following API variations are available:

- `include/libswifft/swifft_avx.h`: Same functions as in `include/libswifft/swifft.h` but with an added suffix `_AVX` and implemented using AVX instruction set.
- `include/libswifft/swifft_avx2.h`: Same functions as in `include/libswifft/swifft.h` but with an added suffix `_AVX2` and implemented using AVX2 instruction set.
- `include/libswifft/swifft_avx512.h`: Same functions as in `include/libswifft/swifft.h` but with an added suffix `_AVX512` and implemented using AVX512 instruction set.
- `include/libswifft/swifft.h`: Selects the implementations using the most advanced instruction set that was built into the library.

The version of LibSWIFFT is provided by the API in `include/libswifft/swifft_ver.h`.

The main LibSWIFFT C++ API is documented in `include/libswifft/swifft.hpp`.

Please refer to:
- the [release checklist document](RELEASE-CHECKLIST.md) for how to generate the documentation for the APIs using doxygen.
- the [code design document](CODE-DESIGN.md) for details on the architecture and design of the LibSWIFFT code.
- the [documentation](http://libswifft.readthedocs.io/en/latest/) for the complete details.

An extended use of the LibSWIFFT API follows the following steps:

- **Allocate buffers**: LibSWIFFT defines 3 types of buffers - input, output and compact. The input buffer `BitSequence input[SWIFFT_INPUT_BLOCK_SIZE]` holds a vector in `Z_2^{2048}` where each element takes 1 bit, the output buffer `BitSequence output[SWIFFT_OUTPUT_BLOCK_SIZE]` holds a vector in `Z_{257}^{64}` where each element takes 16 bits, and the compact buffer `BitSequence compact[SWIFFT_COMPACT_BLOCK_SIZE]` holds a value in `Z_{256}^64` taking 64 bytes.
- **Populate input buffers**: An input buffer is populated in preparation for hashing. This is normally done by directly setting the bits in the input buffer. Each bit corresponds to an element of the vector with a value in `{0,1}`.
- **Populate sign buffers**: A sign buffer is an input buffer whose bits are interpreted as sign bits. A 0-valued (resp. 1-valued) bit corresponds to a positive (resp. negative) sign. When an input buffer and a sign buffer are taken together, they define a vector in `{-1,0,1}^{2048}`.
- **Compute output buffers**: The hash of an input buffer, with or without a sign buffer, is computed into an output buffer. This is normally done using `SWIFFT_Compute` or `SWIFFT_ComputeSigned`.
- **Perform arithmetic operations with output buffers**: LibSWIFFT provides several arithemtic (homomorphic) operations involving output buffers whose result is put into an output buffer. The vectors of output buffers may be added, subtracted, or multiplied element-wise. See below for more details.
- **Compact the output buffer**: The hash in the output buffer is compacted into the compact buffer. This is an optional operation, in that the hash in the output buffer may be sufficient for certain applications.

A more restricted use of the LibSWIFFT API involves only the steps of allocating buffers, populating input buffers, and computing output buffers. 

Typical code using the C API:

```C
#include "libswifft/swifft.h"
/* later, inside a function: */
SWIFFT_ALIGN BitSequence input[SWIFFT_INPUT_BLOCK_SIZE]; /* memory-aligned */
SWIFFT_ALIGN BitSequence output[SWIFFT_OUTPUT_BLOCK_SIZE]; /* memory-aligned */
SWIFFT_ALIGN BitSequence compact[SWIFFT_COMPACT_BLOCK_SIZE]; /* optional, memory-aligned */
SWIFFT_ALIGN BitSequence sign[SWIFFT_INPUT_BLOCK_SIZE]; /* optional, memory-aligned */
/* after input (and optionally sign) is populated (not shown here), it is time to compute the hash: */
SWIFFT_Compute(input, output); /* compute the hash of the input into the output */
SWIFFT_Compact(output, compact); /* optionally, compact the hash */
/* or the signed hash can be computed instead */
SWIFFT_Compute(input, sign, output); /* compute the hash of the signed input into the output */
SWIFFT_Compact(output, compact); /* optionally, compact the hash */
```

Buffers must be memory-aligned in order to avoid a segmentation fault when passed to `LibSWIFFT` functions: statically allocated buffers should be aligned using `SWIFFT_ALIGN`, and dynamically allocated buffers should use an alignment of `SWIFFT_ALIGNMENT`, e.g., via `aligned_alloc` function in `stdlib.h`. The transformation functions `SWIFFT_ComputeMultiple{,Signed}*` and `SWIFFT_CompactMultiple` apply operations to multiple blocks. The arithmetic functions `SWIFFT_{Const,}{Set,Add,Sub,Mul}*` provide vectorized and homomorphic operations on an output block, while `SWIFFT_{Const,}{Set,Add,sub,Mul}Multiple*` provide corresponding operations to multiple blocks.

Typical code using the C++ API:

```C
#include "libswifft/swifft.hpp"
using namespace LibSwifft;
/* later, inside a function: */
SwifftInput input; /* auto-memory-aligned */
SwifftOutput output; /* auto-memory-aligned */
SwifftCompact compact; /* optional, auto-memory-aligned */
SwifftInput sign; /* optional, auto-memory-aligned */
/* after input (and optionally sign) is populated (not shown here), it is time to compute the hash: */
SWIFFT_Compute(input.data, output.data); /* compute the hash of the input into the output */
SWIFFT_Compact(output.data, compact.data); /* optionally, compact the hash */
/* or the signed hash can be computed instead */
SWIFFT_Compute(input.data, sign.data, output.data); /* compute the hash of the signed input into the output */
SWIFFT_Compact(output.data, compact.data); /* optionally, compact the hash */
```

Assignment and equality operators are available for `Swifft{Input,Output,Compact}` instances. Arithemtic and arithmetic-assignment operators, corresponding to the arithmetic functions in the C API, are available for `SwifftOutput` instances.

SWIFFT Object APIs are available since `v1.2.0` of `LibSWIFFT` and are recommended:

```C
#include "libswifft/swifft_object.h"
using namespace LibSwifft;
/* initialize object APIs once, possibly inside a function: */
swifft_object_t swifft;
SWIFFT_InitObject(&swifft);
/* later, inside a function: */
SwifftInput input; /* auto-memory-aligned */
SwifftOutput output; /* auto-memory-aligned */
SwifftCompact compact; /* optional, auto-memory-aligned */
/* arithmetic operations are available via swifft.arith, for example: */
swifft.arith.SWIFFT_ConstSet(input.data, 1);
/* hash operations are available via swifft.hash, for example: */
swifft.hash.SWIFFT_Compute(input.data, output.data); /* compute the hash of the input into the output */
swifft.hash.SWIFFT_Compact(output.data, compact.data); /* optionally, compact the hash */
```

Using the object APIs makes it easy to switch their implementation in the future. For the complete SWIFFT object APIs, refer to the documentation or to `src/swifft_object.inl`.

## Building LibSWIFFT

Currently, LibSWIFFT is implemented to be built using GCC. It has been tested on Linux Ubuntu 20.04 LTS using

- `GCC 9.3.0`: GCC is normally installed using `sudo apt-get install gcc g++`.
- `cmake 3.16.3`: CMake is normally installed using `sudo apt-get install cmake`.
- `Catch 2.13.2`: Catch2 is normally installed (as [documented here](https://github.com/catchorg/Catch2/blob/v2.x/docs/cmake-integration.md#installing-catch2-from-git-repository)) using

```sh
git clone https://github.com/catchorg/Catch2.git
cd Catch2
git checkout v2.13.2
cmake -Bbuild -H. -DBUILD_TESTING=OFF
sudo cmake --build build/ --target install
```

The build is also expected to work on older versions of

- `GCC` supporting C++11 as well as avx, avx2, or avx512f
- `cmake` supporting `target_include_directories`
- `Catch2`

Running the following commands:

```sh
mkdir -p build/release
cd build/release
cmake -DCMAKE_BUILD_TYPE=Release ../..
make
```

will build:

- The static library `src/libswifft.a`.
- The shared library `src/libswifft.so`.
- The tests-executable `test/swifft_catch`.

By default, the build will be for the native machine. To build with different machine settings, set `SWIFFT_MACHINE_COMPILE_FLAGS` on the `cmake` command line, for example:

```sh
cmake -DCMAKE_BUILD_TYPE=Release ../.. -DSWIFFT_MACHINE_COMPILE_FLAGS=-march=skylake
```

To build with OpenMP, in particular for parallelizing multiple-block operations, add `-DSWIFFT_ENABLE_OPENMP=on` to the `cmake` command line, for example:

```sh
cmake -DCMAKE_BUILD_TYPE=Release -DSWIFFT_ENABLE_OPENMP=On ../..
```

After building, run the tests-executable from the `build/release` directory:

```sh
./test/swifft_catch
```

If all tests pass, LibSWIFFT is good to go!

For development with LibSWIFFT, use the headers in the `include` directory and either the static or dynamic library.

## Roadmap

Please see [ROADMAP.md](ROADMAP.md).

## Contributing

Please see [CONTRIBUTING.md](CONTRIBUTING.md).
# Design of LibSWIFFT Code

At a high level, LibSWIFFT includes a layer of C code that provides the C API
and a layer of C++ code on top of it that provides the C++ API. New features are
normally added to the C layer and then wrapped by the C++ layer.

The C layer is composed of the main API and the microarchitecture-specific APIs,
discussed in more detail below. In regular use cases, one would normally use the
main API only, which selects the best microarchitecture-specific API for the
platform the library was built for. In more sophisticated use cases, one could
use a microarchitecture-specific APIs directly, e.g., to provide an optimized
feature. The main C and C++ APIs are respectively provided by the `swifft.h` and
`swifft.hpp` headers.

## File and Directory Structure

LibSWIFFT has the following file and directory structure:

| File or Directory              | Description                                           |
| :----------------------------- | :---------------------------------------------------- |
| - `include`                    | root directory of headers                             |
|  - `libswifft`                 | directory of LibSWIFFT headers                        |
|   - `common.h`                 | LibSWIFFT public C common definitions                 |
|   - `swifft.h`                 | LibSWIFFT public C API                                |
|   - `swifft.hpp`               | LibSWIFFT public C++ API                              |
|   - `swifft_avx.h`             | LibSWIFFT public C API for AVX                        |
|   - `swifft_avx2.h`            | LibSWIFFT public C API for AVX2                       |
|   - `swifft_avx512.h`          | LibSWIFFT public C API for AVX512                     |
|   - `swifft_common.h`          | LibSWIFFT public C definitions                        |
|   - `swifft_iset.inl`          | LibSWIFFT public C API expansion for instruction-sets |
|   - `swifft_ver.h`             | LibSWIFFT public C API                                |
| - `src`                        | directory of LibSWIFFT sources                        |
|  - `swifft.c`                  | LibSWIFFT public C implementation                     |
|  - `swifft.inl`                | LibSWIFFT internal C code expansion                   |
|  - `swifft_avx.c`              | LibSWIFFT public C implementation for AVX             |
|  - `swifft_avx2.c`             | LibSWIFFT public C implementation for AVX2            |
|  - `swifft_avx512.c`           | LibSWIFFT public C implementation for AVX512          |
|  - `swifft_impl.inl`           | LibSWIFFT internal C definitions                      |
|  - `swifft_keygen.cpp`         | LibSWIFFT internal C code generation                  |
|  - `swifft_ops.inl`            | LibSWIFFT internal C code expansion                   |
|  - `transpose_8x8_16_sse2.inl` | LibSWIFFT internal C code for matrix transposing      |

## Main API

The main C API has the following organization:

- **FFT functions**: `SWIFFT_fft{,sum}`. These are the two low-level stages in a
  SWIFFT hash computation and are normally not used directly.
- **Transformation functions**: `SWIFFT_{Compute,Compact}`. These transform from
  input to output and from output to compact forms.
- **Arithmetic functions**: `SWIFFT_{,Const}{Set,Add,Sub,Mul}`. These set, add,
  subtract, or multiply given two output forms or one output form and a constant
  value.
- **Functions for multiple blocks**: These are functions with `Multiple` as part
  of their name. They operate on a number of blocks given as a parameter, rather 
  than one block like the corresponding (i.e., without `Multiple`) single-block
  functions.

The main C++ API has the following organization:

- **Wrapper classes**: `Swifft{Input,Output,Compact}`. These wrap input, output,
  and compact buffers with automatic memory-alignment.
- **Wrapper logical operators**: The wrapper classes provide equality and
  inequality operators.
- **Wrapper arithmetic operators**: The wrapper classes provides operators for
  setting, adding, subtracting, and multiplying the current `SwifftOutput`
  instance with another or with a constant value.

## Microarchitecture-Specific APIs

The microarchitecture-specific APIs have the following organization:

- **Similarity to the main C API**: Each microarchitecure-specific function has
  the same name as a corresponding main C API but with an added suffix, the
  same parameter signature, and the same semantics.
- **Name-suffix depending on microarchitecture feature**: There are 3 sets of
  microarchitecture-specific functions corresponding to the 3 suffixes `_AVX`,
  `_AVX2`, and `_AVX512` that respectively provide implementations optimized
  for a microarchitecture supporting AVX, AVX2, and AVX512F instruction-sets.

## Code Conventions

The major code conventions used in LibSWIFFT are:

- SWIFFT-related symbols of the C API start with `SWIFFT\_`.
- Non-SWIFFT-related symbols of the C API start with `LIBSWIFFT\_`
- Symbols in the C++ API are in namespace `LibSwifft`.
- SWIFFT-related functions are thread-safe and constant-time.
- Arguments are assumed non-overlapping in memory.

## Conditional Compilation

LibSWIFFT code uses conditional compilation directives with the following goals:

- Prevent multiple compilations of the same header.
- Compile microarchitecture-specific code only on a supporting platform.
- Compile C++ or OpenMP code only on a supporting compiler.
# Contributing to LibSWIFFT

We welcome contributions to LibSWIFFT of any type, including the following:

- Issue reports
- Patch proposals
- Porting to other operating systems and C/C++ tool chains
- Wrappers of LibSWIFFT in other programming languages
- In- and out-of-code documentation
- Improvements to source code, testing code, build, and performance

Please see [ROADMAP.md](ROADMAP.md) for development directions of interest.

## Issue Reports

Before reporting issues, please check whether the issue is already reported at the [issues page](https://github.com/gvilitechltd/LibSWIFFT/issues). If not, please open a new issue there.

If the issue is a bug report, please provide any information that could help diagnose the issue, such as version of operating system, version of tool chain used to build the library, version of LibSWIFFT used, and the steps for reproducing the problem.

## Patch Proposals

You are encouraged to propose patches via [a pull request](https://github.com/gvilitechltd/LibSWIFFT/pulls), preferrably with reference to an existing issue or one you created for it.

To help us handle your pull request as efficiently as possible, please ensure that:

- each file include a reference to the license as in existing code.
- the style of new code is similar to that of the existing code.
- the commits in the pull request are to the `main` branch only.
- the purpose and method of a code change is documented in the code.
- there are test cases covering code changes that are invoked via the tests-executable.
- all test cases are using Catch2 with appropriate names and categories.
- the changed code builds and tests using the normal procedures.

## Wrappers of LibSWIFFT

Please note that a wrapper of LibSWIFFT in a different programming language should go in a separate project. If you developed such a wrapper, please let us know about it by creating an issue, so we could add it to a published list.

## Acknowledgements

Contributions accepted to LibSWIFFT will be acknowledged in the contributors list.
---
name: Custom issue template
about: Describe this issue template's purpose here.
title: ''
labels: ''
assignees: ''

---


---
name: Bug report
about: Create a report to help us improve
title: "[BUG]"
labels: ''
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior:
1. Minimal code.
2. Build steps.
3. Run steps.
4. Error observed.

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Environment (please complete the following information):**
 - Platform: [e.g. desktop or mobile]
 - OS and version: [e.g. Ubuntu 20.04 LTS]
 - Tool-chain and version [e.g. GCC 9.3.1]
 - LibSWIFFT version [e.g. v1.0.0]

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
Community Guide
===============

Code of Conduct
---------------

Please refer to the :ref:`code_of_conduct`.

Releases
--------

The release process of LibSWIFFT aims to ensure each released version:

* has fully documented API;
* has passed all tests;
* has good code coverage;
* has performance comparable or better than that of the previous version;
* follows semantic versioning rules (see below).

A release version, which appears with the prefix "v" in a repository tag, has
the structure `x.y.z`:

* `x` is the major version, incrementing on a backward-incompatible change,
  which happen rarely and only with strong justification.
* `y` is the minor version, incrementing on backward-compatible enhancement,
  which happens normally and often involves extending the APIs.
* `z` is the patch version, incrementing on a bug fix, which does not change the
  APIs.

For more details, read about `semantic versioning <http://semver.org>`_.

Support
-------

LibSWIFFT support is accessible via the following channels:

* :ref:`Issue reporting <issue_reports>`.
* :ref:`Discussions <discussions>`.

Community Updates
-----------------

Community updates can be tracked by watching the
`LibSWIFFT repository <https://github.com/gvilitechltd/LibSWIFFT>`_.
User Guide
==========

Getting LibSWIFFT
-----------------

LibSWIFFT is available on `GitHub <https://github.com/gvilitechltd/LibSWIFFT>`_.
The code repository can be cloned using `git <https://git-scm.com>`_:

.. code-block:: sh

    git clone https://github.com/gvilitechltd/LibSWIFFT
    cd LibSWIFFT

To switch to a specific version of LibSWIFFT, such as `v1.2.0`, use:

.. code-block:: sh

    git checkout v1.2.0

As an alternative to getting the repository, a specific version of LibSWIFFT can
be obtained from `GitHub <https://github.com/gvilitechltd/LibSWIFFT>`_ by
switching to the tag of this version and downloading the code as an archive.
Then, unpack the archive and change into the unpacked directory, e.g.:

.. code-block:: sh

    wget https://github.com/gvilitechltd/LibSWIFFT/archive/v1.2.0.zip
    unzip v1.2.0.zip
    cd LibSWIFFT-1.2.0

Building LibSWIFFT
------------------

Currently, LibSWIFFT is implemented to be built using GCC. It has been tested on Linux Ubuntu 20.04 LTS using

- `GCC 9.3.0`: GCC is normally installed using `sudo apt-get install gcc g++`.
- `cmake 3.16.3`: CMake is normally installed using `sudo apt-get install cmake`.
- `Catch 2.13.2`: Catch2 is normally installed (as `documented here <https://github.com/catchorg/Catch2/blob/v2.x/docs/cmake-integration.md#installing-catch2-from-git-repository>`_) using

.. code-block:: sh

    git clone https://github.com/catchorg/Catch2.git
    cd Catch2
    git checkout v2.13.2
    cmake -Bbuild -H. -DBUILD_TESTING=OFF
    sudo cmake --build build/ --target install

The build is also expected to work on older versions of

- `GCC` supporting C++11 as well as avx, avx2, or avx512f
- `cmake` supporting `target_include_directories`
- `Catch2`

Running the following commands:

.. code-block:: sh

    mkdir -p build/release
    cd build/release
    cmake -DCMAKE_BUILD_TYPE=Release ../..
    make

will build:

- The static library `src/libswifft.a`.
- The shared library `src/libswifft.so`.
- The tests-executable `test/swifft_catch`.

By default, the build will be for the native machine. To build with different machine settings, set `SWIFFT_MACHINE_COMPILE_FLAGS` on the `cmake` command line, for example:

.. code-block:: sh

    cmake -DCMAKE_BUILD_TYPE=Release ../.. -DSWIFFT_MACHINE_COMPILE_FLAGS=-march=skylake

To build with OpenMP, in particular for parallelizing multiple-block operations, add `-DSWIFFT_ENABLE_OPENMP=on` to the `cmake` command line, for example:

.. code-block:: sh

    cmake -DCMAKE_BUILD_TYPE=Release -DSWIFFT_ENABLE_OPENMP=On ../..

After building, run the tests-executable from the `build/release` directory:

.. code-block:: sh

    ./test/swifft_catch

If all tests pass, LibSWIFFT is good to go!

For development with LibSWIFFT, use the headers in the `include` directory and either the static or dynamic library.

Using LibSWIFFT
---------------

LibSWIFFT is intended to be used by cryptography researchers and by software developers knowledgeable in cryptography programming. LibSWIFFT is most useful in use cases that require provable-security and
speed on short inputs. It may also be interesting in use cases that take advantage of its uncommon homomorphism property. Future versions of LibSWIFFT may target a larger audience.

The main LibSwifft C API is documented in :libswifft:`swifft.h`. The following API variations are available:

- :libswifft:`swifft_avx.h`: Same functions as in :libswifft:`swifft.h` but with an added suffix `_AVX` and implemented using AVX instruction set.
- :libswifft:`swifft_avx2.h`: Same functions as in :libswifft:`swifft.h` but with an added suffix `_AVX2` and implemented using AVX2 instruction set.
- :libswifft:`swifft_avx512.h`: Same functions as in :libswifft:`swifft.h` but with an added suffix `_AVX512` and implemented using AVX512 instruction set.
- :libswifft:`swifft.h`: Selects the implementations using the most advanced instruction set that was built into the library.

The version of LibSWIFFT is provided by the API in :libswifft:`swifft_ver.h`.

The main LibSWIFFT C++ API is documented in :libswifft:`swifft.hpp`.

An extended use of the LibSWIFFT API follows the following steps:

.. |_| unicode:: 0xA0
   :trim:

- **Allocate buffers**: LibSWIFFT defines 3 types of buffers - input, output and compact. The input buffer ':libswifft:`BitSequence` |_| input[:libswifft:`SWIFFT_INPUT_BLOCK_SIZE`]' holds a vector in `Z_2^{2048}` where each element takes 1 bit, the output buffer ':libswifft:`BitSequence` |_| output[:libswifft:`SWIFFT_OUTPUT_BLOCK_SIZE`]' holds a vector in `Z_{257}^{64}` where each element takes 16 bits, and the compact buffer ':libswifft:`BitSequence` |_| compact[:libswifft:`SWIFFT_COMPACT_BLOCK_SIZE`]' holds a value in `Z_{256}^64` taking 64 bytes.
- **Populate input buffers**: An input buffer is populated in preparation for hashing. This is normally done by directly setting the bits in the input buffer. Each bit corresponds to an element of the vector with a value in `{0,1}`.
- **Populate sign buffers**: A sign buffer is an input buffer whose bits are interpreted as sign bits. A 0-valued (resp. 1-valued) bit corresponds to a positive (resp. negative) sign. When an input buffer and a sign buffer are taken together, they define a vector in `{-1,0,1}^{2048}`.
- **Compute output buffers**: The hash of an input buffer, with or without a sign buffer, is computed into an output buffer. This is normally done using :libswifft:`SWIFFT_Compute` or :libswifft:`SWIFFT_ComputeSigned`.
- **Perform arithmetic operations with output buffers**: LibSWIFFT provides several arithemtic (homomorphic) operations involving output buffers whose result is put into an output buffer. The vectors of output buffers may be added, subtracted, or multiplied element-wise. See below for more details.
- **Compact the output buffer**: The hash in the output buffer is compacted into the compact buffer. This is an optional operation, in that the hash in the output buffer may be sufficient for certain applications.

A more restricted use of the LibSWIFFT API involves only the steps of allocating buffers, populating input buffers, and computing output buffers.

Typical code using the C API:

.. code-block:: c

    #include "libswifft/swifft.h"
    /* later, inside a function: */
    SWIFFT_ALIGN BitSequence input[SWIFFT_INPUT_BLOCK_SIZE]; /* memory-aligned */
    SWIFFT_ALIGN BitSequence output[SWIFFT_OUTPUT_BLOCK_SIZE]; /* memory-aligned */
    SWIFFT_ALIGN BitSequence compact[SWIFFT_COMPACT_BLOCK_SIZE]; /* optional, memory-aligned */
    SWIFFT_ALIGN BitSequence sign[SWIFFT_INPUT_BLOCK_SIZE]; /* optional, memory-aligned */
    /* after input (and optionally sign) is populated (not shown here), it is time to compute the hash: */
    SWIFFT_Compute(input, output); /* compute the hash of the input into the output */
    SWIFFT_Compact(output, compact); /* optionally, compact the hash */
    /* or the signed hash can be computed instead */
    SWIFFT_Compute(input, sign, output); /* compute the hash of the signed input into the output */
    SWIFFT_Compact(output, compact); /* optionally, compact the hash */

Buffers must be memory-aligned in order to avoid a segmentation fault when passed to `LibSWIFFT` functions: statically allocated buffers should be aligned using `SWIFFT_ALIGN`, and dynamically allocated buffers should use an alignment of `SWIFFT_ALIGNMENT`, e.g., via `aligned_alloc` function in `stdlib.h`. The transformation functions :libswifft:`SWIFFT_ComputeMultiple`, :libswifft:`SWIFFT_ComputeMultipleSigned` and :libswifft:`SWIFFT_CompactMultiple` apply operations to multiple blocks. The arithmetic functions :libswifft:`SWIFFT_ConstSet`, :libswifft:`SWIFFT_ConstAdd`, :libswifft:`SWIFFT_ConstSub`, :libswifft:`SWIFFT_ConstMul`, :libswifft:`SWIFFT_Set`, :libswifft:`SWIFFT_Add`, :libswifft:`SWIFFT_Sub`, :libswifft:`SWIFFT_Mul` provide vectorized and homomorphic operations on an output block, while :libswifft:`SWIFFT_ConstSetMultiple`, :libswifft:`SWIFFT_ConstAddMultiple`, :libswifft:`SWIFFT_ConstSubMultiple`, :libswifft:`SWIFFT_ConstMulMultiple`, :libswifft:`SWIFFT_SetMultiple`, :libswifft:`SWIFFT_AddMultiple`, :libswifft:`SWIFFT_SubMultiple`, :libswifft:`SWIFFT_Mul` provide corresponding operations to multiple blocks.

Typical code using the C++ API:

.. code-block:: cpp

    #include "libswifft/swifft.hpp"
    using namespace LibSwifft;
    /* later, inside a function: */
    SwifftInput input; /* auto-memory-aligned */
    SwifftOutput output; /* auto-memory-aligned */
    SwifftCompact compact; /* optional, auto-memory-aligned */
    SwifftInput sign; /* optional, auto-memory-aligned */
    /* after input (and optionally sign) is populated (not shown here), it is time to compute the hash: */
    SWIFFT_Compute(input.data, output.data); /* compute the hash of the input into the output */
    SWIFFT_Compact(output.data, compact.data); /* optionally, compact the hash */
    /* or the signed hash can be computed instead */
    SWIFFT_Compute(input.data, sign.data, output.data); /* compute the hash of the signed input into the output */
    SWIFFT_Compact(output.data, compact.data); /* optionally, compact the hash */

Assignment and equality operators are available for :libswifft:`SwifftInput`, :libswifft:`SwifftOutput`, :libswifft:`SwifftCompact` instances. Arithemtic and arithmetic-assignment operators, corresponding to the arithmetic functions in the C API, are available for :libswifft:`SwifftOutput` instances.

Typical code using the recommended SWIFFT object APIs (since `v1.2.0`):

.. code-block:: cpp

    #include "libswifft/swifft_object.h"
    using namespace LibSwifft;
    /* initialize object APIs once, possibly inside a function: */
    swifft_object_t swifft;
    SWIFFT_InitObject(&swifft);
    /* later, inside a function: */
    SwifftInput input; /* auto-memory-aligned */
    SwifftOutput output; /* auto-memory-aligned */
    SwifftCompact compact; /* optional, auto-memory-aligned */
    /* arithmetic operations are available via swifft.arith, for example: */
    swifft.arith.SWIFFT_ConstSet(input.data, 1);
    /* hash operations are available via swifft.hash, for example: */
    swifft.hash.SWIFFT_Compute(input.data, output.data); /* compute the hash of the input into the output */
    swifft.hash.SWIFFT_Compact(output.data, compact.data); /* optionally, compact the hash */

.. _code_of_conduct:

Contributor Covenant Code of Conduct
====================================

Our Pledge
----------

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, religion, or sexual identity
and orientation.

We pledge to act and interact in ways that contribute to an open, welcoming,
diverse, inclusive, and healthy community.

Our Standards
-------------

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

Enforcement Responsibilities
----------------------------

Community leaders are responsible for clarifying and enforcing our standards of
acceptable behavior and will take appropriate and fair corrective action in
response to any behavior that they deem inappropriate, threatening, offensive,
or harmful.

Community leaders have the right and responsibility to remove, edit, or reject
comments, commits, code, wiki edits, issues, and other contributions that are
not aligned to this Code of Conduct, and will communicate reasons for moderation
decisions when appropriate.

Scope
-----

This Code of Conduct applies within all community spaces, and also applies when
an individual is officially representing the community in public spaces.
Examples of representing our community include using an official e-mail address,
posting via an official social media account, or acting as an appointed
representative at an online or offline event.

Enforcement
-----------

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported to the community leaders responsible for enforcement at
`gvilitechltd@gmail.com <mailto:gvilitechltd@gmail.com>`_.
All complaints will be reviewed and investigated promptly and fairly.

All community leaders are obligated to respect the privacy and security of the
reporter of any incident.

Enforcement Guidelines
----------------------

Community leaders will follow these Community Impact Guidelines in determining
the consequences for any action they deem in violation of this Code of Conduct:

Correction
^^^^^^^^^^

**Community Impact**: Use of inappropriate language or other behavior deemed
unprofessional or unwelcome in the community.

**Consequence**: A private, written warning from community leaders, providing
clarity around the nature of the violation and an explanation of why the
behavior was inappropriate. A public apology may be requested.

Warning
^^^^^^^

**Community Impact**: A violation through a single incident or series
of actions.

**Consequence**: A warning with consequences for continued behavior. No
interaction with the people involved, including unsolicited interaction with
those enforcing the Code of Conduct, for a specified period of time. This
includes avoiding interactions in community spaces as well as external channels
like social media. Violating these terms may lead to a temporary or
permanent ban.

Temporary Ban
^^^^^^^^^^^^^

**Community Impact**: A serious violation of community standards, including
sustained inappropriate behavior.

**Consequence**: A temporary ban from any sort of interaction or public
communication with the community for a specified period of time. No public or
private interaction with the people involved, including unsolicited interaction
with those enforcing the Code of Conduct, is allowed during this period.
Violating these terms may lead to a permanent ban.

Permanent Ban
^^^^^^^^^^^^^

**Community Impact**: Demonstrating a pattern of violation of community
standards, including sustained inappropriate behavior,  harassment of an
individual, or aggression toward or disparagement of classes of individuals.

**Consequence**: A permanent ban from any sort of public interaction within
the community.

Attribution
-----------

This Code of Conduct is adapted from the
`Contributor Covenant <http://contributor-covenant.org>`_
version `2.0 <https://www.contributor-covenant.org/version/2/0/code_of_conduct.html>`_.
Design of LibSWIFFT Code
========================

At a high level, LibSWIFFT includes a layer of C code that provides the C API
and a layer of C++ code on top of it that provides the C++ API. New features are
normally added to the C layer and then wrapped by the C++ layer.

The C layer is composed of the main API and the microarchitecture-specific APIs,
discussed in more detail below. In regular use cases, one would normally use the
main API only, which selects the best microarchitecture-specific API for the
platform the library was built for. In more sophisticated use cases, one could
use a microarchitecture-specific APIs directly, e.g., to provide an optimized
feature. The main C and C++ APIs are respectively provided by the `swifft.h` and
`swifft.hpp` headers.

File and Directory Structure
----------------------------

LibSWIFFT has the following file and directory structure:

.. list-table::
   :widths: 50 50
   :header-rows: 1

   * - File or Directory
     - Description
   * - include
     - root directory of headers
   * - . libswifft
     - directory of LibSWIFFT headers
   * - . . :libswifft:`common.h`
     - LibSWIFFT public C common definitions
   * - . . :libswifft:`swifft.h`
     - LibSWIFFT public C API
   * - . . :libswifft:`swifft.hpp`
     - LibSWIFFT public C++ API
   * - . . :libswifft:`swifft_avx.h`
     - LibSWIFFT public C API for AVX
   * - . . :libswifft:`swifft_avx2.h`
     - LibSWIFFT public C API for AVX2
   * - . . :libswifft:`swifft_avx512.h`
     - LibSWIFFT public C API for AVX512
   * - . . :libswifft:`swifft_common.h`
     - LibSWIFFT public C definitions
   * - . . :libswifft:`swifft_iset.inl`
     - LibSWIFFT public C API expansion for instruction-sets
   * - . . :libswifft:`swifft_ver.h`
     - LibSWIFFT public C API
   * -  src
     - directory of LibSWIFFT sources
   * - . :libswifft:`swifft.c`
     - LibSWIFFT public C implementation
   * - . :libswifft:`swifft.inl`
     - LibSWIFFT internal C code expansion
   * - . :libswifft:`swifft_avx.c`
     - LibSWIFFT public C implementation for AVX
   * - . :libswifft:`swifft_avx2.c`
     - LibSWIFFT public C implementation for AVX2
   * - . :libswifft:`swifft_avx512.c`
     - LibSWIFFT public C implementation for AVX512
   * - . :libswifft:`swifft_impl.inl`
     - LibSWIFFT internal C definitions
   * - . :libswifft:`swifft_keygen.cpp`
     - LibSWIFFT internal C code generation
   * - . :libswifft:`swifft_ops.inl`
     - LibSWIFFT internal C code expansion
   * - . :libswifft:`transpose_8x8_16_sse2.inl`
     - LibSWIFFT internal C code for matrix transposing

Main API
--------

The main C API has the following organization:

- **FFT functions**: :libswifft:`SWIFFT_fft`, :libswifft:`SWIFFT_fftsum`.
  These are the two low-level stages in a SWIFFT hash computation and are
  normally not used directly.
- **Transformation functions**: :libswifft:`SWIFFT_Compute`,
  :libswifft:`SWIFFT_Compact`. These transform from input to output and from
  output to compact forms.
- **Arithmetic functions**: :libswifft:`SWIFFT_Set`, :libswifft:`SWIFFT_Add`,
  :libswifft:`SWIFFT_Sub`, :libswifft:`SWIFFT_Mul`,
  :libswifft:`SWIFFT_ConstSet`, :libswifft:`SWIFFT_ConstAdd`,
  :libswifft:`SWIFFT_ConstSub`, :libswifft:`SWIFFT_ConstMul`. These set, add,
  subtract, or multiply given two output forms or one output form and a constant
  value.
- **Functions for multiple blocks**: These are functions with `Multiple` as part
  of their name. They operate on a number of blocks given as a parameter, rather
  than one block like the corresponding (i.e., without `Multiple`) single-block
  functions.

The main C++ API has the following organization:

- **Wrapper classes**: :libswifft:`SwifftInput`, :libswifft:`SwifftOutput`,
  :libswifft:`SwifftCompact`. These wrap input, output, and compact buffers with
  automatic memory-alignment.
- **Wrapper logical operators**: The wrapper classes provide equality and
  inequality operators.
- **Wrapper arithmetic operators**: The wrapper classes provides operators for
  setting, adding, subtracting, and multiplying the current `SwifftOutput`
  instance with another or with a constant value.

Microarchitecture-Specific APIs
-------------------------------

The microarchitecture-specific APIs have the following organization:

- **Similarity to the main C API**: Each microarchitecure-specific function has
  the same name as a corresponding main C API but with an added suffix, the
  same parameter signature, and the same semantics.
- **Name-suffix depending on microarchitecture feature**: There are 3 sets of
  microarchitecture-specific functions corresponding to the 3 suffixes `_AVX`,
  `_AVX2`, and `_AVX512` that respectively provide implementations optimized
  for a microarchitecture supporting AVX, AVX2, and AVX512F instruction-sets.

Code Conventions
----------------

Major code conventions used in LibSWIFFT are:

- SWIFFT-related symbols of the C API start with `SWIFFT\_`.
- Non-SWIFFT-related symbols of the C API start with `LIBSWIFFT\_`
- Symbols in the C++ API are in namespace `LibSwifft`.
- SWIFFT-related functions are thread-safe and constant-time.
- Arguments are assumed non-overlapping in memory.

Conditional Compilation
-----------------------

Major goals of LibSWIFFT code in using conditional compilation are:

- Prevent a header for being compiled multiple times.
- Compile microarchitecture-specific code only on a supporting platform.
- Compile C++ or OpenMP code only on a supporting compiler.
Rationale for Making LibSWIFFT
==============================

LibSWIFFT is based on academic work from 2007 described further below.

Why now, in early 2021?
-----------------------

In late 2017, NIST has started a `process for standardizing post-quantum cryptography <https://csrc.nist.gov/Projects/post-quantum-cryptography/post-quantum-cryptography-standardization>`_, suggesting that it believes it may not be too long before a practical quantum-computer that threatens critical security standards (including Internet ones) based on classical cryptography will become a reality. About two years later, Google announced it had achieved `quantum supremacy <https://en.wikipedia.org/wiki/Quantum_supremacy>`_, by completing in 200 seconds a task they claimed would have taken a classical supercomputer about 10,000 years to complete. Though IBM, maker of the most powerful supercomputer at the time, disputed this claim and asserted the supercomputer would take only about 2.5 days for the task, it is clear quantum computing technology is advancing quickly. Consequently, post-quantum cryptography is becoming more relevant today and perhaps even urgent to develop.

Why another implementation of SWIFFT?
-------------------------------------

LibSWIFFT is a reliable building block for fast and scalable cryptographic protocols. It is simple to use and maintain, has clean APIs, is well documented and tested, and is at least as fast as other implementations of SWIFFT and often faster. Other implementations of SWIFFT are:

- The `8-bit <https://github.com/anon1985/Swifft-avx2-8>`_ and `16-bit <https://github.com/anon1985/K2SN-MSS/tree/master/swifft16>`_ AVX2 implementations for `K2SN-MSS <https://eprint.iacr.org/2019/442.pdf>`_. Both are arguably not as easy to use nor as well documented and tested as LibSWIFFT. The former is slower and uses less memory while the latter is about as fast as LibSWIFFT for AVX2 yet does not support AVX512.
- `The original implementation <https://github.com/micciancio/SWIFFT>`_ written in 2007. It is minimal non-production code. `The AVX2 implementation for K2SN-MSS <https://eprint.iacr.org/2019/442.pdf>`_ is reported to be 25% faster.

An invocation of the tests-executable of LibSWIFFT running single-threaded using AVX2 on an Intel Skylake microarchitecture (Intel(R) Core(TM) i7-10875H CPU @ 2.30GHz):

.. code-block:: sh

    $ ./test/swifft_catch "swifft takes at most 2000 cycles per call"
    Filters: swifft takes at most 2000 cycles per call
    running 1*10000000 rounds: cycles/rounds=1097.94 cycles/byte=4.28882 Giga-cycles/sec=2.30399 MB/sec=512.322 cycles/rdtsc=16

demonstrates that LibSWIFFT is quite fast on short inputs (here, 256 bytes), often used in practical zero-knowledge proofs and post-quantum digital signatures. This is more than an order of magnitude faster than the `originally reported <https://www.alonrosen.net/PAPERS/lattices/swifft.pdf>`_ 40MB/sec on a 3.2 GHz Intel Pentium 4. This is also faster than `K2SN-MSS's binary 16-bit SWIFFT function implementation <https://github.com/gvilitechltd/K2SN-MSS/tree/swifftperf>`_ (for an input of 128 bytes), which is the fastest one in the K2SN-MSS implementation, for the same executaion settings, i.e. running single-threaded using AVX2 on an Intel Skylake microarchitecture (Intel(R) Core(TM) i7-10875H CPU @ 2.30GHz):

.. code-block:: sh

    $ ./tester
    1000000 SWIFFT16 rounds: cycles/round=737.363098 cycles/byte=5.760649

It also compares well with modern hash functions:

- `Blake3 <https://github.com/BLAKE3-team/BLAKE3>`_ - cryptographic hash function achieving about `3-to-4 cycles/byte using AVX512 on short inputs <https://github.com/BLAKE3-team/BLAKE3-specs/blob/master/blake3.pdf>`_ and are non-homomorphic nor facilitating proofs of knowledge of a preimage.
- `Seahash <https://docs.rs/seahash/4.0.1/seahash/index.html>`_ - a hash function achieving ~0.24 cycles/byte but is non-cryptographic.
.. LibSWIFFT documentation master file, created by
   sphinx-quickstart on Sat Jan 30 16:08:46 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. toctree::
   :maxdepth: 2

.. image:: ../assets/LibSWIFFT-logo.png

`LibSWIFFT v1.2.0 <https://doi.org/10.21105/joss.03040>`_ |josspaper|

.. |josspaper| image:: https://joss.theoj.org/papers/10.21105/joss.03040/status.svg
   :target: https://doi.org/10.21105/joss.03040

LibSWIFFT: Fast Secure Homomorphic Hashing
==========================================

LibSWIFFT is a production-ready C/C++ library providing SWIFFT, one of the
fastest available secure hash functions, which is also collision-resistant as
well as facilitates zero-knowledge proofs of knowledge of a preimage (ZKPoKP)
and post-quantum digital signatures.

Secure Homomorphic Hashing
--------------------------

LibSWIFFT provides a function to hash an input vector of 2048 values in [0..1],
given in a buffer of 256 bytes (packed 8-bits), to an output vector of 64 values
in [0..256], given in a buffer of 64 words (padded 16-bits). It also provides a
function to hash a vector of 2048 values in [-1..1], given in an input and a
sign vector, each of 2048 bits. These are implementations of the SWIFFT secure
hash function homomorphic (i.e., the hash of the sum of vectors is the sum of
their hashes) in modulu-257 arithmetic.

Fast
----

LibSWIFFT is fast. Here is the output from a test-case of LibSwifft v1.0.0 that
demonstrates hashing (for an input of 256 bytes) running at a rate of less than
5 cycles-per-byte single-threaded using AVX2 on an Intel Skylake
microarchitecture (Intel(R) Core(TM) i7-10875H CPU @ 2.30GHz):

.. code-block:: sh

    $ ./test/swifft_catch "swifft takes at most 2000 cycles per call"
    Filters: swifft takes at most 2000 cycles per call
    running 1*10000000 rounds: cycles/rounds=1097.94 cycles/byte=4.28882 Giga-cycles/sec=2.30399 MB/sec=512.322 cycles/rdtsc=16

This is faster than the rate of the recent
`K2SN-MSS's binary 16-bit SWIFFT function implementation <https://github.com/gvilitechltd/K2SN-MSS/tree/swifftperf>`_
(for an input of 128 bytes) on the same execution settings, i.e. single-threaded
using AVX2 on an Intel Skylake microarchitecture (Intel(R) Core(TM) i7-10875H
CPU @ 2.30GHz):

.. code-block:: sh

   $ ./tester
    1000000 SWIFFT16 rounds: cycles/round=737.363098 cycles/byte=5.760649

LibSWIFFT supports microarchitectures with AVX, AVX2, and AVX512 instruction
sets to take advantage of available hardware acceleration.

Production-Ready
----------------

LibSWIFFT v1.2.0 ships with over 30 test-cases including millions of checks that
provide excellent coverage of the library's code and the mathematical properties
of the SWIFFT function it implements:

.. code-block:: sh

    $ ./test/swifft_catch

    ===============================================================================
    All tests passed (6648384 assertions in 32 test cases)

LibSWIFFT also ships with detailed documentation, including references for its
:libswifft:`C API <swifft.h>` and :libswifft:`C++ API <LibSwifft>`. It has a
robust and portable build system based on `cmake <http://cmake.org>`_. All
include files are conveniently organized under a dedicated
`libswifft <api/files.html>`_ directory. Finally,
it requires few dependencies at build-time and run-time.

Easy to use C and C++ APIs
--------------------------

LibSWIFFT is also easy to use. Minimal code using the
:libswifft:`C API <swifft.h>` looks like this:

.. code-block:: c

    #include "libswifft/swifft.h"
    /* later, inside a function: */
    SWIFFT_ALIGN BitSequence input[SWIFFT_INPUT_BLOCK_SIZE];
    SWIFFT_ALIGN BitSequence output[SWIFFT_OUTPUT_BLOCK_SIZE];
    /* after input is populated (not shown here), it is time to compute the hash: */
    SWIFFT_Compute(input, output); /* compute the hash of the input into the output */

The :libswifft:`SWIFFT_ALIGN` macro is needed for the required memory-alignment
of buffers passed to LibSWIFFT functions. The equivalent code using the
:libswifft:`C++ API <LibSwifft>` looks like this:

.. code-block:: cpp

    #include "libswifft/swifft.hpp"
    using namespace LibSwifft;
    /* later, inside a function: */
    SwifftInput input;
    SwifftOutput output;
    /* after input is populated (not shown here), it is time to compute the hash: */
    SWIFFT_Compute(input.data, output.data); /* compute the hash of the input into the output */

:libswifft:`SwifftInput` and :libswifft:`SwifftOutput` are auto-memory-aligned.

The code using the :libswifft:`C API <swifft.h>` for signed hashing in
compact-form looks like this:

.. code-block:: c

    #include "libswifft/swifft.h"
    /* later, inside a function: */
    SWIFFT_ALIGN BitSequence input[SWIFFT_INPUT_BLOCK_SIZE];
    SWIFFT_ALIGN BitSequence output[SWIFFT_OUTPUT_BLOCK_SIZE];
    SWIFFT_ALIGN BitSequence compact[SWIFFT_COMPACT_BLOCK_SIZE];
    SWIFFT_ALIGN BitSequence sign[SWIFFT_INPUT_BLOCK_SIZE];
    /* after input and sign are populated (not shown here), it is time to compute the hash: */
    SWIFFT_Compute(input, sign, output); /* compute the hash of the signed input into the output */
    SWIFFT_Compact(output, compact); /* compact the hash */

And the equivalent code using the :libswifft:`C++ API <LibSwifft>` looks like
this:

.. code-block:: cpp

    #include "libswifft/swifft.hpp"
    using namespace LibSwifft;
    /* later, inside a function: */
    SwifftInput input;
    SwifftOutput output;
    SwifftCompact compact;
    SwifftInput sign;
    /* after input and sign are populated (not shown here), it is time to compute the hash: */
    SWIFFT_Compute(input.data, sign.data, output.data); /* compute the hash of the signed input into the output */
    SWIFFT_Compact(output.data, compact.data);

For future flexibility in implementation, it is recommended to invoke SWIFFT
functions via the :libswifft:`SWIFFT object APIs <swifft_object.h>` (since
`v1.2.0`):

.. code-block:: cpp

    #include "libswifft/swifft_object.h"
    using namespace LibSwifft;
    /* initialize object APIs once, possibly inside a function: */
    swifft_object_t swifft;
    SWIFFT_InitObject(&swifft);
    /* later, inside a function: */
    SwifftInput input; /* auto-memory-aligned */
    SwifftOutput output; /* auto-memory-aligned */
    SwifftCompact compact; /* optional, auto-memory-aligned */
    /* arithmetic operations are available via swifft.arith, for example: */
    swifft.arith.SWIFFT_ConstSet(input.data, 1);
    /* hash operations are available via swifft.hash, for example: */
    swifft.hash.SWIFFT_Compute(input.data, output.data); /* compute the hash of the input into the output */
    swifft.hash.SWIFFT_Compact(output.data, compact.data); /* optionally, compact the hash */

Rationale
=========

Here the rationale for making LibSWIFFT is explained.

.. toctree::
   :maxdepth: 2

   rationale

About SWIFFT and LibSWIFFT
==========================

This part describes SWIFFT and LibSWIFFT.

.. toctree::
   :maxdepth: 2

   about

User Guide
==========

The User Guide explains how to get started with LibSWIFFT and how to make the
most of it.

.. toctree::
   :maxdepth: 2

   user-guide

Community Guide
===============

The Community Guide explains how to participate in LibSWIFFT's community.

.. toctree::
   :maxdepth: 2

   community-guide
   code-design

API Guide
=========

The API Guide leads to the LibSWIFFT API reference.

.. toctree::
   :maxdepth: 2

   api-guide

Contributor Guide
=================

The Contributor Guide explains how to contribute to LibSWIFFT.

.. toctree::
   :maxdepth: 2

   contributor-guide
   code-of-conduct
API Guide
=========

Here are useful entry points to the LibSWIFFT API reference:

* `Main Page <api/index.html>`_
* `Related Pages <api/pages.html>`_
* `Files <api/files.html>`_
* :libswifft:`C API <swifft.h>`
* :libswifft:`C++ API <swifft.hpp>`
About SWIFFT and LibSWIFFT
==========================

On SWIFFT
---------

`SWIFFT <https://en.wikipedia.org/wiki/SWIFFT>`_ is a family of homomorphic hash functions provided as part of a candidate submission to the `NIST hash function competition <https://en.wikipedia.org/wiki/NIST_hash_function_competition>`_ in 2008 by Vadim Lyubashevsky, Daniele Micciancio, Chris Peikert, and Alon Rosen. The family has an interesting set of characteristics:

1. **Provably collision-resistant**: it is provably computationally hard to find a binary-valued `x` given its image under `f` that is chosen randomly.
2. **Universal**: for any `x_1,x_2` in the domain the probability that `f(x_1) = f(x_2)` over the choice of `f` in the family is the inverse of the size of the range.
3. **Regular**: if `x` is distributed uniformly then `f(x)` is distributed uniformly for any `f` in the family.
4. **Randomness extractor**: even if `x` is not distributed uniformly, `f(x)` is distributed uniformly over the choice of `f` in the family.
5. **Homomorphic**: each `f` in the family is homomorphic over the domain.
6. **Facilitates ZKPoKPs**: the resulting hash value can be readily used in ZKPoKPs.
7. **Facilitates post-quantum digital signatures**: its homomorphism property enables short post-quantum hash-based digital signature schemes, such as `K2SN-MSS <https://eprint.iacr.org/2019/442.pdf>`_.
8. **Constant-time**: any `f` in the family is free of data-dependent branching and therefore facilitates avoidance of timing side-channel attacks.

Nevertheless, the family is not pseudorandom, due to the homomorphism property, nor claimed to behave like a random oracle.

On LibSWIFFT
------------

LibSWIFFT was implemented with reference to the `SWIFFTX submission to NIST <https://csrc.nist.gov/projects/hash-functions/sha-3-project>`_ and provides the same SWIFFT hash function that is part of the submission. High speed is achieved using various code optimization techniques, including SIMD instructions that are very natural for the implementation of the SWIFFT function. Compared to the SWIFFT code in the submission, LibSWIFFT adds the following:

1. Automatic library initialization using build-time generation of internal tables.
2. Convenient APIs, including for homomorphic operations, for computing SWIFFT on short inputs.
3. Support for input vectors of either binary-valued (in {0,1}) or trinary-valued (in {-1,0,1}) elements.
4. Bug fixes with respect to the reference submission, in particular related to the homomorphism property.
5. Performance improvements compared to the reference submission.
6. Support for newer CPU instruction sets: AVX, AVX2, and AVX512.
7. Over 20 test-cases providing excellent coverage of the APIs and the mathematical properties of SWIFFT.

Formally, LibSWIFFT provides a single hash function that maps from an input domain `Z_2^{2048}` (taking 256B) to an output domain `Z_{257}^{64}` (taking 128B, at 2B per element) and then to a compact domain `Z_{256}^{65}` (taking 65B). The computation of the first map is done over `Z_{257}`. The homomorphism property applies to the input and output domains, but not to the compact domain, and is revealed when the binary-valued input domain is naturally embedded in `Z_{257}^{2048}`. Generally, it is computationally hard to find a binary-valued pre-image given an output computed as the sum of `N` outputs corresponding to known binary-valued pre-images. On the other hand, it is easy to find a small-valued pre-image (over `Z_{257}^{2048}`) when `N` is small, since it is simply the sum of the known pre-images due to the homomorphism property.
Contributor Guide
=================

Contributing to LibSWIFFT
-------------------------

We welcome contributions to LibSWIFFT of any type, including the following:

- :ref:`Issue reports <issue_reports>`
- :ref:`Discussions <discussions>`
- :ref:`Patch proposals <patch_proposals>`
- :ref:`Porting <porting>` to other operating systems and C/C++ tool chains
- :ref:`Wrappers <wrappers>` of LibSWIFFT in other programming languages
- :ref:`documentation <documentation>` within and outside the code
- :ref:`Improvements <improvements>` to source code, testing code, build system, and performance

Contributors are encouraged to watch the
`LibSWIFFT repository <https://gihub.com/gvilitechltd/LibSWIFFT>`_ for updates.

.. _roadmap:

Roadmap of LibSWIFFT
--------------------

The list below details, in no particular order, a number of directions for future development of LibSWIFFT. We welcome contributions in each of these as well as suggestions for additional directions.

- Support for higher-modulus variations of SWIFFT.
- Improved build modularity: multi-versioned build installation, packaging for availability via `find_package` in cmake.
- Build-support for additional platforms, operating systems and toolchains.
- Improved test coverage: numerical edge cases.
- Automatic runtime selection of the best SWIFFT implementation for the native platform.
- Support for parallel processing using OpenMP.
- SWIFFT based hashing for long inputs.
- GPU kernels for SWIFFT functions.

.. _out_of_scope:

Out of Scope for LibSWIFFT
--------------------------

The list below gives a few examples of relevant directions that are out of scope for LibSWIFFT and should go into separate projects.

- Wrappers of LibSWIFFT in other programming languages.
- Integrations of LibSWIFFT with established cryptography libraries.
- Implementations of proofs of knowledge or digital signatures based on LibSWIFFT.

.. _issue_reports:

Issue Reports
-------------

Before reporting issues, please check whether the issue is already reported at the `issues page <https://github.com/gvilitechltd/LibSWIFFT/issues>`_. If not, please open a new issue there.

If the issue is a bug report, please provide any information that could help diagnose the issue, such as version of operating system, version of tool chain used to build the library, version of LibSWIFFT used, and the steps for reproducing the problem.

.. _discussions:

Discussions
-----------

Before starting new discussions, please consider whether an issue report would be more appropriate. If a discussion is appropriate, check whether you can join an existing discussion that has already been created at the `discussion page <https://github.com/gvilitechltd/LibSWIFFT/discussions>`_. If you cannot find such a discussion, feel free to open a new one there.

Please provide a clear topic for your new discussion and elaborate in the body of the first message. The goal is to keep each discussion focused on one topic.

.. _patch_proposals:

Patch Proposals
---------------

You are encouraged to propose patches via `a pull request <https://github.com/gvilitechltd/LibSWIFFT/pulls>`_, preferrably with reference to an existing issue or one you created for it.

To help us handle your pull request as efficiently as possible, please ensure that:

- each file include a reference to the license as in existing code.
- the style of new code is similar to that of the existing code.
- the commits in the pull request are to the `main` branch only.
- the purpose and method of a code change is documented in the code.
- there are test cases covering code changes that are invoked via the tests-executable.
- all test cases are using Catch2 with appropriate names and categories.
- the changed code builds and tests using the normal procedures.

.. _porting:

Porting
-------

Please note that porting to other operating systems and C/C++ tool chains is expected to enhance the existing `cmake <http://cmake.org>`_-based build system. In case of running into an issue with meeting this expectation for some porting effort, please open a discussion about it.

.. _wrappers:

Wrappers of LibSWIFFT
---------------------

Please note that a wrapper of LibSWIFFT in a different programming language should go in a separate project. If you developed such a wrapper, please let us know about it by creating an issue, so we could add it to a published list.

.. _documentation:

Documentation
-------------

.. _improvements:

Contributions to documentation of various kinds, within and outside of the code, are welcomed, for example:

* Documentation of the build system.
* Explanation of mathematics, algorithms, techniques, and logic.
* Documentation of correctness and performance tests.
* Description of architecture, design, and conventions.
* Discussion or examples of secure and insecure uses of the APIs.
* Guidance for specific types of code contributions, such as multi-threading features.

Improvements
------------

Contributions that improve the source code, testing code, build system, and performance of LibSWIFFT are welcomed. Please refer to the :ref:`roadmap <roadmap>` for possible directions and know what is considered :ref:`out of scope <out_of_scope>` for LibSWIFFT.

Acknowledgements
----------------

Contributions accepted to LibSWIFFT will be acknowledged in the contributors list.

