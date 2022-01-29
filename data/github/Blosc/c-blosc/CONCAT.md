# Blosc: A blocking, shuffling and lossless compression library
| Author | Contact | URL |
|--------|---------|-----|
| Blosc Development Team | blosc@blosc.org | http://www.blosc.org | 

| Gitter | GH Actions | NumFOCUS | Code of Conduct |
|--------|------------|----------|-----------------|
| [![Gitter](https://badges.gitter.im/Blosc/c-blosc.svg)](https://gitter.im/Blosc/c-blosc?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge) | [![CI CMake](https://github.com/Blosc/c-blosc/workflows/CI%20CMake/badge.svg)](https://github.com/Blosc/c-blosc/actions?query=workflow%3A%22CI+CMake%22) | [![Powered by NumFOCUS](https://img.shields.io/badge/powered%20by-NumFOCUS-orange.svg?style=flat&colorA=E1523D&colorB=007D8A)](http://numfocus.org) | [![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4.svg)](code_of_conduct.md) |

## What is it?

Blosc is a high performance compressor optimized for binary data.
It has been designed to transmit data to the processor cache faster
than the traditional, non-compressed, direct memory fetch approach via
a memcpy() OS call.  Blosc is the first compressor (that I'm aware of)
that is meant not only to reduce the size of large datasets on-disk or
in-memory, but also to accelerate memory-bound computations.

It uses the [blocking technique](http://blosc.org/docs/StarvingCPUs-CISE-2010.pdf)
so as to reduce activity in the memory bus as much as possible. In short, this
technique works by dividing datasets in blocks that are small enough
to fit in caches of modern processors and perform compression /
decompression there.  It also leverages, if available, SIMD
instructions (SSE2, AVX2) and multi-threading capabilities of CPUs, in
order to accelerate the compression / decompression process to a
maximum.

See some [benchmarks](http://blosc.org/pages/synthetic-benchmarks/) about Blosc performance.

Blosc is distributed using the BSD license, see LICENSES/BLOSC.txt for
details.

## Meta-compression and other differences over existing compressors

C-Blosc is not like other compressors: it should rather be called a
meta-compressor.  This is so because it can use different compressors
and filters (programs that generally improve compression ratio).  At
any rate, it can also be called a compressor because it happens that
it already comes with several compressor and filters, so it can
actually work like a regular codec.

Currently C-Blosc comes with support of BloscLZ, a compressor heavily
based on FastLZ (http://fastlz.org/), LZ4 and LZ4HC
(https://github.com/Cyan4973/lz4), Snappy
(https://github.com/google/snappy), Zlib (http://www.zlib.net/) and
Zstd (http://www.zstd.net).

C-Blosc also comes with highly optimized (they can use
SSE2 or AVX2 instructions, if available) shuffle and bitshuffle filters
(for info on how and why shuffling works [see here](https://www.slideshare.net/PyData/blosc-py-data-2014/17?src=clipshare)).
However, additional compressors or filters may be added in the future.

Blosc is in charge of coordinating the different compressor and
filters so that they can leverage the 
[blocking technique](http://blosc.org/docs/StarvingCPUs-CISE-2010.pdf)
as well as multi-threaded execution (if several cores are
available) automatically. That makes that every codec and filter
will work at very high speeds, even if it was not initially designed
for doing blocking or multi-threading.

Finally, C-Blosc is specially suited to deal with binary data because
it can take advantage of the type size meta-information for improved
compression ratio by using the integrated shuffle and bitshuffle filters.

When taken together, all these features set Blosc apart from other
compression libraries.

## Compiling the Blosc library

Blosc can be built, tested and installed using CMake_.
The following procedure describes the "out of source" build.

```console

  $ cd c-blosc
  $ mkdir build
  $ cd build
```

Now run CMake configuration and optionally specify the installation
directory (e.g. '/usr' or '/usr/local'):

```console

  $ cmake -DCMAKE_INSTALL_PREFIX=your_install_prefix_directory ..
```

CMake allows to configure Blosc in many different ways, like preferring
internal or external sources for compressors or enabling/disabling
them.  Please note that configuration can also be performed using UI
tools provided by [CMake](http://www.cmake.org) (ccmake or cmake-gui):

```console

  $ ccmake ..      # run a curses-based interface
  $ cmake-gui ..   # run a graphical interface
```

Build, test and install Blosc:


```console

  $ cmake --build .
  $ ctest
  $ cmake --build . --target install
```

The static and dynamic version of the Blosc library, together with
header files, will be installed into the specified
CMAKE_INSTALL_PREFIX.

### Codec support with CMake

C-Blosc comes with full sources for LZ4, LZ4HC, Snappy, Zlib and Zstd
and in general, you should not worry about not having (or CMake
not finding) the libraries in your system because by default the
included sources will be automatically compiled and included in the
C-Blosc library. This effectively means that you can be confident in
having a complete support for all the codecs in all the Blosc deployments
(unless you are explicitly excluding support for some of them).

But in case you want to force Blosc to use external codec libraries instead of
the included sources, you can do that:

``` console

  $ cmake -DPREFER_EXTERNAL_ZSTD=ON ..
```

You can also disable support for some compression libraries:


```console

  $ cmake -DDEACTIVATE_SNAPPY=ON ..  # in case you don't have a C++ compiler
```
 
## Examples

In the [examples/ directory](https://github.com/Blosc/c-blosc/tree/master/examples)
you can find hints on how to use Blosc inside your app.

## Supported platforms

Blosc is meant to support all platforms where a C89 compliant C
compiler can be found.  The ones that are mostly tested are Intel
(Linux, Mac OSX and Windows) and ARM (Linux), but exotic ones as IBM
Blue Gene Q embedded "A2" processor are reported to work too.

### Mac OSX troubleshooting

If you run into compilation troubles when using Mac OSX, please make
sure that you have installed the command line developer tools.  You
can always install them with:

```console

  $ xcode-select --install
```

## Wrapper for Python

Blosc has an official wrapper for Python.  See:

https://github.com/Blosc/python-blosc

## Command line interface and serialization format for Blosc

Blosc can be used from command line by using Bloscpack.  See:

https://github.com/Blosc/bloscpack

## Filter for HDF5

For those who want to use Blosc as a filter in the HDF5 library,
there is a sample implementation in the hdf5-blosc project in:

https://github.com/Blosc/hdf5-blosc

## Mailing list

There is an official mailing list for Blosc at:

blosc@googlegroups.com
http://groups.google.es/group/blosc

## Acknowledgments

See THANKS.rst.


----

  **Enjoy data!**
# Contributing to C-Blosc
We want to make contributing to this project as easy and transparent as
possible.

## Our Development Process
New versions are being developed in the "master" branch,
or in their own feature branch.
When they are deemed ready for a release, they are merged back into "master"
again.

So all contributions must stage first through "master"
or their own feature branch.

## Pull Requests
We actively welcome your pull requests.

1. Fork the repo and create your branch from `master`.
2. If you've added code that should be tested, add tests.
3. If you've changed APIs, update the documentation.
4. Ensure the test suite passes.
5. Make sure your code does not issue new compiler warnings.

## Issues
We use GitHub issues to track public bugs. Please ensure your description is
clear and has sufficient instructions to be able to reproduce the issue.

## Coding Style  
* 2 spaces for indentation rather than tabs.

## License
By contributing to C-Blosc, you agree that your contributions will be licensed
under both the [LICENSE](LICENSES/BLOSC.txt) file of the project.
# Code of Conduct

The Blosc community has adopted a Code of Conduct that we expect project participants to adhere to.
Please read the [full text](https://github.com/Blosc/community/blob/master/code_of_conduct.md)
so that you can understand what actions will and will not be tolerated.
Zstandard library files
================================

The __lib__ directory is split into several sub-directories,
in order to make it easier to select or exclude features.


#### Building

`Makefile` script is provided, supporting [Makefile conventions](https://www.gnu.org/prep/standards/html_node/Makefile-Conventions.html#Makefile-Conventions),
including commands variables, staged install, directory variables and standard targets.
- `make` : generates both static and dynamic libraries
- `make install` : install libraries and headers in target system directories

`libzstd` default scope is pretty large, including compression, decompression, dictionary builder,
and support for decoding legacy formats >= v0.5.0.
The scope can be reduced on demand (see paragraph _modular build_).


#### Multithreading support

When building with `make`, by default the dynamic library is multithreaded and static library is single-threaded (for compatibility reasons).

Enabling multithreading requires 2 conditions :
- set build macro `ZSTD_MULTITHREAD` (`-DZSTD_MULTITHREAD` for `gcc`)
- for POSIX systems : compile with pthread (`-pthread` compilation flag for `gcc`)

For convenience, we provide a build target to generate multi and single threaded libraries:
- Force enable multithreading on both dynamic and static libraries by appending `-mt` to the target, e.g. `make lib-mt`.
- Force disable multithreading on both dynamic and static libraries by appending `-nomt` to the target, e.g. `make lib-nomt`.
- By default, as mentioned before, dynamic library is multithreaded, and static library is single-threaded, e.g. `make lib`.

When linking a POSIX program with a multithreaded version of `libzstd`,
note that it's necessary to invoke the `-pthread` flag during link stage.

Multithreading capabilities are exposed
via the [advanced API defined in `lib/zstd.h`](https://github.com/facebook/zstd/blob/v1.4.3/lib/zstd.h#L351).


#### API

Zstandard's stable API is exposed within [lib/zstd.h](zstd.h).


#### Advanced API

Optional advanced features are exposed via :

- `lib/zstd_errors.h` : translates `size_t` function results
                        into a `ZSTD_ErrorCode`, for accurate error handling.

- `ZSTD_STATIC_LINKING_ONLY` : if this macro is defined _before_ including `zstd.h`,
                          it unlocks access to the experimental API,
                          exposed in the second part of `zstd.h`.
                          All definitions in the experimental APIs are unstable,
                          they may still change in the future, or even be removed.
                          As a consequence, experimental definitions shall ___never be used with dynamic library___ !
                          Only static linking is allowed.


#### Modular build

It's possible to compile only a limited set of features within `libzstd`.
The file structure is designed to make this selection manually achievable for any build system :

- Directory `lib/common` is always required, for all variants.

- Compression source code lies in `lib/compress`

- Decompression source code lies in `lib/decompress`

- It's possible to include only `compress` or only `decompress`, they don't depend on each other.

- `lib/dictBuilder` : makes it possible to generate dictionaries from a set of samples.
        The API is exposed in `lib/dictBuilder/zdict.h`.
        This module depends on both `lib/common` and `lib/compress` .

- `lib/legacy` : makes it possible to decompress legacy zstd formats, starting from `v0.1.0`.
        This module depends on `lib/common` and `lib/decompress`.
        To enable this feature, define `ZSTD_LEGACY_SUPPORT` during compilation.
        Specifying a number limits versions supported to that version onward.
        For example, `ZSTD_LEGACY_SUPPORT=2` means : "support legacy formats >= v0.2.0".
        Conversely, `ZSTD_LEGACY_SUPPORT=0` means "do __not__ support legacy formats".
        By default, this build macro is set as `ZSTD_LEGACY_SUPPORT=5`.
        Decoding supported legacy format is a transparent capability triggered within decompression functions.
        It's also allowed to invoke legacy API directly, exposed in `lib/legacy/zstd_legacy.h`.
        Each version does also provide its own set of advanced API.
        For example, advanced API for version `v0.4` is exposed in `lib/legacy/zstd_v04.h` .

- While invoking `make libzstd`, it's possible to define build macros
        `ZSTD_LIB_COMPRESSION, ZSTD_LIB_DECOMPRESSION`, `ZSTD_LIB_DICTBUILDER`,
        and `ZSTD_LIB_DEPRECATED` as `0` to forgo compilation of the
        corresponding features. This will also disable compilation of all
        dependencies (eg. `ZSTD_LIB_COMPRESSION=0` will also disable
        dictBuilder).

- There are a number of options that can help minimize the binary size of
  `libzstd`.

  The first step is to select the components needed (using the above-described
  `ZSTD_LIB_COMPRESSION` etc.).

  The next step is to set `ZSTD_LIB_MINIFY` to `1` when invoking `make`. This
  disables various optional components and changes the compilation flags to
  prioritize space-saving.

  Detailed options: Zstandard's code and build environment is set up by default
  to optimize above all else for performance. In pursuit of this goal, Zstandard
  makes significant trade-offs in code size. For example, Zstandard often has
  more than one implementation of a particular component, with each
  implementation optimized for different scenarios. For example, the Huffman
  decoder has complementary implementations that decode the stream one symbol at
  a time or two symbols at a time. Zstd normally includes both (and dispatches
  between them at runtime), but by defining `HUF_FORCE_DECOMPRESS_X1` or
  `HUF_FORCE_DECOMPRESS_X2`, you can force the use of one or the other, avoiding
  compilation of the other. Similarly, `ZSTD_FORCE_DECOMPRESS_SEQUENCES_SHORT`
  and `ZSTD_FORCE_DECOMPRESS_SEQUENCES_LONG` force the compilation and use of
  only one or the other of two decompression implementations. The smallest
  binary is achieved by using `HUF_FORCE_DECOMPRESS_X1` and
  `ZSTD_FORCE_DECOMPRESS_SEQUENCES_SHORT` (implied by `ZSTD_LIB_MINIFY`).

  For squeezing the last ounce of size out, you can also define
  `ZSTD_NO_INLINE`, which disables inlining, and `ZSTD_STRIP_ERROR_STRINGS`,
  which removes the error messages that are otherwise returned by
  `ZSTD_getErrorName` (implied by `ZSTD_LIB_MINIFY`).

  Finally, when integrating into your application, make sure you're doing link-
  time optimation and unused symbol garbage collection (via some combination of,
  e.g., `-flto`, `-ffat-lto-objects`, `-fuse-linker-plugin`,
  `-ffunction-sections`, `-fdata-sections`, `-fmerge-all-constants`,
  `-Wl,--gc-sections`, `-Wl,-z,norelro`, and an archiver that understands
  the compiler's intermediate representation, e.g., `AR=gcc-ar`). Consult your
  compiler's documentation.

- While invoking `make libzstd`, the build macro `ZSTD_LEGACY_MULTITHREADED_API=1`
  will expose the deprecated `ZSTDMT` API exposed by `zstdmt_compress.h` in
  the shared library, which is now hidden by default.

- The build macro `DYNAMIC_BMI2` can be set to 1 or 0 in order to generate binaries
  which can detect at runtime the presence of BMI2 instructions, and use them only if present.
  These instructions contribute to better performance, notably on the decoder side.
  By default, this feature is automatically enabled on detecting
  the right instruction set (x64) and compiler (clang or gcc >= 5).
  It's obviously disabled for different cpus,
  or when BMI2 instruction set is _required_ by the compiler command line
  (in this case, only the BMI2 code path is generated).
  Setting this macro will either force to generate the BMI2 dispatcher (1)
  or prevent it (0). It overrides automatic detection.

- The build macro `ZSTD_NO_UNUSED_FUNCTIONS` can be defined to hide the definitions of functions
  that zstd does not use. Not all unused functions are hidden, but they can be if needed.
  Currently, this macro will hide function definitions in FSE and HUF that use an excessive
  amount of stack space.

- The build macro `ZSTD_NO_INTRINSICS` can be defined to disable all explicit intrinsics.
  Compiler builtins are still used.


#### Windows : using MinGW+MSYS to create DLL

DLL can be created using MinGW+MSYS with the `make libzstd` command.
This command creates `dll\libzstd.dll` and the import library `dll\libzstd.lib`.
The import library is only required with Visual C++.
The header file `zstd.h` and the dynamic library `dll\libzstd.dll` are required to
compile a project using gcc/MinGW.
The dynamic library has to be added to linking options.
It means that if a project that uses ZSTD consists of a single `test-dll.c`
file it should be linked with `dll\libzstd.dll`. For example:
```
    gcc $(CFLAGS) -Iinclude/ test-dll.c -o test-dll dll\libzstd.dll
```
The compiled executable will require ZSTD DLL which is available at `dll\libzstd.dll`.


#### Advanced Build options

The build system requires a hash function in order to
separate object files created with different compilation flags.
By default, it tries to use `md5sum` or equivalent.
The hash function can be manually switched by setting the `HASH` variable.
For example : `make HASH=xxhsum`
The hash function needs to generate at least 64-bit using hexadecimal format.
When no hash function is found,
the Makefile just generates all object files into the same default directory,
irrespective of compilation flags.
This functionality only matters if `libzstd` is compiled multiple times
with different build flags.

The build directory, where object files are stored
can also be manually controlled using variable `BUILD_DIR`,
for example `make BUILD_DIR=objectDir/v1`.
In which case, the hash function doesn't matter.


#### Deprecated API

Obsolete API on their way out are stored in directory `lib/deprecated`.
At this stage, it contains older streaming prototypes, in `lib/deprecated/zbuff.h`.
These prototypes will be removed in some future version.
Consider migrating code towards supported streaming API exposed in `zstd.h`.


#### Miscellaneous

The other files are not source code. There are :

 - `BUCK` : support for `buck` build system (https://buckbuild.com/)
 - `Makefile` : `make` script to build and install zstd library (static and dynamic)
 - `README.md` : this file
 - `dll/` : resources directory for Windows compilation
 - `libzstd.pc.in` : script for `pkg-config` (used in `make install`)
# ZSTD Windows binary package

## The package contents

- `zstd.exe` : Command Line Utility, supporting gzip-like arguments
- `dll\libzstd.dll` : The ZSTD dynamic library (DLL)
- `dll\libzstd.lib` : The import library of the ZSTD dynamic library (DLL) for Visual C++
- `example\` : The example of usage of the ZSTD library
- `include\` : Header files required by the ZSTD library
- `static\libzstd_static.lib` : The static ZSTD library (LIB)

## Usage of Command Line Interface

Command Line Interface (CLI) supports gzip-like arguments.
By default CLI takes an input file and compresses it to an output file:

    Usage: zstd [arg] [input] [output]

The full list of commands for CLI can be obtained with `-h` or `-H`. The ratio can
be improved with commands from `-3` to `-16` but higher levels also have slower
compression. CLI includes in-memory compression benchmark module with compression
levels starting from `-b` and ending with `-e` with iteration time of `-i` seconds.
CLI supports aggregation of parameters i.e. `-b1`, `-e18`, and `-i1` can be joined
into `-b1e18i1`.

## The example of usage of static and dynamic ZSTD libraries with gcc/MinGW

Use `cd example` and `make` to build `fullbench-dll` and `fullbench-lib`.
`fullbench-dll` uses a dynamic ZSTD library from the `dll` directory.
`fullbench-lib` uses a static ZSTD library from the `lib` directory.

## Using ZSTD DLL with gcc/MinGW

The header files from `include\` and the dynamic library `dll\libzstd.dll`
are required to compile a project using gcc/MinGW.
The dynamic library has to be added to linking options.
It means that if a project that uses ZSTD consists of a single `test-dll.c`
file it should be linked with `dll\libzstd.dll`. For example:

    gcc $(CFLAGS) -Iinclude\ test-dll.c -o test-dll dll\libzstd.dll

The compiled executable will require ZSTD DLL which is available at `dll\libzstd.dll`.

## The example of usage of static and dynamic ZSTD libraries with Visual C++

Open `example\fullbench-dll.sln` to compile `fullbench-dll` that uses a
dynamic ZSTD library from the `dll` directory. The solution works with Visual C++
2010 or newer. When one will open the solution with Visual C++ newer than 2010
then the solution will upgraded to the current version.

## Using ZSTD DLL with Visual C++

The header files from `include\` and the import library `dll\libzstd.lib`
are required to compile a project using Visual C++.

1. The path to header files should be added to `Additional Include Directories` that can
   be found in project properties `C/C++` then `General`.
2. The import library has to be added to `Additional Dependencies` that can
   be found in project properties `Linker` then `Input`.
   If one will provide only the name `libzstd.lib` without a full path to the library
   the directory has to be added to `Linker\General\Additional Library Directories`.

The compiled executable will require ZSTD DLL which is available at `dll\libzstd.dll`.
