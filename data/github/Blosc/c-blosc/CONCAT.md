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
C-Blosc libraries come with Python-Blosc wheels
===============================================

Starting on version 1.21.0, C-Blosc binary libraries can easily be installed from Python-Blosc (>= 1.10) wheels:

.. code-block:: console

        $ pip install blosc                                                                             (base)
        Collecting blosc
          Downloading blosc-1.10.0-cp37-cp37m-macosx_10_9_x86_64.whl (2.2 MB)
             |████████████████████████████████| 2.2 MB 4.7 MB/s
        Installing collected packages: blosc
          Attempting uninstall: blosc
            Found existing installation: blosc 1.10.0
            Uninstalling blosc-1.10.0:
              Successfully uninstalled blosc-1.10.0
        Successfully installed blosc-1.10.0

As a result, one can easily update to the latest version of C-Blosc binaries without the need to manually compile the thing.  Following are instructions on how to use the libraries in wheels for different platforms.


Compiling C files with Blosc wheels on Windows
----------------------------------------------

- The wheels for Windows have been produced with the Microsoft MSVC compiler, so we recommend that you use it too.  You can get it for free at: https://visualstudio.microsoft.com/es/downloads/.

- In order to check that the MSVC command line is set up correctly, enter ``cl`` in the command prompt window and verify that the output looks something like this:

.. code-block:: console

    > cl
    Microsoft (R) C/C++ Optimizing Compiler Version 19.00.24245 for x64
    Copyright (C) Microsoft Corporation.  All rights reserved.

    usage: cl [ option... ] filename... [ /link linkoption... ]

- Now, install the wheels:

.. code-block:: console

    > pip install blosc
    Collecting blosc
      Using cached blosc-1.10.0-cp37-cp37m-win_amd64.whl (1.5 MB)
    Installing collected packages: blosc
    Successfully installed blosc-1.10.0

- Make the compiler available. Its typical installation location uses to be `C:\\Program files (x86)\\Microsoft Visual Studio`, so change your current directory there. Then, to set up the build architecture environment you can open a command prompt window in the `VC\\Auxiliary\\Build` subdirectory and execute `vcvarsall.bat x64` if your architecture is 64 bits or `vcvarsall.bat x86` if it is 32 bits.

- You will need to know the path where the Blosc wheel has installed its files.  For this we will use the `dir /s` command (but you can use your preferred location method):

.. code-block:: console

    > dir /s c:\blosc.lib
     Volume in drive C is OS
     Volume Serial Number is 7A21-A5D5

     Directory of c:\Users\user\miniconda3\Lib

    14/12/2020  09:56             7.022 blosc.lib
                   1 File(s)          7.022 bytes

         Total list files:
                   1 File(s)          7.022 bytes
                   0 dirs  38.981.902.336 free bytes

- The output shows the path of blosc.lib in your system, but we are rather interested in the parent one:

.. code-block:: console

    > set WHEEL_DIR=c:\Users\user\miniconda3

- Now, it is important to copy the library `blosc.dll` to C:\\Windows\\System32 directory, so it can be found by the executable when it is necessary.

- Finally, to compile C files using Blosc libraries, enter this command:

.. code-block:: console

    > cl <file_name>.c <path_of_blosc.lib> /Ox /Fe<file_name>.exe /I<path_of_blosc.h> /MT /link/NODEFAULTLIB:MSVCRT

- For instance, in the case of blosc "examples/simple.c":

.. code-block:: console

    > cl simple.c %WHEEL_DIR%\lib\blosc.lib /Ox /Fesimple.exe /I%WHEEL_DIR%\include /MT /link/NODEFAULTLIB:MSVCRT

    Microsoft (R) C/C++ Optimizing Compiler Version 19.10.25017 for x86
    Copyright (C) Microsoft Corporation.  All rights reserved.

    simple.c
    Microsoft (R) Incremental Linker Version 14.10.25017.0
    Copyright (C) Microsoft Corporation.  All rights reserved.

    /out:simple.exe
    simple.obj
    /NODEFAULTLIB:MSVCRT
    .\miniconda3\lib\blosc.lib

- And you can run your program:

.. code-block:: console

    > simple
    Blosc version info: 1.20.1 ($Date:: 2020-09-08 #$)
    Compression: 4000000 -> 37816 (105.8x)
    Decompression successful!
    Successful roundtrip!

- Rejoice!


Compiling C files with Blosc wheels on Linux
--------------------------------------------

- Install the wheels:

.. code-block:: console

    $ pip install blosc
    Collecting blosc
      Using cached blosc-1.10.0-cp37-cp37m-manylinux2010_x86_64.whl (2.2 MB)
    Installing collected packages: blosc
    Successfully installed blosc-1.10.0

- Find the path where blosc wheel has installed its files:

.. code-block:: console

    $ find / -name libblosc.so 2>/dev/null
    /home/soscar/miniconda3/lib/libblosc.so

- The output shows the path of libblosc.so, but we are rather interested in the parent one:

.. code-block:: console

    $ WHEEL_DIR=/home/soscar/miniconda3

- To compile C files using blosc you only need to enter the commands:

.. code-block:: console

    $ export LD_LIBRARY_PATH=<path_of_libblosc.so>
    $ gcc <file_name>.c -I<path_of_blosc.h> -o <file_name> -L<path_of_libblosc.so> -lblosc

- For instance, let's compile blosc's "examples/many_compressors.c":

.. code-block:: console

    $ export LD_LIBRARY_PATH=$WHEEL_DIR/lib   # note that you need the LD_LIBRARY_PATH env variable
    $ gcc many_compressors.c -I$WHEEL_DIR/include -o many_compressors -L$WHEEL_DIR/lib -lblosc

- Run your program:

.. code-block:: console

    $ ./many_compressors
    Blosc version info: 1.20.1 ($Date:: 2020-09-08 #$)
    Using 4 threads (previously using 1)
    Using blosclz compressor
    Compression: 4000000 -> 37816 (105.8x)
    Successful roundtrip!
    Using lz4 compressor
    Compression: 4000000 -> 37938 (105.4x)
    Successful roundtrip!
    Using lz4hc compressor
    Compression: 4000000 -> 27165 (147.2x)
    Successful roundtrip!

- Rejoice!


Compiling C files with Blosc wheels on MacOS
--------------------------------------------

- Install the wheels:

.. code-block:: console

        $ pip install blosc                                                                             (base)
        Collecting blosc
          Downloading blosc-1.10.0-cp37-cp37m-macosx_10_9_x86_64.whl (2.2 MB)
             |████████████████████████████████| 2.2 MB 4.7 MB/s
        Installing collected packages: blosc
          Attempting uninstall: blosc
            Found existing installation: blosc 1.10.0
            Uninstalling blosc-1.10.0:
              Successfully uninstalled blosc-1.10.0
        Successfully installed blosc-1.10.0

- Find the path where blosc wheel has installed its files:

.. code-block:: console

    $ find / -name libblosc.dylib 2>/dev/null
    /home/soscar/miniconda3/lib/libblosc.dylib

- The output shows the path of libblosc.dylib, but we are rather interested in the parent one:

.. code-block:: console

    $ WHEEL_DIR=/home/soscar/miniconda3

- To compile C files using blosc you only need to enter the commands:

.. code-block:: console

    $ export LD_LIBRARY_PATH=<path_of_libblosc.dylib>
    $ clang <file_name>.c -I<path_of_blosc.h> -o <file_name> -L<path_of_libblosc.dylib> -lblosc

- For instance, let's compile blosc's "examples/many_compressors.c":

.. code-block:: console

    $ export LD_LIBRARY_PATH=$WHEEL_DIR/lib   # note that you need the LD_LIBRARY_PATH env variable
    $ clang many_compressors.c -I$WHEEL_DIR/include -o many_compressors -L$WHEEL_DIR/lib -lblosc

- Run your program:

.. code-block:: console

    $ ./many_compressors
    Blosc version info: 1.20.1 ($Date:: 2020-09-08 #$)
    Using 4 threads (previously using 1)
    Using blosclz compressor
    Compression: 4000000 -> 37816 (105.8x)
    Successful roundtrip!
    Using lz4 compressor
    Compression: 4000000 -> 37938 (105.4x)
    Successful roundtrip!
    Using lz4hc compressor
    Compression: 4000000 -> 27165 (147.2x)
    Successful roundtrip!

- Rejoice!
===========================
 Release notes for C-Blosc
===========================

Changes from 1.21.1 to 1.21.2
=============================

#XXX version-specific blurb XXX#


Changes from 1.21.0 to 1.21.1
=============================

* Fix pthread flag when linking on ppc64le.  See #318.  Thanks to Axel Huebl.

* Updates in codecs (some bring important performance improvements):
  * BloscLZ updated to 2.5.1
  * Zlib updated to 1.2.11
  * Zstd updated to 1.5.0


Changes from 1.20.1 to 1.21.0
=============================

* Updated zstd codec to 1.4.8.

* Updated lz4 codec to 1.9.3.

* New instructions on how to use the libraries in python-blosc wheels
  so as to compile C-Blosc applications.  See:
  https://github.com/Blosc/c-blosc/blob/master/COMPILING_WITH_WHEELS.rst


Changes from 1.20.0 to 1.20.1
=============================

* Added `<unistd.h>` in vendored zlib 1.2.8 for compatibility with Python 3.8
  in recent Mac OSX.  For details, see:
  https://github.com/Blosc/python-blosc/issues/229


Changes from 1.19.1 to 1.20.0
=============================

* More safety checks have been implemented so that potential flaws
  discovered by new fuzzers in OSS-Fuzzer are fixed now.  Thanks to
  Nathan Moinvaziri (@nmoinvaz).

* BloscLZ updated to 2.3.0. Expect better compression ratios for faster
  codecs.  For details, see our new blog post:
  https://blosc.org/posts/beast-release/

* Fixed the `_xgetbv()` collision. Thanks to Michał Górny (@mgorny).

* The chunk format has been fully described so that 3rd party software
  may come with a different implementation, but still compatible with
  C-Blosc chunks.


Changes from 1.19.0 to 1.19.1
=============================

- pthread_create() errors are now handled and propagated back to the user.
  See https://github.com/Blosc/c-blosc/pull/299.


Changes from 1.18.1 to 1.19.0
=============================

- The length of automatic blocksizes for fast codecs (lz4, blosclz) has
  been incremented quite a bit (up to 256 KB) for better compression ratios.
  The performance in modern CPUs (with at least 256 KB in L2 cache) should
  be better too (for older CPUs the performance should stay roughly the same).

- Continuous integration has been migrated to GitHub actions and much
  more scenarios are tested (specially linking with external codecs).
  Also, a new OSS-Fuzz workflow has been added for increased detection
  of possible vulnerabilities.  Thanks to Nathan Moinvaziri.

- For small buffers that cannot be compressed (typically < 128 bytes),
  `blosc_compress()` returns now a 0 (cannot compress) instead of a negative
  number (internal error).  See https://github.com/Blosc/c-blosc/pull/294.
  Thanks to @kalvdans for providing the initial patch.

- blosclz codec updated to 2.1.0.  Expect better compression ratios and
  performance in a wider variety of scenarios.

- `blosc_decompress_unsafe()`, `blosc_decompress_ctx_unsafe()` and
  `blosc_getitem_unsafe()` have been removed because they are dangerous
  and after latest improvements, they should not be used in production.

- zstd codec updated to 1.4.5.

- Conan packaging has been deprecated (from now on, we should try
  to focus on supporting wheels only).


Changes from 1.17.1 to 1.18.1
=============================

- Fixed the copy of the leftovers of a chunk when its size is not a
  multiple of the typesize.  Although this is a very unusual situation,
  it can certainly happen (e.g.
  https://github.com/Blosc/python-blosc/issues/220).


Changes from 1.17.0 to 1.17.1
=============================

- Zstd codec updated to 1.4.4.

- LZ4 codec updated to 1.9.2.


Changes from 1.16.3 to 1.17.0
=============================

- LZ4 codec updated to 1.9.1.

- Zstd codec updated to 1.4.1.

- BloscLZ codec updated to 2.0.0.  Although this should be fully backward
  compatible, it contains important changes that affects mainly speed, but
  also compression ratios.  Feedback on how it behaves on your own data is
  appreciated.


Changes from 1.16.2 to 1.16.3
=============================

- Fix for building for clang with -march=haswell. See PR #262.

- Fix all the known warnings for GCC/Clang.  Still some work to do for MSVC
  in this front.

- Due to some problems with several CI systems, the check for library symbols
  are deactivated now by default.  If you want to enforce this check, use:
  `cmake .. -DDEACTIVATE_SYMBOLS_CHECK=ON` to re-activate it.


Changes from 1.16.1 to 1.16.2
=============================

- Correct the check for the compressed size when the buffer is memcpyed.  This
  was a regression introduced in 1.16.0.  Fixes #261.


Changes from 1.16.0 to 1.16.1
=============================

- Fixed a regression in 1.16.0 that prevented to compress empty buffers
  (see #260).

- Zstd updated to 1.3.8 (from 1.3.7).


Changes from 1.15.1 to 1.16.0
=============================

- Now the functions that execute Blosc decompressions are safe by default
  for untrusted/possibly corrupted inputs.  The additional checks seem to
  not affect performance significantly (see some benchmarks in #258), so
  this is why they are the default now.

  The previous functions (with less safety) checks are still available with a
  '_unsafe' suffix.  The complete list is:

    - blosc_decompress_unsafe()
    - blosc_decompress_ctx_unsafe()
    - blosc_getitem_unsafe()

  Also, a new API function named blosc_cbuffer_validate(), for validating Blosc
  compressed data, has been added.

  For details, see PR #258.  Thanks to Jeremy Maitin-Shepard.

- Fixed a bug in `blosc_compress()` that could lead to thread deadlock under
  some situations.  See #251.  Thanks to @wenjuno for the report and the fix.

- Fix data race in shuffle.c host_implementation initialization.  Fixes #253.
  Thanks to Jeremy Maitin-Shepard.


Changes from 1.15.0 to 1.15.1
=============================

- Add workaround for Visual Studio 2008's lack of a `stdint.h` file to
  `blosclz.c`.


Changes from 1.14.4 to 1.15.0
=============================

- The `blosc_compress()` and `blosc_decompress()` interfaces are now
  fork-safe, preventing child-process deadlocks in fork-based
  multiprocessing applications. These interfaces with BLOSC_NOLOCK were, and
  continue to be, fork-safe. `_ctx` interface context reuse continues to be
  unsafe in the child process post-fork. See #241.  Thanks to Alex Ford.

- Replaced //-comments with /**/-comments and other improvements for
  compatibility with quite old gcc compilers.  See PR #243.  Thanks to
  Andreas Martin.

- Empty buffers can be compressed again (this was inadvertently prevented while
  fixing #234).  See #247.  Thanks to Valentin Haenel.

- LZ4 internal codec upgraded to 1.8.3 (from 1.8.1.2).

- Zstd internal codec upgraded to 1.3.7 (from 1.3.4).


Changes from 1.14.3 to 1.14.4
=============================

- Added a new `DEACTIVATE_SSE2` option for cmake that is useful for disabling
  SSE2 when doing cross-compilation (see #236).

- New check for detecting output buffers smaller than BLOSC_MAX_OVERHEAD.
  Fixes #234.

- The `complib` and `version` parameters for `blosc_get_complib_info()` can be
  safely set to NULL now.  This allows to call this function even if the user is
  not interested in these parameters (so no need to reserve memory for them).
  Fixes #228.

- In some situations that a supposedly blosc chunk is passed to
  `blosc_decompress()`, one might end with an `Arithmetic exception`.  This
  is probably due to the chunk not being an actual blosc chunk, and divisions
  by zero might occur.  A protection has been added for this. See #237.


Changes from 1.14.2 to 1.14.3
=============================

- Use win32/pthread.c on all Windows builds, even those with GNU compilers.
  Rational: although MinGW provides a more full-featured pthreads replacement,
  it doesn't seem to accomplish anything here since the functionality in
  win32/pthread.c is sufficient for Blosc. Furthermore, using the MinGW
  pthreads adds an additional library dependency to libblosc that is
  annoying for binary distribution. For example, it got in the way of
  distributing cross-compiled Windows binaries for use with Julia, since they
  want the resulting libblosc.dll to be usable on any Windows machine even
  where MinGW is not installed.  See PR #224.  Thanks to Steven G. Johnson.

- Zstd internal sources have been updated to 1.3.4.


Changes from 1.14.1 to 1.14.2
=============================

- Reverted the $Configuration var in CMake configuration for Windows so
  as to restore the compatibility with MS VisualStudio compilers.


Changes from 1.14.0 to 1.14.1
=============================

- Fixed a bug that caused C-Blosc to crash on platforms requiring strict
  alignment (as in some kinds of ARM CPUs).  Fixes #223.  Thanks to Elvis
  Stansvik and Michael Hudson-Doyle for their help.

- Fixed a piece of code that was not C89 compliant.  C89 compliance is
  needed mainly by MS VS2008 which is still used for creating Python 2
  extensions.

- Remove the (spurious) $Configuration var in cmake config for Windows.
  Thanks to Francis Brissette for pointing this out.


Changes from 1.13.7 to 1.14.0
=============================

- New split mode that favors forward compatibility.  That means that,
  from now on, all the buffers created starting with blosc 1.14.0 will
  be forward compatible with any previous versions of the library --at
  least until 1.3.0, when support for multi-codec was introduced.

  So as to select the split mode, a new API function has been introduced:
  https://github.com/Blosc/c-blosc/blob/master/blosc/blosc.h#L500
  Also, the BLOSC_SPLITMODE environment variable is honored when using
  the `blosc_compress()` function.  See
  https://github.com/Blosc/c-blosc/blob/master/blosc/blosc.h#L209

  There is a dedicated blog entry about this at:
  http://blosc.org/posts/new-forward-compat-policy/
  More info in PR #216.

  Caveat Emptor: Note that Blosc versions from 1.11.0 to 1.14.0 *might*
  generate buffers that cannot be read with versions < 1.11.0, so if
  forward compatibility is important to you, an upgrade to 1.14.0 is
  recommended.

- All warnings during cmake build stage are enabled by default now.
  PR #218.  Thanks to kalvdans.

- Better checks on versions of formats inside Blosc.  PR #219.  Thanks
  to kalvdans.

- The BLOSC_PRINT_SHUFFLE_ACCEL environment variable is honored now.
  This is useful for determining *at runtime* whether the different SIMD
  capabilities (only for x86 kind processors) are available to Blosc to get
  better performance during shuffle/bitshuffle operation.  As an example,
  here it is the normal output for the simple.c example::

    $ ./simple
    Blosc version info: 1.14.0.dev ($Date:: 2018-02-15 #$)
    Compression: 4000000 -> 41384 (96.7x)
    Decompression successful!
    Successful roundtrip!

  and here with the BLOSC_PRINT_SHUFFLE_ACCEL environment variable set::

    $ BLOSC_PRINT_SHUFFLE_ACCEL= ./simple
    Blosc version info: 1.14.0.dev ($Date:: 2018-02-15 #$)
    Shuffle CPU Information:
    SSE2 available: True
    SSE3 available: True
    SSSE3 available: True
    SSE4.1 available: True
    SSE4.2 available: True
    AVX2 available: True
    AVX512BW available: False
    XSAVE available: True
    XSAVE enabled: True
    XMM state enabled: True
    YMM state enabled: True
    ZMM state enabled: False
    Compression: 4000000 -> 41384 (96.7x)
    Decompression successful!
    Successful roundtrip!

  Blosc only currently leverages the SSE2 and AVX2 instruction sets, but
  it can recognize all of the above.  This is useful mainly for debugging.


Changes from 1.13.6 to 1.13.7
=============================

- More tests for binaries in https://bintray.com/blosc/Conan.


Changes from 1.13.5 to 1.13.6
=============================

- More tests for binaries in https://bintray.com/blosc/Conan.


Changes from 1.13.4 to 1.13.5
=============================

- New conan binaries publicly accessible in https://bintray.com/blosc/Conan.
  Still experimental, but feedback is appreciated.


Changes from 1.13.3 to 1.13.4
=============================

- Fixed a buffer overrun that happens when compressing small buffers and
  len(destination_buffer) < (len(source_buffer) + BLOSC_MAX_OVERHEAD).
  Reported by Ivan Smirnov.


Changes from 1.13.2 to 1.13.3
=============================

- Tests work now when external compressors are located in non-system locations.
  Fixes #210.  Thanks to Leif Walsh.


Changes from 1.13.1 to 1.13.2
=============================

- C-Blosc can be compiled on CentOS 6 now.

- LZ4 internal codec upgraded to 1.8.1.


Changes from 1.13.0 to 1.13.1
=============================

- Fixed a bug uncovered by the python-blosc test suite: when a buffer is
  to be copied, then we should reserve space for the header, not block pointers.


Changes from 1.12.1 to 1.13.0
=============================

- Serious optimization of memory copy functions (see new `blosc/fastcopy.c`).
  This benefits the speed of all the codecs, but specially the BloscLZ one.

- As a result of the above, the BloscLZ codec received a new adjustment of
  knobs so that you should expect better compression ratios with it too.

- LZ4 internal sources have been updated to 1.8.0.

- Zstd internal sources have been updated to 1.3.3.


Changes from 1.12.0 to 1.12.1
=============================

- Backported BloscLZ parameters that were fine-tuned for C-Blosc2.
  You should expect better compression ratios and faster operation,
  specially on modern CPUs.  See:
  http://blosc.org/posts/blosclz-tuning/


Changes from 1.11.3 to 1.12.0
=============================

- Snappy, Zlib and Zstd codecs are compiled internally now, even if they are
  installed in the machine.  This has been done in order to avoid
  problems in machines having the shared libraries for the codecs
  accessible but not the includes (typical in Windows boxes).  Also,
  the Zstd codec runs much faster when compiled internally.  The
  previous behaviour can be restored by activating the cmake options
  PREFER_EXTERNAL_SNAPPY, PREFER_EXTERNAL_ZLIB and PREFER_EXTERNAL_ZSTD.

- Zstd internal sources have been updated to 1.3.0.


Changes from 1.11.3 to 1.11.4
=============================

- Internal Zstd codec updated to 1.1.4.


Changes from 1.11.2 to 1.11.3
=============================

- Fixed #181: bitshuffle filter for big endian machines.

- Internal Zstd codec updated to 1.1.3.

- New blocksize for complevel 8 in automatic mode.  This should help specially
  the Zstd codec to achieve better compression ratios.


Changes from 1.11.1 to 1.11.2
=============================

- Enabled use as a CMake subproject, exporting shared & static library targets
  for super-projects to use. See PRs #178, #179 and #180.  Thanks to Kevin
  Murray.

- Internal LZ4 codec updated to 1.7.5.

- Internal Zstd codec updated to 1.1.2.


Changes from 1.11.0 to 1.11.1
=============================

- Fixed a bug introduced in 1.11.0 and discovered by pandas test suite. This
  basically prevented to decompress buffers compressed with previous versions of
  C-Blosc. See: https://github.com/Blosc/python-blosc/issues/115


Changes from 1.10.2 to 1.11.0
=============================

- Internal Zstd codec upgraded to 1.0.0.

- New block size computation inherited from C-Blosc2. Benchmarks are saying that
  this benefits mainly to LZ4, LZ4HC, Zlib and Zstd codecs, both in speed and in
  compression ratios (although YMMV for your case).

- Added the @rpath flag in Mac OSX for shared libraries.  Fixes #175.

- Added a fix for VS2008 discovered in: https://github.com/PyTables/PyTables/pull/569/files#diff-953cf824ebfea7208d2a2e312d9ccda2L126

- License changed from MIT to 3-clause BSD style.


Changes from 1.10.1 to 1.10.2
=============================

- Force the use of --std=gnu99 when using gcc.  Fixes #174.


Changes from 1.10.0 to 1.10.1
=============================

- Removed an inconsistent check for C11 (__STDC_VERSION__ >= 201112L and
  _ISOC11_SOURCE) as this seem to pose problems on compilers doing different
  things in this check (e.g. clang). See
  https://github.com/Blosc/bloscpack/issues/50.


Changes from 1.9.3 to 1.10.0
============================

- Initial support for Zstandard (0.7.4). Zstandard (or Zstd for short) is a new
  compression library that allows better compression than Zlib, but that works
  typically faster (and some times much faster), making of it a good match for
  Blosc.

  Although the Zstd format is considered stable
  (http://fastcompression.blogspot.com.es/2016_07_03_archive.html), its API is
  maturing very fast, and despite passing the extreme test suite for C-Blosc,
  this codec should be considered in beta for C-Blosc usage purposes. Please
  test it and report back any possible issues you may get.


Changes from 1.9.2 to 1.9.3
===========================

- Reverted a mistake introduced in 1.7.1.  At that time, bit-shuffling
  was enabled for typesize == 1 (i.e. strings), but the change also
  included byte-shuffling accidentally.  This only affected performance,
  but in a quite bad way (a copy was needed).  This has been fixed and
  byte-shuffling is not active when typesize == 1 anymore.


Changes from 1.9.1 to 1.9.2
===========================

- Check whether Blosc is actually initialized before blosc_init(),
  blosc_destroy() and blosc_free_resources().  This makes the library
  more resistant to different initialization cycles
  (e.g. https://github.com/stevengj/Blosc.jl/issues/19).


Changes from 1.9.0 to 1.9.1
===========================

- The internal copies when clevel=0 are made now via memcpy().  At the
  beginning of C-Blosc development, benchmarks where saying that the
  internal, multi-threaded copies inside C-Blosc were faster than
  memcpy(), but 6 years later, memcpy() made greats strides in terms
  of efficiency.  With this, you should expect an slight speed
  advantage (10% ~ 20%) when C-Blosc is used as a replacement of
  memcpy() (which should not be the most common scenario out there).

- Added a new DEACTIVATE_AVX2 cmake option to explicitly disable AVX2
  at build-time.  Thanks to James Bird.

- The ``make -jN`` for parallel compilation should work now.  Thanks
  to James Bird.


Changes from 1.8.1 to 1.9.0
===========================

* New blosc_get_nthreads() function to get the number of threads that
  will be used internally during compression/decompression (set by
  already existing blosc_set_nthreads()).

* New blosc_get_compressor() function to get the compressor that will
  be used internally during compression (set by already existing
  blosc_set_compressor()).

* New blosc_get_blocksize() function to get the internal blocksize to
  be used during compression (set by already existing
  blosc_set_blocksize()).

* Now, when the BLOSC_NOLOCK environment variable is set (to any
  value), the calls to blosc_compress() and blosc_decompress() will
  call blosc_compress_ctx() and blosc_decompress_ctx() under the hood
  so as to avoid the internal locks.  See blosc.h for details.  This
  allows multi-threaded apps calling the non _ctx() functions to avoid
  the internal locks in C-Blosc.  For the not multi-threaded app
  though, it is in general slower to call the _ctx() functions so the
  use of BLOSC_NOLOCK is discouraged.

* In the same vein, from now on, when the BLOSC_NTHREADS environment
  variable is set to an integer, every call to blosc_compress() and
  blosc_decompress() will call blosc_set_nthreads(BLOSC_NTHREADS)
  before the actual compression/decompression process.  See blosc.h
  for details.

* Finally, if BLOSC_CLEVEL, BLOSC_SHUFFLE, BLOSC_TYPESIZE and/or
  BLOSC_COMPRESSOR variables are set in the environment, these will be
  also honored before calling blosc_compress().

* Calling blosc_init() before any other Blosc call, although
  recommended, is not necessary anymore.  The idea is that you can use
  just the basic blosc_compress() and blosc_decompress() and control
  other parameters (nthreads, compressor, blocksize) by using
  environment variables (see above).


Changes from 1.8.0 to 1.8.1
===========================

* Disable the use of __builtin_cpu_supports() for GCC 5.3.1
  compatibility.  Details in:
  https://lists.fedoraproject.org/archives/list/devel@lists.fedoraproject.org/thread/ZM2L65WIZEEQHHLFERZYD5FAG7QY2OGB/


Changes from 1.7.1 to 1.8.0
===========================

* The code is (again) compatible with VS2008 and VS2010.  This is
  important for compatibility with Python 2.6/2.7/3.3/3.4.

* Introduced a new global lock during blosc_decompress() operation.
  As the blosc_compress() was already guarded by a global lock, this
  means that the compression/decompression is again thread safe.
  However, when using C-Blosc from multi-threaded environments, it is
  important to keep using the *_ctx() functions for performance
  reasons.  NOTE: _ctx() functions will be replaced by more powerful
  ones in C-Blosc 2.0.


Changes from 1.7.0 to 1.7.1
===========================

* Fixed a bug preventing bitshuffle to work correctly on getitem().
  Now, everything with bitshuffle seems to work correctly.

* Fixed the thread initialization for blosc_decompress_ctx().  Issue
  #158.  Thanks to Chris Webers.

* Fixed a bug in the blocksize computation introduced in 1.7.0.  This
  could have been creating segfaults.

* Allow bitshuffle to run on 1-byte typesizes.

* New parametrization of the blocksize to be independent of the
  typesize.  This allows a smoother speed throughout all typesizes.

* lz4 and lz4hc codecs upgraded to 1.7.2 (from 1.7.0).

* When calling set_nthreads() but not actually changing the number of
  threads in the internal pool does not teardown and setup it anymore.
  PR #153.  Thanks to Santi Villalba.


Changes from 1.6.1 to 1.7.0
===========================

* Added a new 'bitshuffle' filter so that the shuffle takes place at a
  bit level and not just at a byte one, which is what it does the
  previous 'shuffle' filter.

  For activating this new bit-level filter you only have to pass the
  symbol BLOSC_BITSHUFFLE to `blosc_compress()`.  For the previous
  byte-level one, pass BLOSC_SHUFFLE.  For disabling the shuffle, pass
  BLOSC_NOSHUFFLE.

  This is a port of the existing filter in
  https://github.com/kiyo-masui/bitshuffle.  Thanks to Kiyo Masui for
  changing the license and allowing its inclusion here.

* New acceleration mode for LZ4 and BloscLZ codecs that enters in
  operation with complevel < 9.  This allows for an important boost in
  speed with minimal compression ratio loss.  Francesc Alted.

* LZ4 codec updated to 1.7.0 (r130).

* PREFER_EXTERNAL_COMPLIBS cmake option has been removed and replaced
  by the more fine grained PREFER_EXTERNAL_LZ4, PREFER_EXTERNAL_SNAPPY
  and PREFER_EXTERNAL_ZLIB.  In order to allow the use of the new API
  introduced in LZ4 1.7.0, PREFER_EXTERNAL_LZ4 has been set to OFF by
  default, whereas PREFER_EXTERNAL_SNAPPY and PREFER_EXTERNAL_ZLIB
  continues to be ON.

* Implemented SSE2 shuffle support for buffers containing a number of
  elements which is not a multiple of (typesize * vectorsize).  Jack
  Pappas.

* Added SSE2 shuffle/unshuffle routines for types larger than 16
  bytes.  Jack Pappas.

* 'test_basic' suite has been split in components for a much better
  granularity on what's a possibly failing test.  Also, lots of new
  tests have been added.  Jack Pappas.

* Fixed compilation on non-Intel archs (tested on ARM).  Zbyszek
  Szmek.

* Modifyied cmake files in order to inform that AVX2 on Visual Studio
  is supported only in 2013 update 2 and higher.

* Added a replacement for stdbool.h for Visual Studio < 2013.

* blosclz codec adds Win64/Intel as a platform supporting unaligned
  addressing.  That leads to a speed-up of 2.2x in decompression.

* New blosc_get_version_string() function for retrieving the version
  of the c-blosc library.  Useful when linking with dynamic libraries
  and one want to know its version.

* New example (win-dynamic-linking.c) that shows how to link a Blosc
  DLL dynamically in run-time (Windows only).

* The `context.threads_started` is initialized now when decompressing.
  This could cause crashes in case you decompressed before compressing
  (e.g. directly deserializing blosc buffers).  @atchouprakov.

* The HDF5 filter has been removed from c-blosc and moved into its own
  repo at: https://github.com/Blosc/hdf5

* The MS Visual Studio 2008 has been tested with c-blosc for ensuring
  compatibility with extensions for Python 2.6 and up.


Changes from 1.6.0 to 1.6.1
===========================

* Support for *runtime* detection of AVX2 and SSE2 SIMD instructions.
  These changes make it possible to compile one single binary that
  runs on a system that supports SSE2 or AVX2 (or neither), so the
  redistribution problem is fixed (see #101).  Thanks to Julian Taylor
  and Jack Pappas.

* Added support for MinGW and TDM-GCC compilers for Windows.  Thanks
  to yasushima-gd.

* Fixed a bug in blosclz that could potentially overwrite an area
  beyond the output buffer.  See #113.

* New computation for blocksize so that larger typesizes (> 8 bytes)
  would benefit of much better compression ratios.  Speed is not
  penalized too much.

* New parametrization of the hash table for blosclz codec.  This
  allows better compression in many scenarios, while slightly
  increasing the speed.


Changes from 1.5.4 to 1.6.0
===========================

* Support for AVX2 is here!  The benchmarks with a 4-core Intel
  Haswell machine tell that both compression and decompression are
  accelerated around a 10%, reaching peaks of 9.6 GB/s during
  compression and 26 GB/s during decompression (memcpy() speed for
  this machine is 7.5 GB/s for writes and 11.7 GB/s for reads).  Many
  thanks to @littlezhou for this nice work.

* Support for HPET (high precision timers) for the `bench` program.
  This is particularly important for microbenchmarks like bench is
  doing; since they take so little time to run, the granularity of a
  less-accurate timer may account for a significant portion of the
  runtime of the benchmark itself, skewing the results.  Thanks to
  Jack Pappas.


Changes from 1.5.3 to 1.5.4
===========================

* Updated to LZ4 1.6.0 (r128).

* Fix resource leak in t_blosc.  Jack Pappas.

* Better checks during testing.  Jack Pappas.

* Dynamically loadable HDF5 filter plugin. Kiyo Masui.


Changes from 1.5.2 to 1.5.3
===========================

* Use llabs function (where available) instead of abs to avoid
  truncating the result.  Jack Pappas.

* Use C11 aligned_alloc when it's available.  Jack Pappas.

* Use the built-in stdint.h with MSVC when available.  Jack Pappas.

* Only define the __SSE2__ symbol when compiling with MS Visual C++
  and targeting x64 or x86 with the correct /arch flag set. This
  avoids re-defining the symbol which makes other compilers issue
  warnings.  Jack Pappas.

* Reinitializing Blosc during a call to set_nthreads() so as to fix
  problems with contexts.  Francesc Alted.



Changes from 1.5.1 to 1.5.2
===========================

* Using blosc_compress_ctx() / blosc_decompress_ctx() inside the HDF5
  compressor for allowing operation in multiprocess scenarios.  See:
  https://github.com/PyTables/PyTables/issues/412

  The drawback of this quick fix is that the Blosc filter will be only
  able to use a single thread until another solution can be devised.


Changes from 1.5.0 to 1.5.1
===========================

* Updated to LZ4 1.5.0.  Closes #74.

* Added the 'const' qualifier to non SSE2 shuffle functions. Closes #75.

* Explicitly call blosc_init() in HDF5 blosc_filter.c, fixing a
  segfault.

* Quite a few improvements in cmake files for HDF5 support.  Thanks to
  Dana Robinson (The HDF Group).

* Variable 'class' caused problems compiling the HDF5 filter with g++.
  Thanks to Laurent Chapon.

* Small improvements on docstrings of c-blosc main functions.


Changes from 1.4.1 to 1.5.0
===========================

* Added new calls for allowing Blosc to be used *simultaneously*
  (i.e. lock free) from multi-threaded environments.  The new
  functions are:

  - blosc_compress_ctx(...)
  - blosc_decompress_ctx(...)

  See the new docstrings in blosc.h for how to use them.  The previous
  API should be completely unaffected.  Thanks to Christopher Speller.

* Optimized copies during BloscLZ decompression.  This can make BloscLZ
  to decompress up to 1.5x faster in some situations.

* LZ4 and LZ4HC compressors updated to version 1.3.1.

* Added an examples directory on how to link apps with Blosc.

* stdlib.h moved from blosc.c to blosc.h as suggested by Rob Lathm.

* Fix a warning for {snappy,lz4}-free compilation.  Thanks to Andrew Schaaf.

* Several improvements for CMakeLists.txt (cmake).

* Fixing C99 compatibility warnings.  Thanks to Christopher Speller.


Changes from 1.4.0 to 1.4.1
===========================

* Fixed a bug in blosc_getitem() introduced in 1.4.0.  Added a test for
  blosc_getitem() as well.


Changes from 1.3.6 to 1.4.0
===========================

* Support for non-Intel and non-SSE2 architectures has been added.  In
  particular, the Raspberry Pi platform (ARM) has been tested and all
  tests pass here.

* Architectures requiring strict access alignment are supported as well.
  Due to this, architectures with a high penalty in accessing unaligned
  data (e.g. Raspberry Pi, ARMv6) can compress up to 2.5x faster.

* LZ4 has been updated to r119 (1.2.0) so as to fix a possible security
  breach.


Changes from 1.3.5 to 1.3.6
===========================

* Updated to LZ4 r118 due to a (highly unlikely) security hole.  For
  details see:

  http://fastcompression.blogspot.fr/2014/06/debunking-lz4-20-years-old-bug-myth.html


Changes from 1.3.4 to 1.3.5
===========================

* Removed a pointer from 'pointer from integer without a cast' compiler
  warning due to a bad macro definition.


Changes from 1.3.3 to 1.3.4
===========================

* Fixed a false buffer overrun condition.  This bug made c-blosc to
  fail, even if the failure was not real.

* Fixed the type of a buffer string.


Changes from 1.3.2 to 1.3.3
===========================

* Updated to LZ4 1.1.3 (improved speed for 32-bit platforms).

* Added a new `blosc_cbuffer_complib()` for getting the compression
  library for a compressed buffer.


Changes from 1.3.1 to 1.3.2
===========================

* Fix for compiling Snappy sources against MSVC 2008.  Thanks to Mark
  Wiebe!

* Version for internal LZ4 and Snappy are now supported.  When compiled
  against the external libraries, this info is not available because
  they do not support the symbols (yet).


Changes from 1.3.0 to 1.3.1
===========================

* Fixes for a series of issues with the filter for HDF5 and, in
  particular, a problem in the decompression buffer size that made it
  impossible to use the blosc_filter in combination with other ones
  (e.g. fletcher32).  See
  https://github.com/PyTables/PyTables/issues/21.

  Thanks to Antonio Valentino for the fix!


Changes from 1.2.4 to 1.3.0
===========================

A nice handful of compressors have been added to Blosc:

* LZ4 (http://code.google.com/p/lz4/): A very fast
  compressor/decompressor.  Could be thought as a replacement of the
  original BloscLZ, but it can behave better is some scenarios.

* LZ4HC (http://code.google.com/p/lz4/): This is a variation of LZ4
  that achieves much better compression ratio at the cost of being
  much slower for compressing.  Decompression speed is unaffected (and
  sometimes better than when using LZ4 itself!), so this is very good
  for read-only datasets.

* Snappy (http://code.google.com/p/snappy/): A very fast
  compressor/decompressor.  Could be thought as a replacement of the
  original BloscLZ, but it can behave better is some scenarios.

* Zlib (http://www.zlib.net/): This is a classic.  It achieves very
  good compression ratios, at the cost of speed.  However,
  decompression speed is still pretty good, so it is a good candidate
  for read-only datasets.

With this, you can select the compression library with the new
function::

  int blosc_set_complib(char* complib);

where you pass the library that you want to use (currently "blosclz",
"lz4", "lz4hc", "snappy" and "zlib", but the list can grow in the
future).

You can get more info about compressors support in you Blosc build by
using these functions::

  char* blosc_list_compressors(void);
  int blosc_get_complib_info(char *compressor, char **complib, char **version);


Changes from 1.2.2 to 1.2.3
===========================

- Added a `blosc_init()` and `blosc_destroy()` so that the global lock
  can be initialized safely.  These new functions will also allow other
  kind of initializations/destructions in the future.

  Existing applications using Blosc do not need to start using the new
  functions right away, as long as they calling `blosc_set_nthreads()`
  previous to anything else.  However, using them is highly recommended.

  Thanks to Oscar Villellas for the init/destroy suggestion, it is a
  nice idea!


Changes from 1.2.1 to 1.2.2
===========================

- All important warnings removed for all tested platforms.  This will
  allow less intrusiveness compilation experiences with applications
  including Blosc source code.

- The `bench/bench.c` has been updated so that it can be compiled on
  Windows again.

- The new web site has been set to: http://www.blosc.org


Changes from 1.2 to 1.2.1
=========================

- Fixed a problem with global lock not being initialized.  This
  affected mostly to Windows platforms.  Thanks to Christoph
  Gohlke for finding the cure!


Changes from 1.1.5 to 1.2
=========================

- Now it is possible to call Blosc simultaneously from a parent threaded
  application without problems.  This has been solved by setting a
  global lock so that the different calling threads do not execute Blosc
  routines at the same time.  Of course, real threading work is still
  available *inside* Blosc itself.  Thanks to Thibault North.

- Support for cmake is now included.  Linux, Mac OSX and Windows
  platforms are supported.  Thanks to Thibault North, Antonio Valentino
  and Mark Wiebe.

- Fixed many compilers warnings (specially about unused variables).

- As a consequence of the above, as minimal change in the API has been
  introduced.  That is, the previous API::

    void blosc_free_resources(void)

  has changed to::

    int blosc_free_resources(void)

  Now, a return value of 0 means that the resources have been released
  successfully.  If the return value is negative, then it is not
  guaranteed that all the resources have been freed.

- Many typos were fixed and docs have been improved.  The script for
  generating nice plots for the included benchmarks has been improved
  too.  Thanks to Valetin Haenel.


Changes from 1.1.4 to 1.1.5
===========================

- Fix compile error with msvc compilers (Christoph Gohlke)


Changes from 1.1.3 to 1.1.4
===========================

- Redefinition of the BLOSC_MAX_BUFFERSIZE constant as (INT_MAX -
  BLOSC_MAX_OVERHEAD) instead of just INT_MAX.  This prevents to produce
  outputs larger than INT_MAX, which is not supported.

- `exit()` call has been replaced by a ``return -1`` in blosc_compress()
  when checking for buffer sizes.  Now programs will not just exit when
  the buffer is too large, but return a negative code.

- Improvements in explicit casts.  Blosc compiles without warnings
  (with GCC) now.

- Lots of improvements in docs, in particular a nice ascii-art diagram
  of the Blosc format (Valentin Haenel).

- Improvements to the plot-speeds.py (Valentin Haenel).

- [HDF5 filter] Adapted HDF5 filter to use HDF5 1.8 by default
  (Antonio Valentino).

- [HDF5 filter] New version of H5Z_class_t definition (Antonio Valentino).


Changes from 1.1.2 to 1.1.3
===========================

- Much improved compression ratio when using large blocks (> 64 KB) and
  high compression levels (> 6) under some circumstances (special data
  distribution).  Closes #7.


Changes from 1.1.1 to 1.1.2
===========================

- Fixes for small typesizes (#6 and #1 of python-blosc).


Changes from 1.1 to 1.1.1
=========================

- Added code to avoid calling blosc_set_nthreads more than necessary.
  That will improve performance up to 3x or more, specially for small
  chunksizes (< 1 MB).


Changes from 1.0 to 1.1
=======================

- Added code for emulating pthreads API on Windows.  No need to link
  explicitly with pthreads lib on Windows anymore.  However, performance
  is a somewhat worse because the new emulation layer does not support
  the `pthread_barrier_wait()` call natively.  But the big improvement
  in installation easiness is worth this penalty (most specially on
  64-bit Windows, where pthreads-win32 support is flaky).

- New BLOSC_MAX_BUFFERSIZE, BLOSC_MAX_TYPESIZE and BLOSC_MAX_THREADS
  symbols are available in blosc.h.  These can be useful for validating
  parameters in clients.  Thanks to Robert Smallshire for suggesting
  that.

- A new BLOSC_MIN_HEADER_LENGTH symbol in blosc.h tells how many bytes
  long is the minimum length of a Blosc header.  `blosc_cbuffer_sizes()`
  only needs these bytes to be passed to work correctly.

- Removed many warnings (related with potentially dangerous type-casting
  code) issued by MSVC 2008 in 64-bit mode.

- Fixed a problem with the computation of the blocksize in the Blosc
  filter for HDF5.

- Fixed a problem with large datatypes.  See
  http://www.pytables.org/trac/ticket/288 for more info.

- Now Blosc is able to work well even if you fork an existing process
  with a pool of threads.  Bug discovered when PyTables runs in
  multiprocess environments.  See http://pytables.org/trac/ticket/295
  for details.

- Added a new `blosc_getitem()` call to allow the retrieval of items in
  sizes smaller than the complete buffer.  That is useful for the carray
  project, but certainly for others too.


Changes from 0.9.5 to 1.0
=========================

- Added a filter for HDF5 so that people can use Blosc outside PyTables,
  if they want to.

- Many small improvements, specially in README files.

- Do not assume that size_t is uint_32 for every platform.

- Added more protection for large buffers or in allocation memory
  routines.

- The src/ directory has been renamed to blosc/.

- The `maxbytes` parameter in `blosc_compress()` has been renamed to
  `destsize`.  This is for consistency with the `blosc_decompress()`
  parameters.


Changes from 0.9.4 to 0.9.5
===========================

- Now, compression level 0 is allowed, meaning not compression at all.
  The overhead of this mode will be always BLOSC_MAX_OVERHEAD (16)
  bytes.  This mode actually represents using Blosc as a basic memory
  container.

- Supported a new parameter `maxbytes` for ``blosc_compress()``.  It
  represents a maximum of bytes for output.  Tests unit added too.

- Added 3 new functions for querying different metadata on compressed
  buffers.  A test suite for testing the new API has been added too.


Changes from 0.9.3 to 0.9.4
===========================

- Support for cross-platform big/little endian compatibility in Blosc
  headers has been added.

- Fixed several failures exposed by the extremesuite.  The problem was a
  bad check for limits in the buffer size while compressing.

- Added a new suite in bench.c called ``debugsuite`` that is
  appropriate for debugging purposes.  Now, the ``extremesuite`` can be
  used for running the complete (and extremely long) suite.


Changes from 0.9.0 to 0.9.3
===========================

- Fixed several nasty bugs uncovered by the new suites in bench.c.
  Thanks to Tony Theodore and Gabriel Beckers for their (very)
  responsive beta testing and feedback.

- Added several modes (suites), namely ``suite``, ``hardsuite`` and
  ``extremehardsuite`` in bench.c so as to allow different levels of
  testing.


Changes from 0.8.0 to 0.9
=========================

- Internal format version bumped to 2 in order to allow an easy way to
  indicate that a buffer is being saved uncompressed.  This is not
  supported yet, but it might be in the future.

- Blosc can use threads now for leveraging the increasing number of
  multi-core processors out there.  See README-threaded.txt for more
  info.

- Added a protection for MacOSX so that it has to not link against
  posix_memalign() function, which seems not available in old versions of
  MacOSX (for example, Tiger).  At nay rate, posix_memalign() is not
  necessary on Mac because 16 bytes alignment is ensured by default.
  Thanks to Ivan Vilata.  Fixes #3.
================
Releasing Blosc
================

:Author: Francesc Alted
:Contact: francesc@blosc.org
:Date: 2014-01-15


Preliminaries
-------------

- Switch to master branch::

    $ git switch master

- Make sure that ``RELEASE_NOTES.rst`` and ``ANNOUNCE.rst`` are up to
  date with the latest news in the release.

- Check that *VERSION* symbols in blosc/blosc.h contains the correct info.

- Commit the changes::

    $ git commit -a -m"Getting ready for X.Y.Z release"


Testing
-------

Create a new build/ directory, change into it and issue::

  $ cmake ..
  $ cmake --build .
  $ ctest

To actually test Blosc the hard way, look at the end of:

http://blosc.org/synthetic-benchmarks.html

where instructions on how to intensively test (and benchmark) Blosc
are given.

Forward compatibility testing
-----------------------------

First, go to the compat/ directory and generate a file with the current
version::

  $ cd ../compat
  $ export LD_LIBRARY_PATH=../build/blosc
  $ gcc -o filegen filegen.c -L$LD_LIBRARY_PATH -lblosc -I../blosc
  $ ./filegen compress lz4 blosc-1.y.z-lz4.cdata

In order to make sure that we are not breaking forward compatibility,
link and run the `compat/filegen` utility against different versions of
the Blosc library (suggestion: 1.3.0, 1.7.0, 1.11.1, 1.14.1).

You can compile the utility with different blosc shared libraries with::

  $ export LD_LIBRARY_PATH=shared_blosc_library_path
  $ gcc -o filegen filegen.c -L$LD_LIBRARY_PATH -lblosc -Iblosc.h_include_path

Then, test the file created with the new version with::

  $ ./filegen decompress blosc-1.y.z-lz4.cdata

If that works and you want to keep track of this for future compatibility checks
just add the new file to the suite::

  $ git add blosc-1.y.z-lz4.cdata
  $ git commit -m"Add a new cdata file for compatibility checks"

Repeat this for every codec shipped with Blosc (blosclz, lz4, lz4hc, snappy,
zlib and zstd).

Tagging
-------

- Create a tag ``X.Y.Z`` from ``master``::

    $ git tag -a vX.Y.Z -m "Tagging version X.Y.Z"

- Push the previous commits and tag to the github repo::

    $ git push
    $ git push --tags


Announcing
----------

- Send an announcement to the blosc, pytables-dev, bcolz and
  comp.compression lists.  Use the ``ANNOUNCE.rst`` file as skeleton
  (possibly as the definitive version).


Post-release actions
--------------------

- Edit *VERSION* symbols in blosc/blosc.h in master to increment the
  version to the next minor one (i.e. X.Y.Z --> X.Y.(Z+1).dev).

- Create new headers for adding new features in ``RELEASE_NOTES.rst``
  and add this place-holder instead:

  #XXX version-specific blurb XXX#

- Commit the changes::

    $ git commit -a -m"Post X.Y.Z release actions done"
    $ git push


That's all folks!


.. Local Variables:
.. mode: rst
.. coding: utf-8
.. fill-column: 70
.. End:
Blosc supports threading
========================

Threads are the most efficient way to program parallel code for
multi-core processors, but also the more difficult to program well.
Also, they has a non-negligible start-up time that does not fit well
with a high-performance compressor as Blosc tries to be.

In order to reduce the overhead of threads as much as possible, I've
decided to implement a pool of threads (the workers) that are waiting
for the main process (the master) to send them jobs (basically,
compressing and decompressing small blocks of the initial buffer).

Despite this and many other internal optimizations in the threaded
code, it does not work faster than the serial version for buffer sizes
around 64/128 KB or less.  This is for Intel Quad Core2 (Q8400 @ 2.66
GHz) / Linux (openSUSE 11.2, 64 bit), but your mileage may vary (and
will vary!) for other processors / operating systems.

In contrast, for buffers larger than 64/128 KB, the threaded version
starts to perform significantly better, being the sweet point at 1 MB
(again, this is with my setup).  For larger buffer sizes than 1 MB,
the threaded code slows down again, but it is probably due to a cache
size issue and besides, it is still considerably faster than serial
code.

This is why Blosc falls back to use the serial version for such a
'small' buffers.  So, you don't have to worry too much about deciding
whether you should set the number of threads to 1 (serial) or more
(parallel).  Just set it to the number of cores in your processor and
your are done!

Francesc Alted
Blosc Chunk Format
==================

The chunk is composed by a header and a blocks / splits section::

    +---------+--------+---------+
    |  header | blocks / splits  |
    +---------+--------+---------+

These are described below.

The header section
------------------

Blosc (as of Version 1.0.0) has the following 16 byte header that stores
information about the compressed buffer::

    |-0-|-1-|-2-|-3-|-4-|-5-|-6-|-7-|-8-|-9-|-A-|-B-|-C-|-D-|-E-|-F-|
      ^   ^   ^   ^ |     nbytes    |   blocksize   |    cbytes     |
      |   |   |   |
      |   |   |   +--typesize
      |   |   +------flags
      |   +----------versionlz
      +--------------version

Datatypes of the header entries
-------------------------------

All entries are little endian.

:version:
    (``uint8``) Blosc format version.
:versionlz:
    (``uint8``) Version of the internal compressor used.
:flags and compressor enumeration:
    (``bitfield``) The flags of the buffer

    :bit 0 (``0x01``):
        Whether the byte-shuffle filter has been applied or not.
    :bit 1 (``0x02``):
        Whether the internal buffer is a pure memcpy or not.
    :bit 2 (``0x04``):
        Whether the bit-shuffle filter has been applied or not.
    :bit 3 (``0x08``):
        Reserved, must be zero.
    :bit 4 (``0x10``):
        If set, the blocks will not be split in sub-blocks during compression.
    :bit 5 (``0x20``):
        Part of the enumeration for compressors.
    :bit 6 (``0x40``):
        Part of the enumeration for compressors.
    :bit 7 (``0x80``):
        Part of the enumeration for compressors.

    The last three bits form an enumeration that allows to use alternative
    compressors.

    :``0``:
        ``blosclz``
    :``1``:
        ``lz4`` or ``lz4hc``
    :``2``:
        ``snappy``
    :``3``:
        ``zlib``
    :``4``:
        ``zstd``

:typesize:
    (``uint8``) Number of bytes for the atomic type.
:nbytes:
    (``uint32``) Uncompressed size of the buffer (this header is not included).
:blocksize:
    (``uint32``) Size of internal blocks.
:cbytes:
    (``uint32``) Compressed size of the buffer (including this header).

The blocks / splits section
---------------------------

After the header, there come the blocks / splits section.  Blocks are equal-sized parts of the chunk, except for the last block that can be shorter or equal than the rest.

At the beginning of the blocks section, there come a list of `int32_t bstarts` to indicate where the different encoded blocks starts (counting from the end of this `bstarts` section)::

    +=========+=========+========+=========+
    | bstart0 | bstart1 |   ...  | bstartN |
    +=========+=========+========+=========+

Finally, it comes the actual list of compressed blocks / splits data streams.  It turns out that a block may optionally (see bit 4 in `flags` above) be further split in so-called splits which are the actual data streams that are transmitted to codecs for compression.  If a block is not split, then the split is equivalent to a whole block.  Before each split in the list, there is the compressed size of it, expressed as an `int32_t`::

    +========+========+========+========+========+========+========+
    | csize0 | split0 | csize1 | split1 |   ...  | csizeN | splitN |
    +========+========+========+========+========+========+========+


*Note*: all the integers are stored in little endian.

I'd like to thank the PyTables community that have collaborated in the
exhaustive testing of Blosc.  With an aggregate amount of more than
300 TB of different datasets compressed *and* decompressed
successfully, I can say that Blosc is pretty safe now and ready for
production purposes.

Other important contributions:

* Valentin Haenel did a terrific work implementing the support for the
  Snappy compression, fixing typos and improving docs and the plotting
  script.

* Thibault North, with ideas from Oscar Villellas, contributed a way
  to call Blosc from different threads in a safe way.  Christopher
  Speller introduced contexts so that a global lock is not necessary
  anymore.

* The CMake support was initially contributed by Thibault North, and
  Antonio Valentino and Mark Wiebe made great enhancements to it.

* Christopher Speller also introduced the two new '_ctx' calls to
  avoid the use of the blosc_init() and blosc_destroy().

* Jack Pappas contributed important portability enhancements,
  specially runtime and cross-platform detection of SSE2/AVX2 as well
  as high precision timers (HPET) for the benchmark program.

* @littlezhou implemented the AVX2 version of shuffle routines.

* Julian Taylor contributed a way to detect AVX2 in runtime and
  calling the appropriate routines only if the underlying hardware
  supports it.

* Kiyo Masui for relicensing his bitshuffle project for allowing the
  inclusion of part of his code in Blosc.
===============================================================
 Announcing C-Blosc 1.21.1
 A blocking, shuffling and lossless compression library for C
===============================================================

What is new?
============

This is a maintenance release.  Fix pthread flag when linking on ppc64le.
Vendored BloscLZ, Zlib and Zstd codecs have been updated to their latest
versions too; this can bring important performance improvements, so if
speed is a priority to you, an upgrade is recommended.

For more info, please see the release notes in:

https://github.com/Blosc/c-blosc/blob/master/RELEASE_NOTES.rst


What is it?
===========

Blosc (http://www.blosc.org) is a high performance meta-compressor
optimized for binary data.  It has been designed to transmit data to
the processor cache faster than the traditional, non-compressed,
direct memory fetch approach via a memcpy() OS call.

Blosc has internal support for different compressors like its internal
BloscLZ, but also LZ4, LZ4HC, Snappy, Zlib and Zstd.  This way these can
automatically leverage the multithreading and pre-filtering
(shuffling) capabilities that comes with Blosc.


Download sources
================

The github repository is over here:

https://github.com/Blosc

Blosc is distributed using the BSD license, see LICENSES/BLOSC.txt for
details.


Mailing list
============

There is an official Blosc mailing list at:

blosc@googlegroups.com
http://groups.google.es/group/blosc


Enjoy Data!
Compressed datafiles for testing backward/forward compatibility
===============================================================

The files here have been created with different versions of the C-Blosc library and are meant to test backward/forward compatibility among different versions of the library.
Examples on how to add Blosc support for your programs
======================================================

In this directory you can find a series of examples on how to link
your apps with the Blosc library:

* simple.c -- The simplest way to add Blosc to your app
* multithread.c -- Add multithreading into the equation
* many_compressors.c -- Use different compressors inside Blosc

For more info, please visit the `official API documentation
<https://github.com/Blosc/c-blosc/blob/master/blosc/blosc.h>`_.
