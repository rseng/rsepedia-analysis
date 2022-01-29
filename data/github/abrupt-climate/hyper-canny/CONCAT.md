[![Build Status](https://travis-ci.org/abrupt-climate/hyper-canny.svg?branch=master)](https://travis-ci.org/abrupt-climate/hyper-canny)
[![codecov](https://codecov.io/gh/abrupt-climate/hyper-canny/branch/master/graph/badge.svg)](https://codecov.io/gh/abrupt-climate/hyper-canny)
[![DOI](https://zenodo.org/badge/86834996.svg)](https://zenodo.org/badge/latestdoi/86834996)

HyperCanny project
==================

Doing Canny edge detection in higher dimensions.

Requirements
------------

HyperCanny is written in C++17, using NetCDF as data container. Dependencies:

- GCC >= 6 or CLANG >= 4.0
- [NetCDF-cxx4](https://github.com/Unidata/netcdf-cxx4)

Building
--------

HyperCanny uses [Meson](http://mesonbuild.com/), and [Ninja](https://ninja-build.org/). To build HyperCanny, run

```bash
> meson build           # you only need to run this once
> ninja -C build
```

To run the unit tests, follow up with `ninja test` in the `build` directory.

Alternatively, there is a `Makefile` that runs these commands for you; this also allows for easier integration with some IDEs.

Python module
-------------

The code also compiles into a Python module. From the project root run:

```bash
> pip install .
```

Testing / Development
---------------------

This project uses [Google Test](https://github.com/google/googletest) for testing. The source code of `gtest` is included in the repository, so you don't need to install it separately.

### coverage
To run tests with coverage you may want to use the `scripts/coverage.sh` script. It recompiles the code with `--coverage`, runs tests and then either generates a HTML report using the `lcov` tool (`./scripts/coverage.sh -html`), or sends the coverage report to [codecov](https://codecov.io/) (`./scripts/coverage.sh -codecov ${CODECOV_TOKEN}`).

### code quality
- [CppCheck](http://cppcheck.sourceforge.net/) is run through the `scripts/run_cppcheck.sh` script.
- To check that all headers start with license info and a `#pragma once` directive, run the `scripts/check-headers.sh` script. This also checks that all headers have a Doxygen `\file` directive.

### development
- If you add `.cc` files somewhere, be sure to add them to the appropriate `meson.build` file. This is done automatically using the `scripts/premeson.sh` script.

## Google Mock ##

The Google C++ mocking framework.

### Overview ###

Google's framework for writing and using C++ mock classes.
It can help you derive better designs of your system and write better tests.

It is inspired by:

  * [jMock](http://www.jmock.org/),
  * [EasyMock](http://www.easymock.org/), and
  * [Hamcrest](http://code.google.com/p/hamcrest/),

and designed with C++'s specifics in mind.

Google mock:

  * lets you create mock classes trivially using simple macros.
  * supports a rich set of matchers and actions.
  * handles unordered, partially ordered, or completely ordered expectations.
  * is extensible by users.

We hope you find it useful!

### Features ###

  * Provides a declarative syntax for defining mocks.
  * Can easily define partial (hybrid) mocks, which are a cross of real
    and mock objects.
  * Handles functions of arbitrary types and overloaded functions.
  * Comes with a rich set of matchers for validating function arguments.
  * Uses an intuitive syntax for controlling the behavior of a mock.
  * Does automatic verification of expectations (no record-and-replay needed).
  * Allows arbitrary (partial) ordering constraints on
    function calls to be expressed,.
  * Lets a user extend it by defining new matchers and actions.
  * Does not use exceptions.
  * Is easy to learn and use.

Please see the project page above for more information as well as the
mailing list for questions, discussions, and development.  There is
also an IRC channel on OFTC (irc.oftc.net) #gtest available.  Please
join us!

Please note that code under [scripts/generator](scripts/generator/) is
from [cppclean](http://code.google.com/p/cppclean/) and released under
the Apache License, which is different from Google Mock's license.

## Getting Started ##

If you are new to the project, we suggest that you read the user
documentation in the following order:

  * Learn the [basics](../../master/googletest/docs/Primer.md) of
    Google Test, if you choose to use Google Mock with it (recommended).
  * Read [Google Mock for Dummies](../../master/googlemock/docs/ForDummies.md).
  * Read the instructions below on how to build Google Mock.

You can also watch Zhanyong's [talk](http://www.youtube.com/watch?v=sYpCyLI47rM) on Google Mock's usage and implementation.

Once you understand the basics, check out the rest of the docs:

  * [CheatSheet](../../master/googlemock/docs/CheatSheet.md) - all the commonly used stuff
    at a glance.
  * [CookBook](../../master/googlemock/docs/CookBook.md) - recipes for getting things done,
    including advanced techniques.

If you need help, please check the
[KnownIssues](docs/KnownIssues.md) and
[FrequentlyAskedQuestions](docs/FrequentlyAskedQuestions.md) before
posting a question on the
[discussion group](http://groups.google.com/group/googlemock).


### Using Google Mock Without Google Test ###

Google Mock is not a testing framework itself.  Instead, it needs a
testing framework for writing tests.  Google Mock works seamlessly
with [Google Test](http://code.google.com/p/googletest/), but
you can also use it with [any C++ testing framework](../../master/googlemock/docs/ForDummies.md#using-google-mock-with-any-testing-framework).

### Requirements for End Users ###

Google Mock is implemented on top of [Google Test](
http://github.com/google/googletest/), and depends on it.
You must use the bundled version of Google Test when using Google Mock.

You can also easily configure Google Mock to work with another testing
framework, although it will still need Google Test.  Please read
["Using_Google_Mock_with_Any_Testing_Framework"](
    ../../master/googlemock/docs/ForDummies.md#using-google-mock-with-any-testing-framework)
for instructions.

Google Mock depends on advanced C++ features and thus requires a more
modern compiler. The following are needed to use Google Mock:

#### Linux Requirements ####

  * GNU-compatible Make or "gmake"
  * POSIX-standard shell
  * POSIX(-2) Regular Expressions (regex.h)
  * C++98-standard-compliant compiler (e.g. GCC 3.4 or newer)

#### Windows Requirements ####

  * Microsoft Visual C++ 8.0 SP1 or newer

#### Mac OS X Requirements ####

  * Mac OS X 10.4 Tiger or newer
  * Developer Tools Installed

### Requirements for Contributors ###

We welcome patches. If you plan to contribute a patch, you need to
build Google Mock and its tests, which has further requirements:

  * Automake version 1.9 or newer
  * Autoconf version 2.59 or newer
  * Libtool / Libtoolize
  * Python version 2.3 or newer (for running some of the tests and
    re-generating certain source files from templates)

### Building Google Mock ###

If you have CMake available, it is recommended that you follow the
[build instructions][gtest_cmakebuild]
as described for Google Test. If are using Google Mock with an
existing CMake project, the section
[Incorporating Into An Existing CMake Project][gtest_incorpcmake]
may be of particular interest. Otherwise, the following sections
detail how to build Google Mock without CMake.

#### Preparing to Build (Unix only) ####

If you are using a Unix system and plan to use the GNU Autotools build
system to build Google Mock (described below), you'll need to
configure it now.

To prepare the Autotools build system:

    cd googlemock
    autoreconf -fvi

To build Google Mock and your tests that use it, you need to tell your
build system where to find its headers and source files.  The exact
way to do it depends on which build system you use, and is usually
straightforward.

This section shows how you can integrate Google Mock into your
existing build system.

Suppose you put Google Mock in directory `${GMOCK_DIR}` and Google Test
in `${GTEST_DIR}` (the latter is `${GMOCK_DIR}/gtest` by default).  To
build Google Mock, create a library build target (or a project as
called by Visual Studio and Xcode) to compile

    ${GTEST_DIR}/src/gtest-all.cc and ${GMOCK_DIR}/src/gmock-all.cc

with

    ${GTEST_DIR}/include and ${GMOCK_DIR}/include

in the system header search path, and

    ${GTEST_DIR} and ${GMOCK_DIR}

in the normal header search path.  Assuming a Linux-like system and gcc,
something like the following will do:

    g++ -isystem ${GTEST_DIR}/include -I${GTEST_DIR} \
        -isystem ${GMOCK_DIR}/include -I${GMOCK_DIR} \
        -pthread -c ${GTEST_DIR}/src/gtest-all.cc
    g++ -isystem ${GTEST_DIR}/include -I${GTEST_DIR} \
        -isystem ${GMOCK_DIR}/include -I${GMOCK_DIR} \
        -pthread -c ${GMOCK_DIR}/src/gmock-all.cc
    ar -rv libgmock.a gtest-all.o gmock-all.o

(We need -pthread as Google Test and Google Mock use threads.)

Next, you should compile your test source file with
${GTEST\_DIR}/include and ${GMOCK\_DIR}/include in the header search
path, and link it with gmock and any other necessary libraries:

    g++ -isystem ${GTEST_DIR}/include -isystem ${GMOCK_DIR}/include \
        -pthread path/to/your_test.cc libgmock.a -o your_test

As an example, the make/ directory contains a Makefile that you can
use to build Google Mock on systems where GNU make is available
(e.g. Linux, Mac OS X, and Cygwin).  It doesn't try to build Google
Mock's own tests.  Instead, it just builds the Google Mock library and
a sample test.  You can use it as a starting point for your own build
script.

If the default settings are correct for your environment, the
following commands should succeed:

    cd ${GMOCK_DIR}/make
    make
    ./gmock_test

If you see errors, try to tweak the contents of
[make/Makefile](make/Makefile) to make them go away.

### Windows ###

The msvc/2005 directory contains VC++ 2005 projects and the msvc/2010
directory contains VC++ 2010 projects for building Google Mock and
selected tests.

Change to the appropriate directory and run "msbuild gmock.sln" to
build the library and tests (or open the gmock.sln in the MSVC IDE).
If you want to create your own project to use with Google Mock, you'll
have to configure it to use the `gmock_config` propety sheet.  For that:

 * Open the Property Manager window (View | Other Windows | Property Manager)
 * Right-click on your project and select "Add Existing Property Sheet..."
 * Navigate to `gmock_config.vsprops` or `gmock_config.props` and select it.
 * In Project Properties | Configuration Properties | General | Additional
   Include Directories, type <path to Google Mock>/include.

### Tweaking Google Mock ###

Google Mock can be used in diverse environments.  The default
configuration may not work (or may not work well) out of the box in
some environments.  However, you can easily tweak Google Mock by
defining control macros on the compiler command line.  Generally,
these macros are named like `GTEST_XYZ` and you define them to either 1
or 0 to enable or disable a certain feature.

We list the most frequently used macros below.  For a complete list,
see file [${GTEST\_DIR}/include/gtest/internal/gtest-port.h](
../googletest/include/gtest/internal/gtest-port.h).

### Choosing a TR1 Tuple Library ###

Google Mock uses the C++ Technical Report 1 (TR1) tuple library
heavily.  Unfortunately TR1 tuple is not yet widely available with all
compilers.  The good news is that Google Test 1.4.0+ implements a
subset of TR1 tuple that's enough for Google Mock's need.  Google Mock
will automatically use that implementation when the compiler doesn't
provide TR1 tuple.

Usually you don't need to care about which tuple library Google Test
and Google Mock use.  However, if your project already uses TR1 tuple,
you need to tell Google Test and Google Mock to use the same TR1 tuple
library the rest of your project uses, or the two tuple
implementations will clash.  To do that, add

    -DGTEST_USE_OWN_TR1_TUPLE=0

to the compiler flags while compiling Google Test, Google Mock, and
your tests.  If you want to force Google Test and Google Mock to use
their own tuple library, just add

    -DGTEST_USE_OWN_TR1_TUPLE=1

to the compiler flags instead.

If you want to use Boost's TR1 tuple library with Google Mock, please
refer to the Boost website (http://www.boost.org/) for how to obtain
it and set it up.

### As a Shared Library (DLL) ###

Google Mock is compact, so most users can build and link it as a static
library for the simplicity.  Google Mock can be used as a DLL, but the
same DLL must contain Google Test as well.  See
[Google Test's README][gtest_readme]
for instructions on how to set up necessary compiler settings.

### Tweaking Google Mock ###

Most of Google Test's control macros apply to Google Mock as well.
Please see [Google Test's README][gtest_readme] for how to tweak them.

### Upgrading from an Earlier Version ###

We strive to keep Google Mock releases backward compatible.
Sometimes, though, we have to make some breaking changes for the
users' long-term benefits.  This section describes what you'll need to
do if you are upgrading from an earlier version of Google Mock.

#### Upgrading from 1.1.0 or Earlier ####

You may need to explicitly enable or disable Google Test's own TR1
tuple library.  See the instructions in section "[Choosing a TR1 Tuple
Library](../googletest/#choosing-a-tr1-tuple-library)".

#### Upgrading from 1.4.0 or Earlier ####

On platforms where the pthread library is available, Google Test and
Google Mock use it in order to be thread-safe.  For this to work, you
may need to tweak your compiler and/or linker flags.  Please see the
"[Multi-threaded Tests](../googletest#multi-threaded-tests
)" section in file Google Test's README for what you may need to do.

If you have custom matchers defined using `MatcherInterface` or
`MakePolymorphicMatcher()`, you'll need to update their definitions to
use the new matcher API (
[monomorphic](http://code.google.com/p/googlemock/wiki/CookBook#Writing_New_Monomorphic_Matchers),
[polymorphic](http://code.google.com/p/googlemock/wiki/CookBook#Writing_New_Polymorphic_Matchers)).
Matchers defined using `MATCHER()` or `MATCHER_P*()` aren't affected.

### Developing Google Mock ###

This section discusses how to make your own changes to Google Mock.

#### Testing Google Mock Itself ####

To make sure your changes work as intended and don't break existing
functionality, you'll want to compile and run Google Test's own tests.
For that you'll need Autotools.  First, make sure you have followed
the instructions above to configure Google Mock.
Then, create a build output directory and enter it.  Next,

    ${GMOCK_DIR}/configure  # try --help for more info

Once you have successfully configured Google Mock, the build steps are
standard for GNU-style OSS packages.

    make        # Standard makefile following GNU conventions
    make check  # Builds and runs all tests - all should pass.

Note that when building your project against Google Mock, you are building
against Google Test as well.  There is no need to configure Google Test
separately.

#### Contributing a Patch ####

We welcome patches.
Please read the [Developer's Guide](docs/DevGuide.md)
for how you can contribute. In particular, make sure you have signed
the Contributor License Agreement, or we won't be able to accept the
patch.

Happy testing!

[gtest_readme]: ../googletest/README.md "googletest"
[gtest_cmakebuild]:  ../googletest/README.md#using-cmake "Using CMake"
[gtest_incorpcmake]: ../googletest/README.md#incorporating-into-an-existing-cmake-project "Incorporating Into An Existing CMake Project"

### Generic Build Instructions ###

#### Setup ####

To build Google Test and your tests that use it, you need to tell your
build system where to find its headers and source files.  The exact
way to do it depends on which build system you use, and is usually
straightforward.

#### Build ####

Suppose you put Google Test in directory `${GTEST_DIR}`.  To build it,
create a library build target (or a project as called by Visual Studio
and Xcode) to compile

    ${GTEST_DIR}/src/gtest-all.cc

with `${GTEST_DIR}/include` in the system header search path and `${GTEST_DIR}`
in the normal header search path.  Assuming a Linux-like system and gcc,
something like the following will do:

    g++ -isystem ${GTEST_DIR}/include -I${GTEST_DIR} \
        -pthread -c ${GTEST_DIR}/src/gtest-all.cc
    ar -rv libgtest.a gtest-all.o

(We need `-pthread` as Google Test uses threads.)

Next, you should compile your test source file with
`${GTEST_DIR}/include` in the system header search path, and link it
with gtest and any other necessary libraries:

    g++ -isystem ${GTEST_DIR}/include -pthread path/to/your_test.cc libgtest.a \
        -o your_test

As an example, the make/ directory contains a Makefile that you can
use to build Google Test on systems where GNU make is available
(e.g. Linux, Mac OS X, and Cygwin).  It doesn't try to build Google
Test's own tests.  Instead, it just builds the Google Test library and
a sample test.  You can use it as a starting point for your own build
script.

If the default settings are correct for your environment, the
following commands should succeed:

    cd ${GTEST_DIR}/make
    make
    ./sample1_unittest

If you see errors, try to tweak the contents of `make/Makefile` to make
them go away.  There are instructions in `make/Makefile` on how to do
it.

### Using CMake ###

Google Test comes with a CMake build script (
[CMakeLists.txt](CMakeLists.txt)) that can be used on a wide range of platforms ("C" stands for
cross-platform.). If you don't have CMake installed already, you can
download it for free from <http://www.cmake.org/>.

CMake works by generating native makefiles or build projects that can
be used in the compiler environment of your choice.  You can either
build Google Test as a standalone project or it can be incorporated
into an existing CMake build for another project.

#### Standalone CMake Project ####

When building Google Test as a standalone project, the typical
workflow starts with:

    mkdir mybuild       # Create a directory to hold the build output.
    cd mybuild
    cmake ${GTEST_DIR}  # Generate native build scripts.

If you want to build Google Test's samples, you should replace the
last command with

    cmake -Dgtest_build_samples=ON ${GTEST_DIR}

If you are on a \*nix system, you should now see a Makefile in the
current directory.  Just type 'make' to build gtest.

If you use Windows and have Visual Studio installed, a `gtest.sln` file
and several `.vcproj` files will be created.  You can then build them
using Visual Studio.

On Mac OS X with Xcode installed, a `.xcodeproj` file will be generated.

#### Incorporating Into An Existing CMake Project ####

If you want to use gtest in a project which already uses CMake, then a
more robust and flexible approach is to build gtest as part of that
project directly. This is done by making the GoogleTest source code
available to the main build and adding it using CMake's
`add_subdirectory()` command. This has the significant advantage that
the same compiler and linker settings are used between gtest and the
rest of your project, so issues associated with using incompatible
libraries (eg debug/release), etc. are avoided. This is particularly
useful on Windows. Making GoogleTest's source code available to the
main build can be done a few different ways:

* Download the GoogleTest source code manually and place it at a
  known location. This is the least flexible approach and can make
  it more difficult to use with continuous integration systems, etc.
* Embed the GoogleTest source code as a direct copy in the main
  project's source tree. This is often the simplest approach, but is
  also the hardest to keep up to date. Some organizations may not
  permit this method.
* Add GoogleTest as a git submodule or equivalent. This may not
  always be possible or appropriate. Git submodules, for example,
  have their own set of advantages and drawbacks.
* Use CMake to download GoogleTest as part of the build's configure
  step. This is just a little more complex, but doesn't have the
  limitations of the other methods.

The last of the above methods is implemented with a small piece
of CMake code in a separate file (e.g. `CMakeLists.txt.in`) which
is copied to the build area and then invoked as a sub-build
_during the CMake stage_. That directory is then pulled into the
main build with `add_subdirectory()`. For example:

New file `CMakeLists.txt.in`:

    cmake_minimum_required(VERSION 2.8.2)
 
    project(googletest-download NONE)
 
    include(ExternalProject)
    ExternalProject_Add(googletest
      GIT_REPOSITORY    https://github.com/google/googletest.git
      GIT_TAG           master
      SOURCE_DIR        "${CMAKE_BINARY_DIR}/googletest-src"
      BINARY_DIR        "${CMAKE_BINARY_DIR}/googletest-build"
      CONFIGURE_COMMAND ""
      BUILD_COMMAND     ""
      INSTALL_COMMAND   ""
      TEST_COMMAND      ""
    )
    
Existing build's `CMakeLists.txt`:

    # Download and unpack googletest at configure time
    configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
    execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
      RESULT_VARIABLE result
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
    if(result)
      message(FATAL_ERROR "CMake step for googletest failed: ${result}")
    endif()
    execute_process(COMMAND ${CMAKE_COMMAND} --build .
      RESULT_VARIABLE result
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
    if(result)
      message(FATAL_ERROR "Build step for googletest failed: ${result}")
    endif()

    # Prevent overriding the parent project's compiler/linker
    # settings on Windows
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    
    # Add googletest directly to our build. This defines
    # the gtest and gtest_main targets.
    add_subdirectory(${CMAKE_BINARY_DIR}/googletest-src
                     ${CMAKE_BINARY_DIR}/googletest-build)

    # The gtest/gtest_main targets carry header search path
    # dependencies automatically when using CMake 2.8.11 or
    # later. Otherwise we have to add them here ourselves.
    if (CMAKE_VERSION VERSION_LESS 2.8.11)
      include_directories("${gtest_SOURCE_DIR}/include")
    endif()

    # Now simply link against gtest or gtest_main as needed. Eg
    add_executable(example example.cpp)
    target_link_libraries(example gtest_main)
    add_test(NAME example_test COMMAND example)

Note that this approach requires CMake 2.8.2 or later due to
its use of the `ExternalProject_Add()` command. The above
technique is discussed in more detail in 
[this separate article](http://crascit.com/2015/07/25/cmake-gtest/)
which also contains a link to a fully generalized implementation
of the technique.


### Legacy Build Scripts ###

Before settling on CMake, we have been providing hand-maintained build
projects/scripts for Visual Studio, Xcode, and Autotools.  While we
continue to provide them for convenience, they are not actively
maintained any more.  We highly recommend that you follow the
instructions in the above sections to integrate Google Test
with your existing build system.

If you still need to use the legacy build scripts, here's how:

The msvc\ folder contains two solutions with Visual C++ projects.
Open the `gtest.sln` or `gtest-md.sln` file using Visual Studio, and you
are ready to build Google Test the same way you build any Visual
Studio project.  Files that have names ending with -md use DLL
versions of Microsoft runtime libraries (the /MD or the /MDd compiler
option).  Files without that suffix use static versions of the runtime
libraries (the /MT or the /MTd option).  Please note that one must use
the same option to compile both gtest and the test code.  If you use
Visual Studio 2005 or above, we recommend the -md version as /MD is
the default for new projects in these versions of Visual Studio.

On Mac OS X, open the `gtest.xcodeproj` in the `xcode/` folder using
Xcode.  Build the "gtest" target.  The universal binary framework will
end up in your selected build directory (selected in the Xcode
"Preferences..." -> "Building" pane and defaults to xcode/build).
Alternatively, at the command line, enter:

    xcodebuild

This will build the "Release" configuration of gtest.framework in your
default build location.  See the "xcodebuild" man page for more
information about building different configurations and building in
different locations.

If you wish to use the Google Test Xcode project with Xcode 4.x and
above, you need to either:

 * update the SDK configuration options in xcode/Config/General.xconfig.
   Comment options `SDKROOT`, `MACOS_DEPLOYMENT_TARGET`, and `GCC_VERSION`. If
   you choose this route you lose the ability to target earlier versions
   of MacOS X.
 * Install an SDK for an earlier version. This doesn't appear to be
   supported by Apple, but has been reported to work
   (http://stackoverflow.com/questions/5378518).

### Tweaking Google Test ###

Google Test can be used in diverse environments.  The default
configuration may not work (or may not work well) out of the box in
some environments.  However, you can easily tweak Google Test by
defining control macros on the compiler command line.  Generally,
these macros are named like `GTEST_XYZ` and you define them to either 1
or 0 to enable or disable a certain feature.

We list the most frequently used macros below.  For a complete list,
see file [include/gtest/internal/gtest-port.h](include/gtest/internal/gtest-port.h).

### Choosing a TR1 Tuple Library ###

Some Google Test features require the C++ Technical Report 1 (TR1)
tuple library, which is not yet available with all compilers.  The
good news is that Google Test implements a subset of TR1 tuple that's
enough for its own need, and will automatically use this when the
compiler doesn't provide TR1 tuple.

Usually you don't need to care about which tuple library Google Test
uses.  However, if your project already uses TR1 tuple, you need to
tell Google Test to use the same TR1 tuple library the rest of your
project uses, or the two tuple implementations will clash.  To do
that, add

    -DGTEST_USE_OWN_TR1_TUPLE=0

to the compiler flags while compiling Google Test and your tests.  If
you want to force Google Test to use its own tuple library, just add

    -DGTEST_USE_OWN_TR1_TUPLE=1

to the compiler flags instead.

If you don't want Google Test to use tuple at all, add

    -DGTEST_HAS_TR1_TUPLE=0

and all features using tuple will be disabled.

### Multi-threaded Tests ###

Google Test is thread-safe where the pthread library is available.
After `#include "gtest/gtest.h"`, you can check the `GTEST_IS_THREADSAFE`
macro to see whether this is the case (yes if the macro is `#defined` to
1, no if it's undefined.).

If Google Test doesn't correctly detect whether pthread is available
in your environment, you can force it with

    -DGTEST_HAS_PTHREAD=1

or

    -DGTEST_HAS_PTHREAD=0

When Google Test uses pthread, you may need to add flags to your
compiler and/or linker to select the pthread library, or you'll get
link errors.  If you use the CMake script or the deprecated Autotools
script, this is taken care of for you.  If you use your own build
script, you'll need to read your compiler and linker's manual to
figure out what flags to add.

### As a Shared Library (DLL) ###

Google Test is compact, so most users can build and link it as a
static library for the simplicity.  You can choose to use Google Test
as a shared library (known as a DLL on Windows) if you prefer.

To compile *gtest* as a shared library, add

    -DGTEST_CREATE_SHARED_LIBRARY=1

to the compiler flags.  You'll also need to tell the linker to produce
a shared library instead - consult your linker's manual for how to do
it.

To compile your *tests* that use the gtest shared library, add

    -DGTEST_LINKED_AS_SHARED_LIBRARY=1

to the compiler flags.

Note: while the above steps aren't technically necessary today when
using some compilers (e.g. GCC), they may become necessary in the
future, if we decide to improve the speed of loading the library (see
<http://gcc.gnu.org/wiki/Visibility> for details).  Therefore you are
recommended to always add the above flags when using Google Test as a
shared library.  Otherwise a future release of Google Test may break
your build script.

### Avoiding Macro Name Clashes ###

In C++, macros don't obey namespaces.  Therefore two libraries that
both define a macro of the same name will clash if you `#include` both
definitions.  In case a Google Test macro clashes with another
library, you can force Google Test to rename its macro to avoid the
conflict.

Specifically, if both Google Test and some other code define macro
FOO, you can add

    -DGTEST_DONT_DEFINE_FOO=1

to the compiler flags to tell Google Test to change the macro's name
from `FOO` to `GTEST_FOO`.  Currently `FOO` can be `FAIL`, `SUCCEED`,
or `TEST`.  For example, with `-DGTEST_DONT_DEFINE_TEST=1`, you'll
need to write

    GTEST_TEST(SomeTest, DoesThis) { ... }

instead of

    TEST(SomeTest, DoesThis) { ... }

in order to define a test.

## Developing Google Test ##

This section discusses how to make your own changes to Google Test.

### Testing Google Test Itself ###

To make sure your changes work as intended and don't break existing
functionality, you'll want to compile and run Google Test's own tests.
For that you can use CMake:

    mkdir mybuild
    cd mybuild
    cmake -Dgtest_build_tests=ON ${GTEST_DIR}

Make sure you have Python installed, as some of Google Test's tests
are written in Python.  If the cmake command complains about not being
able to find Python (`Could NOT find PythonInterp (missing:
PYTHON_EXECUTABLE)`), try telling it explicitly where your Python
executable can be found:

    cmake -DPYTHON_EXECUTABLE=path/to/python -Dgtest_build_tests=ON ${GTEST_DIR}

Next, you can build Google Test and all of its own tests.  On \*nix,
this is usually done by 'make'.  To run the tests, do

    make test

All tests should pass.

Normally you don't need to worry about regenerating the source files,
unless you need to modify them.  In that case, you should modify the
corresponding .pump files instead and run the pump.py Python script to
regenerate them.  You can find pump.py in the [scripts/](scripts/) directory.
Read the [Pump manual](docs/PumpManual.md) for how to use it.
HyperCanny Documentation
========================

