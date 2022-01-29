# Clingo: A grounder and solver for logic programs

Clingo is part of the [Potassco](https://potassco.org) project for *Answer Set
Programming* (ASP).  ASP offers a simple and powerful modeling language to
describe combinatorial problems as *logic programs*.  The *clingo* system then
takes such a logic program and computes *answer sets* representing solutions to
the given problem.  To get an idea, check our [Getting
Started](https://potassco.org/doc/start/) page and the [online
version](https://potassco.org/clingo/run/) of clingo.

Please consult the following resources for further information:

  - [**Downloading source and binary releases**](https://github.com/potassco/clingo/releases)
  - [**Installation and software requirements**](INSTALL.md)
  - [Changes between releases](CHANGES.md)
  - [Documentation](https://github.com/potassco/guide/releases)
  - [Potassco clingo page](https://potassco.org/clingo/)

Clingo is distributed under the [MIT License](LICENSE.md).
# clingo 5
- address `_` in negated literals
  - use `_` or maybe `*`
- add sort-constraint
  - `order(B,A) :- (A, B) = #sort{ X : p(X) }.`
  - `order(A,B) :- ((_,A), (_,B)) = #sort{ K,X : p(X), key(X,K) }.`
- profiling
- sorting via conditional literals became less efficient with the latest implementation in some cases
- projection is disabled in non-monotone constructs for now
  it could be enabled again if equivalences are used for affected atoms
- remove CSP support
- shifting of disjunctions

# constraints
- integrate constraint variables tighter into the gringo language
  - any term can contain csp variables
  - use [HT\_LC](http://www.cs.uni-potsdam.de/wv/pdfformat/cakaossc16a.pdf) as basis
- csp terms could be supported in aggregates too
- when the value of a csp term is required rewriting has to happen
  - `p(f($X+3))` becomes `p(f(Aux)), Aux = $X+3`
- supporting csp terms in functions symbols in build-ins is more tricky
  - `f($X+3) < g($Y+1)`
  - a translation as above is unsafe
  - the literal is a tautology

# memory management
- I think that the proper way to handle symbols is to use reference counting
  - the c api would have to put the burden of reference counting on the user
  - the APIs can provide reference counted symbols
- symbols should only be owned in key locations where the reference count is
  increased
  - for example a domain owns a symbol and can also take care of flyweighting
    it
- touching the counters should be avoided as much as possible during grounding
  - during the backtraking based instantiation a very limited set of
    intermediate symbols is created
  - these intermediate symbols can be freed upon backtracking (if they have not
    been added to a domain)
  - the temporary symbols may refer to owned symbols without increasing their
    reference count
  - this has the advantage that insertion into hash tables and touching of
    reference counts only happens when a symbol is commited into a domain
- implementing this should not be terribly difficult but affects *a lot* of
  code

# misc
- **enlarge test suites**
- incremental programs
  - atm indexes have to be cleared and recreated afterwards
    - it might be a good idea to optimize this and reuse indices later on
    - for now just clear them to not have them dangling around
- missing features in view of the ASP standard
  - queries
- assignment rewriting
  - enqueue: `expr(X,Z,Val):-expr(X,Y,Val_1)?,sing_term(Y,Z,Val_2)?,Val=(Val_1+Val_2),#X0=(Val_1+Val_2),#X0=Val.`
  - handle assignments in a more clever way...
- it would be nice to block grounding of rules if one index is empty
  (maybe even delaying the filling of indices if one index is empty)
- indices could be specialized to handle zero-ary predicates more efficiently
  - there could be one domain for all zero-ary predicates
- detect ground rules and implement more clever groundig
- domains
  - using a value as representation is wasteful
  - uses one unordered\_map to much
- predicate indices
  - using a valvec as key is wasteful
  - uses one unordered\_map too much
- on large instances both optimizations should safe a lot of memory
# Changes

## clingo 5.5.1
  * extend theory class to get version information
  * improve performance of `Model.symbol` (#296)
  * tidy up `clingo.hh` header regarding C++17 deprecations (#344)
  * fix error handling while solving in Python API (#334)
  * fix various outher bugs

## clingo 5.5.0
  * allow for using `not` as a theory operator (#193)
  * improve parsing of disjunctions (#171)
  * rename `max_size` to `size` in APIs and remove previous `size` method
    (breaks backward compatibility but makes things consistent)
  * refine cmake configuration variables (#283)
    (package maintainers might have to pay attention here)
  * extend clingo API (#236, #231, #228, #227, #187, #171, #174, #183)
  * add type annotations and add stub files for type checkers (#212)
  * reimplement Python API using CFFI (#253)
    (breaks backward compatibility regarding the AST and some details)
  * add `remove_watch` and `freeze_literal` to `propagate_init` (#285)
  * add support for Lua 5.4
  * add options to select specific Lua/Python versions
  * add single-shot solving option (#269)
    (see installation instructions)
  * rename `clingo.Tuple` to `clingo.Tuple_`
    (to avoid name clashes with `typing.Tuple`)
  * fix propagator initialization (#166)
  * fix cleanup function and handling of theory (#169)
  * fix shifting of disjunctions (in clasp) (#173)
  * fix handling of pools in externals (#179)
  * fix logger in Python API (#175)
  * fix memory bugs regarding strings (#196)
  * fix undo function (#191)
  * fix adding literals in propagate init (#192)
  * fix handling of undefined arithmetics (#218)
  * fix incremental grounding (#248)
  * fix/improve handling of classical negation (#268)
  * update to clasp version 3.3.6 fixing various issues
## clingo 5.4.0
  * add extension to implement custom heuristics
  * add const modifiers to C API
  * add flags to external and const statements to match API functions
  * fix python memleaks
  * make compatible with msvc 15
  * C ABI changes
    * extended propagators
  * C++ API changes
    * unify usage of ProgramBuilder and Backend
  * python API changes
    * `TruthValue.{_True,_False}` becomes `TruthValue.{True_,False_}`
    * `HeuristicType.{True,False}` becomes `TruthValue.{True_,False_}`
    * `async` and yield becomes `async_` and `yield_`
  * improve python API documentation
  * use cmakes `find_python` module
  * update to clasp version 3.3.5
## clingo 5.3.0
  * change C API to use numeric instead of symbolic literals
    * affects assumptions and assigning/releasing externals
      (breaks backward compatibility)
    * added overloads to C++, python and lua API to support
      both numeric and symbolic version
      (preserves backward compatibility for most code)
  * the python, C and C++ APIs now allow for customizing clingo by implementing
    a custom main function but reusing the rest of the application including
    the standard output
  * add API function to detect conflicting programs
  * add message logger to python and lua interface
  * add support for primes in the beginning of identifiers and variable names
  * add per solver registration of watches during propagator initialization
  * add a directive to selectivel suppress undefined atom warnings
  * add support for user defined statistics
  * add _to_c functions for python API to be able to call C functions from
    python
  * only create ground representations for requested program parts when
    grounding (#71)
  * improve program observer (#19)
  * support for binary, octal, and hexadecimal numbers (#65)
  * the backend has to be opened/closed now
  * release python's GIL while grounding (#82)
  * TruthValue.{True,False} becomes TruthValue.{\_True,\_False} in python API
  * improve API and it's documentation
## clingo 5.2.3
  * update to clasp version 3.3.4
  * really fix --pre option
  * fix link to potassco guide (#74)
  * fix location printing (#78)
  * fix linking problems (#79)
  * fix modulo zero (#100)
  * fix web builds (#103)
  * fix addding clauses after a model has been found (#104)
  * fix python program observer bindings (#105)
  * expose exponentiation operator in APIs
  * improve python docstrings (#101, #102)
  * add option to build python and lua modules against an existing libclingo
## clingo 5.2.2
  * update to clasp version 3.3.3
  * use GNUInstallDirs in cmake files to simplify packaging
  * fix --pre option
  * fix swapped clingo\_assignment\_size and clingo\_assignment\_max\_size
  * fix docstrings
  * fix incremental mode
  * fix sup and inf in python/lua bindings
  * fix reified format term tuples
  * fix wrong use of python API (causing trouble with python 3.6)
  * fix compilation problems on 32bit linux (missing libatomic)
## clingo 5.2.1
  * update to clasp version 3.3.2
  * fix handling of istop in incmode programs
  * fix handling of undefined ** operations
  * fix preprocessing of disjunctions with undefined operations
    (regression in clingo-5)
  * fix segfault during preprocessing
    (regression in clingo-5)
## clingo 5.2.0
  * switch to MIT license
  * improve compatibility with abstract gringo
  * switch build system from scons to cmake
  * improve windows compatibility
  * make tests and examples python 3 compatible
  * bison and re2c are no longer required to build source releases
  * update to clasp 3.3.0
  * the CLINGOPATH environment variable can be set
    to control from where to include files in logic programs
  * propagators can add variables while solving now
  * refactor interfaces (breaking backward compatibility)
    * there is just one solve function now
    * in the C API do not pass structs by value to functions
      because FFIs of some languages do not support this
  * fix cleanup function
  * numerous other bugfixes not listed here
## clingo 5.1.1
  * fix thread id to start with one in propagator.undo in lua
  * fix version macro in clingo.h
  * fix added missing methods to get thread id to model in lua/python
  * fix child\_key property in python ast
## clingo 5.1.0
  * update to clasp 3.2.1
  * add interface to add variables during propagation
  * add interface to inspect ground rules (C/C++ only)
  * add experimental interface to access clasp facade (C/C++ only)
  * fixed smodels output (--output=smodels)
## clingo 5.0.0
  * cleanup of python and lua API (breaks backwards compatibility)
  * added new aspif output format replacing the old smodels format
  * added input language support for clasp features
    * #edge directives to add acyclicity constraints
    * #project directives for enumeration of projected models
    * #heuristic directives to steer clasp's search
  * added theory atoms to write aggregate like constructs
  * added stable C API documented with doxygen
  * added experimental C++ API based on C API
  * added theory propagator interface to clingo APIs
  * added support for compilation with Visual Studio 2015
  * improved data structures to reduce memory consumption on typical input
  * updated to clasp version 3.2.0 + patches
## gringo/clingo 4.5.4
  * fixed bug when creating multiple Control objects
    (affects lua only)
  * fixed bug when trying to configure more solvers than in portfolio
    (affects python only)
  * fixed #disjoint constraints
  * improved build scripts
  * added option to keep facts in normal rules
## gringo/clingo 4.5.3
  * fixed regression w.r.t gringo 4.4 in translation of conditional literals
  * fixed projection in incremental programs
  * fixed bug with (double) negative literals in minimize constraints
## gringo/clingo 4.5.2
  * fixed memory leak in python API when enumerating models
  * updated to clasp version 3.1.3
## gringo/clingo 4.5.1
  * ground term parser returns None/nil for undefined terms now
  * added warning if a global variable occurs in a tuple of an aggregate element
  * added auto detection of libraries
  * changed option --update-domains into API function Control:cleanup\_domains
  * fixed domain cleanup when used with minimize constraints
  * fixed grounding of recursive disjunctions (regression in 4.5.0)
  * fixed Control.stats in lua bindings
  * fixed a bug in clingo that would print 0-ary classically negated atoms wrongly
## gringo/clingo 4.5.0
  * fixed grounding of recursive aggregates
  * fixed usage of lua\_next
  * fixed bug when applying constant definitions
  * updated underlying clasp to version 3.1.1
  * added support for negation in front of relation literals
  * added option --update-domains to cleanup gringo's domains
    using clasp's top-level assignment when solving incrementally
  * added domain inspection to scripting interface
  * added term parser to scripting interface
  * added support for python 3 (experimental)
  * added support for one elementary tuples
  * added support for unary - operator in front of functions and symbols
  * added support for recursive nonmonotone aggregate via translation
  * added program reify to reify logic programs
  * added option to rewrite minimize constaints for use with reify
  * changed inbuilt iclingo mode
    (breaks backwards compatibility)
  * changed handling of pools, intervals, and undefined operations according to AG
    (breaks backwards compatibility)
  * changed handling of ==, it is treated like = now
  * changed SolveFuture.interrupt to SolveFuture.cancel
    (breaks backwards compatibility)
## gringo/clingo 4.4.0
  * updated underlying clasp to version 3.1.0
    * this version brings numerous fixes regarding incremental solving
  * scripting API changes
    * ground takes a list of programs to ground now and immediately starts
      grounding (breaks backwards compatibility)
    * asolve has been renamed to solveAsync
      (breaks backwards compatibility)
    * the solver configuration is better integrated now
      (breaks backwards compatibility)
    * solver statistics are a property now
      (breaks backwards compatibility)
    * added a method to add clauses during solving
    * added load method to load files
    * added solveIter method to iterate over methods without using a callback
    * added optional assumptions to solve/solveAsync/solveIter method
    * enableEnumAssumption became a property
  * added library that can be imported in python
  * rules with fact heads where not simplified in all cases
  * fixed grounding of recursive aggregates
  * fixed translation of aggregates with multiple guards
## gringo/clingo 4.3.0
  * fixed bug with incremental parameters in minimize constraints
  * fixed handling of empty tuples
  * fixed translation of conditional literals
  * fixed translation of factual body aggregates
  * fixed bug not properly recognizing aggregates as non-monotone
  * fixed bug not properly grounding recursive head aggregates
  * fixed bug with recursive negated aggregates
  * fixed bug with head aggregates with multiple elements
  * improved handling of conditional literals
  * added method to get optimization values of model in scripting language
  * clingo uses clasp 3.0 now
## gringo/clingo 4.2.1
  * fixed bug in simplification of aggregates
  * fixed bug with raw strings in macros
  * fixed compilation issues with older glibc versions
  * fixed output for enumeration of cautious consequences
  * fixed bugs in clasp library
    * fixed race in parallel model enumeration
    * fixed incremental optimization
    * fixed cleanup up of learnt constraints during incremental solving
  * workaround for libstdc++'s bad choice for hash<uint64_t> on 32bit arches
## gringo/clingo 4.2
  * added clingo 
    * supports very flexible scripting support
    * can cover iclingo and oclingo functionality now
  * added stack traces to lua error messages
  * added support for incremental optimization
  * improved python error messages
  * renamed gringo.Function to gringo.Fun
  * removed luabind dependency
  * removed boost-python dependency
  * consistently use not instead of #not as keyword for negation
  * fixed translation of conditions in head aggregates
  * fixed replacement of constants
  * fixed grounding of recursive head aggregates
  * fixed translation of head aggregates
  * fixed show statements for CSP variables (condition was ignored)
  * fixed plain text output of body aggregates
  * added a ton of new bugs
## gringo 4.1
  * added scripting languages python and lua
  * added -c option to define constants
  * added constraints over integer variables
    * linear constraints
    * disjoint constraints
    * show statements for constraint variables
    * (experimental and subject to change)
  * improved translation of disjunctions
  * fixed include directives
  * fixed preprocessing of definitions
  * fixed lparse translation of optimization constructs

# Table of Contents

- [Requirements](#requirements)
  - [Development Dependencies](#development-dependencies)
  - [Optional Dependencies](#optional-dependencies)
- [Build, Install, and Test](#build-install-and-test)
  - [Build Options](#build-options)
    - [Generic Options](#generic-options)
    - [Python Support](#python-support)
    - [Lua Support](#lua-support)
  - [Compilation to JavaScript](#compilation-to-javascript)
- [Troubleshooting](#troubleshooting)
  - [Notes for Windows Users](#notes-for-windows-users)

# Requirements

This document is about installing clingo from source. We also provide
precompiled packages for various package managers:
<https://potassco.org/clingo/#packages>.

- a C++14 conforming compiler
  - *at least* [GCC](https://gcc.gnu.org/) version 4.9
  - [Clang](http://clang.llvm.org/) version 3.1 (using either libstdc++
    provided by gcc 4.9 or libc++)
  - *at least* MSVC 14.0 ([Visual Studio](https://www.visualstudio.com/) 2015
    Update 3)
  - other compilers might work
- the [cmake](https://www.cmake.org/) build system
  - at least version 3.18 is recommended
  - at least version 3.1 is *required*

## Development Dependencies

The following dependencies are only required when compiling a development
branch. Releases already include the necessary generated files.

- the [bison](https://www.gnu.org/software/bison/) parser generator
  - *at least* version 2.5
  - version 3.0 produces harmless warnings
    (to stay backwards-compatible)
- the [re2c](https://re2c.org/) lexer generator
  - *at least* version 0.13 is required

## Optional Dependencies

- the [Python](https://www.python.org/) script language and the
  [CFFI](cffi.readthedocs.io/) package
  - *at least* Python version 3.6
  - at least CFFI 1.14 is required (earlier versions have not been tested)
- the [Lua](https://www.lua.org/) script language
  - *at least* Lua 5.1 is required

# Build, Install, and Test

When cloning the git repository, do not forget to update the submodules (with
source releases, you can skip this step):

    git submodule update --init --recursive

To build gringo, clingo, and reify in their default configurations in release
mode, run:

    cmake -H<SOURCE_DIR> -B<BUILD_DIR> -DCMAKE_BUILD_TYPE=Release
    cmake --build <BUILD_DIR>

The resulting binaries and shared libraries will be in `<BUILD_DIR>/bin` and
are ready to use.

To install all binaries and development files under cmake's install
prefix (see the [build options](#build-options)), run:

    cmake --build <BUILD_DIR> --target install

To run the tests, enable option `CLINGO_BUILD_TESTS` (see [build
options](#build-options)) and run:

    cmake --build <BUILD_DIR> --target test

## Build Options

Cmake's `-L` option can be used to get an overview over the variables that can
be set for building gringo/clingo. To get gringo/clingo specific options, run

    cmake -H<SOURCE_DIR> -B<BUILD_DIR> -DCMAKE_BUILD_TYPE=Release -LH

or, to also print important cmake specific configuration variables

    cmake -H<SOURCE_DIR> -B<BUILD_DIR> -DCMAKE_BUILD_TYPE=Release -LAH

Options and variables can be passed to
cmake on the command line using `-D<VARIABLE>=<VALUE>` or by editing
`<BUILD_DIR>/CMakeCache.txt` after running cmake.

The build scripts by default try to detect optional dependencies, like Python
and Lua scripting support.

Clingo uses [libpotassco](https://github.com/potassco/libpotassco) and
[clasp](https://github.com/potassco/potassco).  Both components have their own
sets of configuration variables:
- [building libpotassco](https://github.com/potassco/libpotassco#installation)
- [building clasp](https://github.com/potassco/clasp#building--installing)

In the following, the most important options to control the build are listed.

### Generic Options

- Variable `CMAKE_BUILD_TYPE` should be set to `Release`.
- Variable `CMAKE_INSTALL_PREFIX` controls where to install clingo.
- Option `CLINGO_MANAGE_RPATH` controls how to find libraries on platforms
  where this is supported, like Linux, macOS, or BSD but not Windows. This
  option should be enabled if clingo is installed in a non-default location,
  like the users home directory; otherwise it has no effect.
  (Default: `ON`)
- Option `CLINGO_BUILD_APPS` controls whether to build the applications gringo,
  clingo, and reify.
  (Default: `ON`)
- Option `CLINGO_BUILD_EXAMPLES` controls whether to build the clingo API
  examples.
  (Default: `OFF`)
- Option `CLINGO_BUILD_TESTS` controls whether to build the clingo tests and
  enable the test target running unit as well as acceptance tests.
  (Default: `OFF`)

### Python Support

With the default configuration, Python support will be auto-detected if the
Python development packages are installed.

- Varibale `CLINGO_BUILD_WITH_PYTHON` can be set to `ON` to enable Python
  support, `OFF` to disable Python support, `auto` to enable Python support if
  available, or `pip` for advanced configuration to build a Python module
  exporting clingo symbols.
  (Default: `auto`)
- Variable `CLINGO_PYTHON_VERSION` can be used to select a specific Python
  version. For example, passing `-DCLINGO_PYTHON_VERSION:LIST="3.8;EXACT"` to
  cmake requires Python version 3.8 to be available and not just a compatbile
  Python version. Starting with cmake 3.15, variable `Python_ROOT` can be used
  to specify where to search for a Python installation (the documentation of
  cmake's `FindPython` module has further information).
  (Default: `3.6`)
- Variable `PYCLINGO_INSTALL` controls where to install the Python module. It
  can be set to `user` to install in the user prefix, `system` to install in
  the system prefix, or `prefix` to install into the installation prefix.
  (Default: `prefix`)
- Variable `PYCLINGO_SUFFIX` can be used to customize which suffix to use for
  the Python module.
  (Default: automatically detected)
- Variable `PYCLINGO_INSTALL_DIR` can be used to customize where to install the
  Python module.
  (Default: automatically detected)

With cmake versions before 3.15, it can happen that the found Python
interpreter does not match the found Python libraries if the development
headers for the interpreter are not installed. Make sure to install them before
running cmake (or remove or adjust the `CMakeCache.txt` file).

### Lua Support

With the default configuration, Lua support will be auto-detected if the Lua
development packages are installed.

- Varibale `CLINGO_BUILD_WITH_LUA` can be set to `ON` to enable Lua support,
  `OFF` to disable Lua support, `auto` to enable Lua support if available.
  (Default: `auto`)
- Variable `CLINGO_LUA_VERSION` can be used to select a specific Lua
  version. For example, passing `-DCLINGO_LUA_VERSION:LIST="5.3;EXACT"` to
  cmake requires Lua version 5.3 to be available and not just a compatbile
  Lua version.
  (Default: `5.0`)
- Variable `LUACLINGO_SUFFIX` can be used to customize which suffix to use for
  the Lua module.
  (Default: automatically detected)
- Variable `LUACLINGO_INSTALL_DIR` can be used to customize where to install
  the Lua module.
  (Default: automatically detected)

## Compilation to JavaScript

Clingo can be compiled to JavaScript with Empscripten. The following notes
assume that [Emscripten](https://kripken.github.io/emscripten-site/) has been
installed. Only the web target and a subset of clingo's configuration are
supported when compiling to JavaScript:

    emcmake cmake -H<SOURCE_DIR> -B<BUILD_DIR> \
        -DCLINGO_BUILD_WEB=On \
        -DCLINGO_BUILD_WITH_PYTHON=Off \
        -DCLINGO_BUILD_WITH_LUA=Off \
        -DCLINGO_BUILD_SHARED=Off \
        -DCLASP_BUILD_WITH_THREADS=Off \
        -DCMAKE_VERBOSE_MAKEFILE=On \
        -DCMAKE_BUILD_TYPE=release \
        -DCMAKE_CXX_FLAGS="-std=c++11 -Wall -s DISABLE_EXCEPTION_CATCHING=0" \
        -DCMAKE_CXX_FLAGS_RELEASE="-Os -DNDEBUG" \
        -DCMAKE_EXE_LINKER_FLAGS="" \
        -DCMAKE_EXE_LINKER_FLAGS_RELEASE=""
    cmake --build <BUILD_DIR> --target web

Note that is is possible to enable Lua support. Therefore Lua has to be
compiled with emscripten, too. See [Lua Support](#lua-support) for information
about pointing clingo to a custom Lua installation.

For examples how to use the resulting JavaScript code, check out one of the
following:
- [webclingo example by Lucas Bourneuf](https://github.com/Aluriak/webclingo-example), or
- [the source of our website](https://github.com/potassco/potassco.github.io)

# Troubleshooting

After installing the required packages clingo should compile on most \*nixes.
If a dependency is missing or a software version too old, then there are
typically community repositories that provide the necessary packages. To list a
few:
- the [ToolChain](https://wiki.ubuntu.com/ToolChain) repository for Ubuntu
  versions before 18.04.
- the [Developer
  Toolset](https://wiki.centos.org/SpecialInterestGroup/SCLo/CollectionsList)
  for CentOS
- the [Cygwin](http://cygwin.org) project under Windows (re2c must be compiled
  by hand)
- both [Homebrew](https://brew.sh/) and [MacPorts](https://www.macports.org/)
  provide all the software necessary to compile clingo

And, well, you can compile a recent gcc version yourself. Even on ancient Linux
systems. ;)

## Notes for Windows Users

clingo can be compiled using the
[Mingw-w64](https://mingw-w64.sourceforge.net/) compiler, the Cygwin project,
or Visual Studio 2015 Update 3. For development, The [re2c](https://re2c.org/)
and [winflexbison](https://github.com/lexxmark/winflexbison) (providing
[bison](https://www.gnu.org/software/bison/) for Windows) packages can be
installed using the [chocolatey](https://chocolatey.org/) package manager, too.
# Clingo: A grounder and solver for logic programs

Clingo is part of the [Potassco](https://potassco.org) project for *Answer Set
Programming* (ASP).  ASP offers a simple and powerful modeling language to
describe combinatorial problems as *logic programs*.  The *clingo* system then
takes such a logic program and computes *answer sets* representing solutions to
the given problem.  To get an idea, check our [Getting
Started](https://potassco.org/doc/start/) page and the [online
version](https://potassco.org/clingo/run/) of clingo.

Clingo is distributed under the [MIT License](LICENSE.md).
# Sorting Terms

Sorting a set of terms in ASP is quite difficult to implement efficiently.
Nicely readable programs like

    next(X,Z) :- p(X), #false : X<Y, p(Y), Y<Z ; p(Z), X<Z.

or

    next(X,Y) :- p(X), Y = #min { Z : p(Z), Z > X }, Y != #sup.

do not scale well enough for large input (qubic or quadratic for the programs
above).

In practice it might perform best to already sort the terms in the instance. Or
the scripting API can be used to obtain better performing code.

Note that in the example there are no next predicates but an enumerate
predicate, which is more flexible in practice. For example, to encode the
above:

    next(X,Y) :- enumerate(N,X), enumerate(N+1,Y).

## Example calls

    gringo --text sort-lua.lp encoding.lp
    gringo --text sort-py.lp encoding.lp
# Solving the Towers of Hanoi Problem

This example shows how to incrementally solve the towers of hanoi problem using
clingo's inbuild incremental solving mode.

Note that, with a fixed bound, the problem can also be grounded by gringo and
then solved by clasp.

## Example Calls

    clingo tohI.lp tohE.lp
    gringo -c imax=16 tohE.lp tohI.lp | clasp
## Meta Encodings

Generic meta encodings used by some of the examples in the reify folder.
# Austere logic programs

Computing stable models using a meta-encoding where default negation only
occurs in the constraints, as in austere logic programs [1].


## Example Call

Finding all stable models of `example.lp`:

    $ clingo --output=reify example.lp | clingo -Wno-atom-undefined - encoding.lp 0


## References

[1] Jorge Fandinno, Seemran Mishra, Javier Romero, Torsten Schaub: Answer Set Programming Made Easy. Submitted for publication
# Many diverse stable models

Computing many diverse stable models of logic programs


## Example Calls

Finding 3 different (or 1-diverse) stable models of `example.lp`:

    $ clingo --output=reify example.lp | clingo -Wno-atom-undefined - encoding.lp -c m=3 -c option=1 -c k=1

Finding 3 6-diverse stable models of `example.lp`:

    $ clingo --output=reify example.lp | clingo -Wno-atom-undefined - encoding.lp -c m=3 -c option=1 -c k=6

Finding 3 most diverse stable models of `example.lp`:

    $ clingo --output=reify example.lp | clingo -Wno-atom-undefined - encoding.lp -c m=3 -c option=2 --quiet=1,2,2

Example Calls
=============

The encodings here have been developed within the [metasp] project allowing for
handling complex optimization criteria, e.g., inclusion-based minimization or
Pareto efficiency.

Non-reified + Cardinality Minimization
--------------------------------------

    $ clingo --opt-mode=optN -q1 0 example1.lp

Reified + Cardinality Minimization
----------------------------------

    $ clingo --rewrite-minimize --output=reify --reify-sccs example1.lp |\
      clingo -Wno-atom-undefined - encoding.lp \
      <(echo "optimize(0,1,card).") 0

Reified + Subset Minimization
-----------------------------

    $ clingo --rewrite-minimize --output=reify --reify-sccs example1.lp |\
      clingo -Wno-atom-undefined - encoding.lp \
      <(echo "optimize(0,1,incl).") 0

Reified + Subset Minimization + Query to Solve the Conformant Planning Problem
------------------------------------------------------------------------------

    $ clingo example2.lp --output=reify |\
      clingo -Wno-atom-undefined - encoding.lp --project 0


Improving Performance
=====================

The above calls use clingo's reification backend. Better performance can be
obtained by first preprocessing the program and then reifying a potentially
much smaller program:

    $ clingo --pre --rewrite-minimize ... |\
      reify --sccs |\
      clingo -Wno-atom-undefined - encoding.lp ...


[metasp]: https://potassco.org/labs/metasp/
# Computing supported models

Using meta-programming for computing supported models of logic programs.


## Example Call

Finding all supported models of `example1.lp`:

    $ clingo --output=reify example1.lp | clingo -Wno-atom-undefined - encoding.lp 0


## Warning

The call returns supported models of the ground logic program that results
after grounding the input program (`example1.lp`). You can run the following
call to see that ground logic program:

    $ clingo --text example1.lp

The grounding process applies simplifications that may eliminate some atoms of
the input program. To avoid this, in `example1.lp` we have added an external
declaration. You can comment it and run the call again to see what happens.
# Guess and Check Programming

An implementation of Guess and Check Programming [1] in clingo using
meta-programming (described in [2]).

A Guess and Check program is a pair of logic programs `<G,C>`. A set of atoms
`X` is a stable model of `<G,C>` if `X` is a stable model of `G`, and `C`
together with `H` is unsatisfiable, where `H` contains the facts `guess(x).`
for all atoms of the form `guess(x)` in `X`.

The implementation translates a Guess and Check program into a disjunctive
logic program, that is solved by clingo.


## Usage

The following call computes solutions for the given `[guess_programs]` and
`[check_programs]`:

    $ clingo --output=reify domain.lp [guess_programs]             | \
    grep "output(guess (.*))"                                      | \
    clingo --output=reify --reify-sccs - guess.lp [check_programs] | \
    clingo -Wno-atom-undefined - glue.lp [guess_programs]

Note that the `[guess_programs]` have to be given twice. The first time the
programs are passed is to compute the domain to ground the `[check_programs]`.

There is also a small script to simplify the above call:

    $ run.sh [guess_programs] -- [check_programs] -- [options]

See the readmes in examples subfolder for further informations.


## Notes

Predicate `guess/1` should not appear in any head of the `check_programs`, and
predicate `output/1` should not appear in the `check_programs`.


## References

[1] Thomas Eiter, Axel Polleres: Towards automated integration of guess and check programs in answer set programming: a meta-interpreter and applications. TPLP 6(1-2): 23-60 (2006).
[2] Roland Kaminski, Javier Romero, Torsten Schaub, Philipp Wanko: How to build your own ASP system?! Submitted for publication
# Conformant planning

Conformant planning can be represented in `QBF` style as:

* There exists an initial situation `I`, and a sequence of actions `S` that
  achieve the goal starting from `I`, such that for all initial situations
  `I'`, `S` achieves the goal starting from `I'`.

This is equivalent to:

* There exists an initial situation `I`, and a sequence of actions `S` that
  achieve the goal starting from `I`, such that there is no `I'` such that `S`
  does not achieve the goal starting from `I'`.

In our approach the guess program represents the first exists part:

* find an initial situation `I`, and a sequence of actions `S` that achieve the
  goal starting from `I`.

The check program represents the second part:

* find an initial situation `I'`, such that `S` (given by the guess program)
  does not achieve the goal starting from `I'`.

    $ ../../run.sh base.lp instance.lp guess.lp -- base.lp instance.lp check.lp
    clingo version 5.5.0
    Reading from - ...
    Solving...
    Answer: 1
    occurs(cpa_go_down(cpa_e0,cpa_f1,cpa_f0),1) occurs(cpa_step_in(cpa_e0,cpa_f0,cpa_p0),2)
    occurs(cpa_go_up(cpa_e0,cpa_f0,cpa_f1),3) occurs(cpa_step_out(cpa_e0,cpa_f1,cpa_p0),4)
    occurs(cpa_collect(cpa_c1,cpa_f1,cpa_p0),5) occurs(cpa_collect(cpa_c0,cpa_f1,cpa_p0),6)
    occurs(cpa_move_right(cpa_f1,cpa_p0,cpa_p1),7) occurs(cpa_collect(cpa_c1,cpa_f1,cpa_p1),8)
    occurs(cpa_collect(cpa_c0,cpa_f1,cpa_p1),9)
    SATISFIABLE

    Models       : 1+
    Calls        : 1
    Time         : 264.631s (Solving: 264.46s 1st Model: 264.46s Unsat: 0.00s)
    CPU Time     : 264.527s
# Preferences 

Call to compute all superset maximal stable models of `guess.lp`:

    $ ../../run.sh base.lp guess.lp -- base.lp check_superset.lp -- 0
    clingo version 5.5.0
    Reading from - ...
    Solving...
    Answer: 1
    a(1) a(2)
    SATISFIABLE

    Models       : 1
    Calls        : 1
    Time         : 0.009s (Solving: 0.00s 1st Model: 0.00s Unsat: 0.00s)
    CPU Time     : 0.012s

Call to compute all subset minimal stable models of `guess.lp`:

    $ ../../run.sh base.lp guess.lp -- base.lp check_subset.lp -- 0
    clingo version 5.5.0
    Reading from - ...
    Solving...
    Answer: 1
    a(2)
    Answer: 2
    a(1)
    SATISFIABLE

    Models       : 2
    Calls        : 1
    Time         : 0.009s (Solving: 0.00s 1st Model: 0.00s Unsat: 0.00s)
    CPU Time     : 0.012s
# Simple Tic-Tac-Toe

In this example there is a 3x3 Tic-Tac-Toe square. The guessing player has to
place her 3 tokens in a winning position such that afterwards the checking
player cannot place her tokens in a winning position:

    $ ../../run.sh base.lp guess.lp -- base.lp check.lp -- 0
    clingo version 5.5.0
    Reading from - ...
    Solving...
    Answer: 1
    m(1,1) m(2,2) m(3,3)
    Answer: 2
    m(1,3) m(2,2) m(3,1)
    SATISFIABLE

    Models       : 2
    Calls        : 1
    Time         : 0.013s (Solving: 0.00s 1st Model: 0.00s Unsat: 0.00s)
    CPU Time     : 0.012s

# Solving 2QBF

Given the following 2QBF:

    EXISTS {x(1),x(2)} FORALL {y(1),y(2)} (-x(1) OR -y(1)) AND (-x(2) v -y(2))

the next call returns the unique solution, where both `x(1)` and `x(2)` are false:

    $ ../../run.sh base.lp guess.lp -- check.lp -- 0
    Reading from - ...
    Solving...
    Answer: 1

    SATISFIABLE

    Models       : 1
    Calls        : 1
    Time         : 0.008s (Solving: 0.00s 1st Model: 0.00s Unsat: 0.00s)
    CPU Time     : 0.008s
# Simple Example

This is the call for a simple example:

    $ ../../run.sh base.lp guess.lp -- check.lp -- 0
    clingo version 5.5.0
    Reading from - ...
    Solving...
    Answer: 1
    a(2)
    SATISFIABLE

    Models       : 1
    Calls        : 1
    Time         : 0.007s (Solving: 0.00s 1st Model: 0.00s Unsat: 0.00s)
    CPU Time     : 0.008s
# Here-and-There models

Computing Here-and-There models of logic programs


## Example Calls

Finding all Here-and-There models of `example1.lp`:

    $ clingo --output=reify example1.lp | clingo -Wno-atom-undefined - encoding.lp 0 -c option=1

Finding all Here-and-There models minimizing the Here world of `example1.lp`:

    $ clingo --output=reify example1.lp | clingo -Wno-atom-undefined - encoding.lp 0 -c option=2

Finding all Equilibrium models of `example1.lp`:

    $ clingo --output=reify example1.lp | clingo -Wno-atom-undefined - encoding.lp 0 -c option=3
# Simple Meta Encoding

This example simply reproduces the answer sets of the reified program.


## Example Call

Finding all stable models of `example.lp`:

    $ clingo --output=reify example.lp | clingo -Wno-atom-undefined - ../common/meta.lp 0
# Computing Classical Models

Using meta-programming for computing classical models of logic programs.


## Example Call

Finding all classical models of `example1.lp`:

    $ clingo --output=reify example1.lp | clingo -Wno-atom-undefined - encoding.lp 0

## Warning


The call returns classical models of the ground logic program that results
after grounding the input program (`example1.lp`). You can run the following
call to see that ground logic program:

    $ clingo --text example1.lp

The grounding process applies simplifications that may eliminate some atoms of
the input program. To avoid this, in `example1.lp` we have added an external
declaration. You can comment it and run the call again to see what happens.

# Compte the well-founded model of a program

This examples computes the well-founded model of a normal logic program. It
prints true and unknown atoms before solving. Note that it does not print false
atoms because there are in general infinitely many of them.

Disclaimer: Writing an efficient algorithm to compute the well-founded model is
actually quite tricky. This example is not very well tested - there might be
bugs.

## Dependencies

The example needs the clingox and networkx python packages in addition to the
clingo package. The packages can be installed using:

    python3 -m pip install --user --upgrade --extra-index-url https://test.pypi.org/simple/ clingo-cffi clingox
    python3 -m pip install --user --upgrade networkx

# Example Calls

    ➜ python well-founded.py example.lp
    level version 1.0
    Reading from example.lp
    Facts:
    r s
    Unknown:
    u v x y
    Solving...
    UNSATISFIABLE

    Models       : 0
    Calls        : 1
    Time         : 0.001s (Solving: 0.00s 1st Model: 0.00s Unsat: 0.00s)
    CPU Time     : 0.001s
# Example that Shows how to Export an Answer Set to Excel

This example exports each predicate to different sheet in the file
`excel.xlsx`. To run the program, you need a clingo build with Python support
and the Python modules: `pandas`, `xlsxwriter`, `openpyxl`

## Example calls

    clingo excel-py.lp example.lp
# Modeling Transition Systems

This example is similar to clingo's incremental mode but additionally takes
care of attaching a time parameter to atoms.

# Example Calls

    python tmode.py example.lp
# Extending Models for Printing

This example shows how to add symbols to a model before printing.

# Example Calls

    clingo extend-model-py.lp 0
    clingo extend-model-lua.lp 0
# 15 Puzzle

This is an example how to solve the [15 puzzle] problem.  It supports two
solving modes: a non-consecutive mode where each tile is moved individually and
a consecutive mode where multiple tiles are moved at once.  In the latter mode,
much harder instances can be solved but the resulting solutions are not
necessarily optimal w.r.t. the number of single tile moves.

## Example Calls

    clingo encoding.lp instance1.lp -c consecutive=1
    clingo encoding.lp instance1.lp -c consecutive=0
    clingo encoding.lp instance2.lp -c consecutive=1 -t8,split --config=crafty

[15 puzzle]: https://en.wikipedia.org/wiki/15_puzzle
Exmaple to show different usages of the clingo module.

Examples
========

Embedding python code to implement a function evaluated during grounding:

    $ clingo embedded.lp example.lp 
    clingo version 5.5.0
    Reading from embedded.lp ...
    Solving...
    Answer: 1
    num(3) num(6) div(3,1) div(3,3) div(6,1) div(6,2) div(6,3) div(6,6)
    SATISFIABLE

    Models       : 1
    Calls        : 1
    Time         : 0.028s (Solving: 0.00s 1st Model: 0.00s Unsat: 0.00s)
    CPU Time     : 0.028s

Do something similar as above but use clingo's python module:

    $ python module.py
    num(3) num(6) div(3,1) div(3,3) div(6,1) div(6,2) div(6,3) div(6,6)

Writing a custom application:

    $ python app.py example.lp
    example version 1.0
    Reading from example.lp
    Solving...
    Answer: 1
    num(3) num(6) div(3,1) div(3,3) div(6,1) div(6,2) div(6,3) div(6,6)
    SATISFIABLE

    Models       : 1+
    Calls        : 1
    Time         : 0.003s (Solving: 0.00s 1st Model: 0.00s Unsat: 0.00s)
    CPU Time     : 0.001s
Example show casing multi-shot solving. It implements branch-and-bound-based
optimization and incremental solving.

Examples
========

The branch-and-bound example:

    $ python opt.py tohE.lp tohI.lp tohB.lp -c n=20
    opt-example version 1.0
    Reading from tohE.lp ...
    Solving...
    Answer: 1
    move(3,c,2)  move(4,b,1)  move(4,c,3)  move(2,b,4)  move(4,b,5)  move(4,a,6)  move(3,b,7)
    move(3,c,8)  move(4,b,9)  move(4,a,10) move(3,b,11) move(4,b,12) move(1,c,13) move(4,c,14)
    move(3,a,15) move(4,a,16) move(2,c,17) move(4,b,18) move(3,c,19) move(4,c,20)
    Found new bound: 20
    Solving...
    ...
    Solving...
    Answer: 1
    move(3,c,2)  move(4,b,1)  move(4,c,3)  move(2,b,4)  move(4,a,5)  move(3,b,6)  move(4,b,7)
    move(1,c,8)  move(4,c,9)  move(3,a,10) move(4,a,11) move(2,c,12) move(4,b,13) move(3,c,14)
    move(4,c,15)
    Found new bound: 15
    Solving...
    Optimum found
    UNSATISFIABLE
    
    Models       : 5
    Calls        : 6
    Time         : 0.015s (Solving: 0.01s 1st Model: 0.01s Unsat: 0.00s)
    CPU Time     : 0.015s

The incremental solving example:

    $ python inc.py tohE.lp tohI.lp
    inc-example version 1.0
    Reading from tohE.lp ...
    Solving...
    ...
    Solving...
    Answer: 1
    move(4,b,1)  move(3,c,2)  move(4,c,3)  move(2,b,4)  move(4,a,5)  move(3,b,6)  move(4,b,7)
    move(1,c,8)  move(4,c,9)  move(3,a,10) move(4,a,11) move(2,c,12) move(4,b,13) move(3,c,14)
    move(4,c,15)
    SATISFIABLE
    
    Models       : 1+
    Calls        : 16
    Time         : 0.020s (Solving: 0.00s 1st Model: 0.00s Unsat: 0.00s)
    CPU Time     : 0.020s
# Incremental Grounding and Solving

With clingo-4 there is no iclingo version dedicated to incremental solving
anymore.  This example shows how the same functionality can be implemented
using clingo's API.

For convenience, this functionality has been build into clingo and can be used
without requiring additional files:

    #include <incmode>.

## Example calls

    clingo incmode-int.lp example.lp
    clingo incmode-lua.lp example.lp
    clingo incmode-py.lp example.lp
# Incremantally Solving the n-Queens Problem

In this example, we calculate solutions for the n-Queens problem for different
board sizes, which can be given as a list of intervals on the command line.

## Example Calls

    clingo incqueens.lp incqueens-py.lp -c calls="list((1,1),(3,5),(8,9))"
    clingo incqueens.lp incqueens-lua.lp -c calls="list((1,1),(3,5),(8,9))"

This example shows how to implement a guess and check based propagator.

Examples
========

    $ python app.py example.lp
    guess-and-check version 1.0
    Reading from example.lp
    Solving...
    Answer: 1
    a(2)
    SATISFIABLE
    
    Models       : 1+
    Calls        : 1
    Time         : 0.001s (Solving: 0.00s 1st Model: 0.00s Unsat: 0.00s)
    CPU Time     : 0.001s
# A Ricochet Robots Solver

Alex Randolph's board game Ricochet Robots offers a rich and versatile
benchmark for ASP. As it stands, it represents a simple multi-agent planning
problem in which each agent, i.e., robot, has limited sensing capacities (that
is, only bumps are detected).

This little script allows you to first select a target and then move the
robots. If you run out of ideas, you can click on solve! to get a little
help from ASP. To start, call

    python visualize.py

Simple example to show how to use program parts and externals.

Examples
========

    $ python app.py
    Exmaple 1:
      loading files:
      - chemistry.lp
      grounding:
      - base
      solutions:
      - a(1) a(2)
    
    Exmaple 2:
      loading files:
      - chemistry.lp
      grounding:
      - acid(42)
      solutions:
      - b(42)
    
    Exmaple 3:
      loading files:
      - chemistry.lp
      - external.lp
      grounding:
      - base
      - acid(42)
      assigning externals:
      - d(1,42)=True
      solutions:
      - a(1) a(2) b(42) c(1,42) c(2,42) d(1,42) e(1,42)
# Peg Solitaire

## The Problem

Given the board below, find a solution consisting of successive moves such that
only one peg (o) remains in the center of the board.  A move is an orthogonal
jump of length two above an adjacent peg to a free position on the board.  The
jumped peg has to be removed after the move.


        -------
        |o|o|o|
        -------
        |o|o|o|
    ---------------
    |o|o|o|o|o|o|o|
    ---------------
    |o|o|o| |o|o|o|
    ---------------
    |o|o|o|o|o|o|o|
    ---------------
        |o|o|o|
        -------
        |o|o|o|
        -------

## Setup

The visualizer requires urwid.  To install under Debian it, run as root:

    aptitude install python-urwid

## Finding a Solution

The problem can be solved using the incremental encoding in solitaire.lp.  To
find a solution, call (takes about 10s on a fast machine):

    clingo solitaire.lp instance.lp

To visualize the solution, call:

    ./visualize.py solitaire.lp instance.lp

## Notes

The encoding allows for parallel moves, i.e., moves that do not influence each
other may occur at the same time step.  This is essential to solve the problem.
When visualizing the solution the moves are sequentialized again.
# User Heuristic

This is an example how to implement a simple domain-specific heuristic for the
graph coloring problem.

## Example Calls

    clingo encoding-py.lp instance.lp
    clingo encoding-lua.lp instance.lp
# Solving the Towers of Hanoi Problem

This example solves Towers of Hanoi problems.  There are two variants.  First,
there is an incremental encoding.  Second, there is an optimizing version that
tries to find a shortest plan within a given planning horizon.

# Example Calls

Using the incremental version:

    $ clingo tohE.lp tohI.lp

Using the bounded version:

    $ clingo opt.lp tohB.lp tohI.lp -c n=30
This example implements a scaled down version of clingo-dl supporting
constraints of the following form:

    &diff{u-v} <= d :- [rule body].

Examples
========

Enumerating all solutions of a flowshop problem.

    $ python app.py fsE.lp fsI.lp 0
    clingo-dl version 1.0
    Reading from fsE.lp ...
    Solving...
    Answer: 1
    dl((a,1),1) dl(bound,16) dl((a,2),7) dl((b,1),0) dl((b,2),1) dl((c,1),4) dl((c,2),11) permutation(b,a) permutation(a,c)
    Answer: 2
    dl((a,1),6) dl(bound,20) dl((a,2),16) dl((b,1),5) dl((b,2),10) dl((c,1),0) dl((c,2),5) permutation(b,a) permutation(c,b)
    Answer: 3
    dl((a,1),0) dl(bound,19) dl((a,2),3) dl((b,1),8) dl((b,2),13) dl((c,1),3) dl((c,2),8) permutation(c,b) permutation(a,c)
    Answer: 4
    dl((a,1),6) dl(bound,16) dl((a,2),12) dl((b,1),0) dl((b,2),1) dl((c,1),1) dl((c,2),7) permutation(b,c) permutation(c,a)
    Answer: 5
    dl((a,1),0) dl(bound,18) dl((a,2),3) dl((b,1),3) dl((b,2),7) dl((c,1),4) dl((c,2),13) permutation(b,c) permutation(a,b)
    Answer: 6
    dl((a,1),5) dl(bound,20) dl((a,2),10) dl((b,1),8) dl((b,2),14) dl((c,1),0) dl((c,2),5) permutation(c,a) permutation(a,b)
    SATISFIABLE
    
    Models       : 6
    Calls        : 1
    Time         : 0.006s (Solving: 0.00s 1st Model: 0.00s Unsat: 0.00s)
    CPU Time     : 0.006s

Finding the optimal solution of a flowshop problem.

    $ python app.py --minimize-variable bound fsE.lp fsI.lp
    clingo-dl version 1.0
    Reading from examples/fsE.lp ...
    Solving...
    Answer: 1
    dl((a,1),1) dl(bound,16) dl((a,2),7) dl((b,1),0) dl((b,2),1) dl((c,1),4) dl((c,2),11) permutation(b,a) permutation(a,c)
    Found new bound: 16
    Solving...
    Optimum found
    UNSATISFIABLE
    
    Models       : 1
    Calls        : 2
    Time         : 0.006s (Solving: 0.00s 1st Model: 0.00s Unsat: 0.00s)
    CPU Time     : 0.006s
# Using Propagators and Theory Language to Ease Debugging

To ease development of encodings, `&cannot` atoms can be used in rule heads to
get useful debugging messages.  The idea is to use them as rule heads where
otherwise integrity constraints would be used.  If an encoding produces unsat
results, this can be used to relax integrity constraints.  The script takes
care of minimizing the number of messages printed.  To get an idea, take a look
at `example.lp` and run one of the following calls:

    $ clingo cannot-py.lp example.lp
    $ clingo cannot-lua.lp example.lp
---
layout: home
---
