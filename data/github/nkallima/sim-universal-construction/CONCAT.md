Changelog
---------

v3.1.1
------
- Bug fix for the pool objects on the recycling functionality.

v3.1.0
------
- A new high-performant, optimized implementation for the flat-combining synchronization technique. This new implementation provided by the Synch framework is written from the scratch.
- A few API fixes.

v3.0.1
------
- README.md updated.
- Adding the JOSS paper to the repository.

v3.0.0
------
- Header-files cleanup in order to provide a better and consistent API.
- Adding the `Synch` prefix to all structs provided by the runtime/primitives to the end-user.
- Adding the `synch` prefix to all functions provided by the runtime/primitives to the end-user.
- Extensive documentation for the API.
- Create a `doxygen` configuration file for auto-generating man-pages/documentation.
- API documentation is now provided at GitHub pages: https://nkallima.github.io/sim-universal-construction/index.html.
- A code coverage report is provided through codecov.io for the validation script.
- Support for installing the framework.
- Improvements on the build environment.
- Expansion of the CONTRIBUTING.md.
- Providing a short discussion (see PERFORMANCE.md) of the expected performance for the various objects provided by the framework.
- The `bench.sh` and `validate.sh` scripts use similar arguments as the provided benchmarks.

v2.4.0
------
- Adding the LCRQ Queue implementation (Adam Morrison and Yehuda Afek, PPoPP 2013, http://mcg.cs.tau.ac.il/projects/lcrq).
- Memory reclamation for SimStack is supported.
- Adding support for 128 bit Compare&Swap.
- Documentation improvements; a memory reclamation section was added in `README.md`.

v2.3.1
------
- Improvements on `validate.sh` for better error reporting.
- Enhancements on `README.md`.
- Improved logo.

v2.3.0
------
- Introducing `validate.sh`, a smoke/validation script that verifies many of the provided data structures.
- A few improvements on some benchmarks in order to be conformed with `validate.sh`.
- Introducing the logo of the repository.
- Enhancements on `README.md`.

v2.2.1
------
- Executing `Pause` instructions on spinning-loops for the X86 architecture. Spin-locks and some lock-free algorithms may have performance benefits in SMT architectures.  
- Documentation improvements.

v2.2.0
------
- A new implementation of the pool.c functionality. Most stack and queue implementations support memory reclamation, the only exception are the simstack, simqueue, lfstack and msqueue implementations.
- Code clean-up for the Sim family of algorithms.
- Improvements on the HSynch family of algorithms.
- Enhanced NUMA support for the NUMA-aware data-structures.
- Homogenized code for stack and queue implementations.
- Better output for bench.sh script.
- Numerous performance optimizations and bug-fixes, especially in machines with weak memory models.
- Documentation improvements.

v2.1.1
------
- Improving the output messages of build-system.
- Adding support for clang in Makefile (i.e.,use `make clang` for building the sources with the clang compile).
- A few bug-fixes and minor improvements.
- Dropped support for Intel icc compiler and added support for Intel icx compiler. Tested with Intel(R) oneAPI DPC++ Compiler 2021.1.2.

v2.1.0
------
- Significant enhancements on bench.sh.
- New folder/file library structure.
- ftime is replaced by clock_gettime.

v2.0.1
------
- Removing dead-code from the unsupported sparc/solaris architecture.
- PAPI and NUMA libraries automatically used depending user's settings in config.h.
- Pid fix for the Sync family of algorithms.

v2.0
----
- Initial support for ARM-V8 and RISC-V machine architectures.
- Significant performance optimizations for the Synch family of algorithms on modern NUMA machines.
- The produced binaries can run with any number of threads/fibers without re-compilation.
- All the produced benchmark binaries can accept arguments, such as the number of threads, the number of benchmark's iterations, etc.
- Enhancements on how the runtime handles the NUMA characteristics of modern multiprocessors.
- Better thread affinity policy in NUMA machines.
- A new script called run_all.sh is provided for automatically testing all the produced binaries.
- A .clang-format file is provided in order to maintain the styling-consistency of the source code.
- Numerous performance optimizations and bug-fixes.

v1.9.1
------
- Testing has been performed in RISC-V and aarch64 machine architectures.
- Bug fixes and performance improvements.

v1.9
----
- Support for ARM-V8 and RISC-V machine architectures.

v1.8
----
- The last version that supports solaris/sparc architecture. From now on sparc/solaris architecture is unsupported.
[![check-build](https://github.com/nkallima/sim-universal-construction/actions/workflows/check-build.yml/badge.svg)](https://github.com/nkallima/sim-universal-construction/actions/workflows/check-build.yml) [![validate](https://github.com/nkallima/sim-universal-construction/actions/workflows/validate.yml/badge.svg)](https://github.com/nkallima/sim-universal-construction/actions/workflows/validate.yml) [![codecov](https://codecov.io/gh/nkallima/sim-universal-construction/branch/codecov/graph/badge.svg?token=1V8A6BOABM)](https://codecov.io/gh/nkallima/sim-universal-construction) [![status](https://joss.theoj.org/papers/07bf35ba1bd72c38cc8076fee6864409/status.svg)](https://joss.theoj.org/papers/07bf35ba1bd72c38cc8076fee6864409)

<p align="center">
    <img src="resources/logo_synch.png" alt="The Synch Framework" width="80%">
</p>

# Summary

This is an open-source framework for concurrent data-structures and benchmarks. The provided framework contains a substantial set of concurrent data-structures such as `queues`, `stacks`, `combining-objects`,
`hash-tables`, `locks`, etc. This framework also provides a user-friendly runtime for developing and benchmarking concurrent data-structures. Among other features, this runtime provides functionality for creating threads easily (both Posix and user-level threads), tools for measuring performance, etc. The provided concurrent data-structures and the runtime are highly optimized for contemporary NUMA multiprocessors such as AMD Epyc and Intel Xeon.

The current version of this code is optimized for x86_64 machine architecture, but the code is also successfully tested in other machine architectures, such as ARM-V8 and RISC-V. Some of the benchmarks perform much better in architectures that natively support Fetch&Add instructions (e.g., x86_64, etc.).


# Collection

The Synch framework provides a large set of highly efficient concurrent data-structures, such as combining-objects, concurrent queues and stacks, concurrent hash-tables and locks. The cornerstone of the Synch framework are the combining objects. A Combining object is a concurrent object/data-structure that is able to simulate any other concurrent object, e.g. stacks, queues, atomic counters, barriers, etc. The Synch framework provides the PSim wait-free combining object [2,10], the blocking combining objects CC-Synch, DSM-Synch and H-Synch [1], and the blocking combining object based on the technique presented in [4]. Moreover, the Synch framework provides the Osci blocking, combining technique [3] that achieves good performance using user-level threads. Since v3.1.0, the Synch framework offers a new high performant implementation of flat-combining synchronization technique [14]. This novel version is implemented from the scratch and is not just an optimized version of the original code provided in [15].

In terms of concurrent queues, the Synch framework provides the SimQueue [2,10] wait-free queue implementation that is based on the PSim combining object, the CC-Queue, DSM-Queue and H-Queue [1] blocking queue implementations based on the CC-Synch, DSM-Synch and H-Synch combining objects. A blocking queue implementation based on the CLH locks [5,6] and the lock-free implementation presented in [7] are also provided.
Since v2.4.0, the Synch framework provides the LCRQ [11,12] queue implementation. In terms of concurrent stacks, the Synch framework provides the SimStack [2,10] wait-free stack implementation that is based on the PSim combining object, the CC-Stack, DSM-Stack and H-Stack [1] blocking stack implementations based on the CC-Synch, DSM-Synch and H-Synch combining objects. Moreover, the lock-free stack implementation of [8] and the blocking implementation based on the CLH locks [5,6] are provided. The Synch framework also provides concurrent queue and stacks implementations (i.e. OsciQueue and OsciStack implementations) that achieve very high performance using user-level threads [3]. Since v3.1.0, the Synch framework provides stack and queue implementations (i.e. FC-Stack and FC-Queue) based on the  implementation of flat-combining provided by the Synch framework.

Furthermore, the Synch framework provides a few scalable lock implementations, i.e. the MCS queue-lock presented in [9] and the CLH queue-lock presented in [5,6]. Finally, the Synch framework provides two example-implementations of concurrent hash-tables. More specifically, it provides a simple implementation based on CLH queue-locks [5,6] and an implementation based on the DSM-Synch [1] combining technique.

The following table presents a summary of the concurrent data-structures offered by the Synch framework.
| Concurrent  Object    |                Provided Implementations                           |
| --------------------- | ----------------------------------------------------------------- |
| Combining Objects     | CC-Synch, DSM-Synch and H-Synch [1]                               |
|                       | PSim [2,10]                                                       |
|                       | Osci [3]                                                          |
|                       | Oyama [4]                                                         |
|                       | FC: a new implementation of flat-combining [14]                   |
| Concurrent Queues     | CC-Queue, DSM-Queue and H-Queue [1]                               |
|                       | SimQueue [2,10]                                                   |
|                       | OsciQueue [3]                                                     |
|                       | CLH-Queue [5,6]                                                   |
|                       | MS-Queue [7]                                                      |
|                       | LCRQ [11,12]                                                      |
|                       | FC-Queue [14]                                                     |
| Concurrent Stacks     | CC-Stack, DSM-Stack and H-Stack [1]                               |
|                       | SimStack [2,10]                                                   |
|                       | OsciStack [3]                                                     |
|                       | CLH-Stack [5,6]                                                   |
|                       | LF-Stack [8]                                                      |
|                       | FC-Stack [14]                                                     |
| Locks                 | CLH [5,6]                                                         |
|                       | MCS [9]                                                           |
| Hash Tables           | CLH-Hash [5,6]                                                    |
|                       | A hash-table based on DSM-Synch [1]                               |


# Requirements

- A modern 64-bit machine. Currently, 32-bit architectures are not supported. The current version of this code is optimized for the x86_64 machine architecture, but the code is also successfully tested in other machine architectures, such as ARM-V8 and RISC-V. Some of the benchmarks perform much better in architectures that natively support Fetch&Add instructions (e.g., x86_64, etc.).
- A recent Linux distribution. The Synch environment may also build/run in some other Unix-like systems, (i.e. BSD, etc.). In this case the result is not guaranteed, since the environment is not tested in systems other than Linux.
- As a compiler, gcc of version 4.8 or greater is recommended, but you may also try to use icx or clang.
- Building requires the following development packages:
    - `libatomic`
    - `libnuma`
    - `libpapi` in case that the `SYNCH_TRACK_CPU_COUNTERS` flag is enabled in `libconcurrent/config.h`.
- For building the documentation (i.e. man-pages), `doxygen` is required.


# Configuring, compiling and installing the framework

In the `libconcurrent/config.h` file, the user can configure some basic options for the framework, such as:
- Enable/disable debug mode.
- Support for Numa machines.
- Enable performance statistics, etc.

The provided default configuration should work well in many cases. However, the default configuration may not provide the best performance. For getting the best performance, modifying the `libconcurrent/config.h` may be needed (see more on Performance/Optimizations Section).

In case that you want to compile the library that provides all the implemented concurrent algorithms just execute `make` in the root directory of the source files. This step is necessary in case that you want to run benchmarks. However, some extra make options are provided in case the user wants to compile the framework with other than system's default compiler, clean the binary files, etc. The following table provides the list with all the available make options.

|     Command             |                                Description                                                                    |
| ----------------------- | ------------------------------------------------------------------------------------------------------------- |
|  `make`                 |  Auto-detects the current architecture and compiles the source-code for it (this should work for most users). |
|  `make CC=cc`           |  Compiles the source-code for the current architecture using the `cc` compiler.                               |
|  `make clang`           |  Compiles the source-code using the clang compiler.                                                           |
|  `make icx`             |  Compiles the source-code using the Intel icx compiler.                                                       |
|  `make unknown`         |  Compiles the source-code for architectures other than X86_64, e.g. RISC-V, ARM, etc.                         |
|  `make clean`           |  Cleaning-up all the binary files.                                                                            |
|  `make docs`            |  Creating the documentation (i.e. man-pages).                                                                 |
|  `make install`         |  Installing the framework on the default location (i.e. `/opt/Synch/`).                                       |
|  `make install DIR=dir` |  Installing the framework on the `dir/Synch/` location.                                                       |
|  `make uninstall`       |  Uninstalling the framework.                                                                                  |

For building the documentation (i.e. man-pages), the user should execute `make docs`. Notice that for building the documentation the system should be equipped with `doxygen` documentation tool.

For installing the framework, the user should execute `make install`. In this case, the framework will be installed in the default location which is `/opt/Synch/`. Notice that in this case, the user should have write access on the `/opt` directory or sudo access. The `make install DIR=dir` command installs the framework in the `dir/Synch` path, while the `make uninstall` uninstalls the framework. For accessing the man pages, the user should manually setup the `MANPATH` environmental variable appropriately (e.g. `export MANPATH=$MANPATH:/opt/Synch/docs/man`).


# Running Benchmarks

For running benchmarks use the `bench.sh` script file that is provided in the main directory of this source tree.

Example usage: `./bench.sh FILE.run OPTION1 VALUE1 OPTION2 VALUE2 ...`

Each benchmark reports the time that needs to be completed, the average throughput of operations performed and some performance statistics if `DEBUG` option is enabled during framework build. The `bench.sh` script measures the strong scaling of the benchmark that is executed.

The following options are available:

|     Option              |                       Description                                                     |
| ----------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------ |
|  `-t`, `--max_threads`  |  set the maximum number number of POSIX threads to be used in the last set of iterations of the benchmark, default is the number of system cores |
|  `-s`, `--step`         |  set the step (extra number of threads to be used) in each set of iterations of the benchmark, default is number of processors/8 or 1            |
|  `-f`, `--fibers`       |  set the number of user-level threads per posix thread                                                                                           |
|  `-r`, `--repeat`       |  set the total number of operations executed by the benchmark, default is 1000000                                                                |
|  `-i`, `--iterations`   |  set the number of times that the benchmark should be executed, default is 10                                                                    |
|  `-w`, `--workload`     |  set the amount of workload (i.e. dummy loop iterations among two consecutive operations of the benchmarked object), default is 64               |
|  `-l`, `--list`         |  displays the list of the available benchmarks                                                                                                   |
|  `-n`, `--numa_nodes`   |  set the number of numa nodes (which may differ with the actual hw numa nodes) that hierarchical algorithms should take account                  |
|  `-b`, `--backoff`, `--backoff_high` |  set an upper backoff bound for lock-free and Sim-based algorithms                                                                  |
|  `-bl`, `--backoff_low` |  set a lower backoff bound (only for msqueuebench, lfstackbench and lfuobjectbench benchmarks)                                                                  |
|  `-h`, `--help`         |  displays this help and exits                                                                                                                    |

The framework provides the `validate.sh` validation/smoke script. The `validate.sh` script compiles the sources in `DEBUG` mode and runs a big set of benchmarks with various numbers of threads. After running each of the benchmarks, the script evaluates the `DEBUG` output and in case of success it prints `PASS`. In case of a failure, the script simply prints `FAIL`. In order to see all the available options of the validation/smoke script, execute `validate.sh -h`. Given that the `validate.sh` validation/smoke script depends on binaries that are compiled in `DEBUG` mode, it is not installed while using `make install`. The following image shows the execution and the default behavior of `validate.sh`.

![](resources/validate_example.gif)

The framework provides another simple fast smoke test: `./run_all.sh`. This will quickly run all available benchmarks with default options and store the results in the `results.txt` file.


# Performance/Optimizations

Getting the best performance from the provided benchmarks is not always an easy task. For getting the best performance, some modifications in Makefiles may be needed (compiler flags, etc.). Important parameters for the benchmarks and/or library are placed in the `libconcurrent/config.h` file. A useful guide to consider in order to get better performance in a modern multiprocessor follows.

- In case that the target machine is a NUMA machine make sure `SYNCH_NUMA_SUPPORT` is enabled in `libconcurrent/config.h`. Usually, when this option is enabled, it gives much better performance in NUMA machines. However, in some older machines this option may induce performance overheads.
- Whenever the `SYNCH_NUMA_SUPPORT` option is enabled, the runtime will detect the system's number of NUMA nodes and will setup the environment appropriately. However, significant performance benefits have been observed by manually setting-up the number of NUMA nodes manually (see the `--numa_nodes` option). For example, the performance of the H-Synch family algorithms on an AMD EPYC machine consisting of 2x EPYC 7501 processors (i.e., 128 hardware threads) is much better by setting `--numa_nodes` equal to `2`. Notice that the runtime successfully reports that the available NUMA nodes are `8`, but this value is not optimal for H-Synch in this configuration. An experimental analysis for different values of `--numa_nodes` may be needed.
- Check the performance impact of the `SYNCH_COMPACT_ALLOCATION` option in `libconcurrent/config.h`. In modern AMD multiprocessors (i.e., equipped with EPYC processors) this option gives tremendous performance boost. In contrast to AMD processors, this option introduces serious performance overheads in Intel Xeon processors. Thus, a careful experimental analysis is needed in order to show the possible benefits of this option.
- Check the cache line size (`CACHE_LINE_SIZE` and `S_CACHE_LINE` options in includes/system.h). These options greatly affect the performance in all modern processors. Most Intel machines behave better with `CACHE_LINE_SIZE` equal or greater than `128`, while most modern AMD machine achieve better performance with a value equal to `64`. Notice that `CACHE_LINE_SIZE` and `S_CACHE_LINE` depend on the `SYNCH_COMPACT_ALLOCATION` option (see includes/system.h).
- Use backoff if it is available. Many of the provided algorithms could use backoff in order to provide better performance (e.g., sim, LF-Stack, MS-Queue, SimQueue, SimStack, etc.). In this case, it is of crucial importance to use `-b` (and in some cases `-bl` arguments) in order to get the best performance. 
- Ensure that you are using a recent gcc-compatible compiler, e.g. a `gcc` compiler of version `7.0` or greater is highly recommended.
- Check the performance impact of the different available compiler optimizations. In most cases, gcc's `-Ofast` option gives the best performance. In addition, some algorithms (i.e., sim, osci, simstack, oscistack, simqueue and osciqueue) benefit by enabling the `-mavx` option (in case that AVX instructions are supported by the hardware).
- Check if system oversubscription with user-level fibers enhances the performance. Many algorithms (i.e., the Sim and Osci families of algorithms) show tremendous performance boost by using oversubscription with user-level threads [3]. In this case, use the `--fibers` option.

## Expected performance

The expected performance of the Synch framework is discussed in the [PERFORMANCE.md](PERFORMANCE.md) file.

# Memory reclamation (stacks and queues)

The Synch framework provides a pool mechanism (see `includes/pool.h`) that efficiently allocates and de-allocates memory for the provided concurrent stack and queue implementations. The allocation mechanism of this pool implementation is low-overhead.  All the provided stack and queue implementations use the functionality of this pool mechanism. In order to support memory reclamation in a safe manner, a concurrent object should guarantee that each memory object that is going to de-allocated should be accessed only by the thread that is going to free it. Generally, de-allocating and thus reclaiming memory is easy in many blocking objects, since there is a lock that protects the de-allocated memory object. Currently, the Synch framework supports memory reclamation for the following concurrent stack and queue implementations:
- Concurrent Queues:
    - CC-Queue, DSM-Queue and H-Queue [1]
    - OsciQueue [3]
    - CLH-Queue [5,6]
- Concurrent Stacks:
    - CC-Stack, DSM-Stack and H-Stack [1]
    - OsciStack [3]
    - CLH-Stack [5,6]
    - SimStack [2,10] (since v2.4.0)

Note that de-allocating and thus recycling memory in lock-free and wait-free objects is not an easy task. Since v2.4.0, SimStack supports memory reclamation using the functionality of `pool.h` and a technique that is similar to that presented by Blelloch and Weiin in [13]. Notice that the MS-Queue [7], LCRQ [11,12] queue implementations and the LF-Stack [8] stack implementation support memory reclamation through hazard-pointers. However, the current version of the Synch framework does not provide any implementation of hazard-pointers. In case that a user wants to use memory reclamation in these objects, a custom hazard-pointers implementation should be integrated in the environment.

By default, memory-reclamation is enabled. In case that there is need to disable memory reclamation, the `SYNCH_POOL_NODE_RECYCLING_DISABLE` option should be enabled in `config.h`.

The following table shows the memory reclamation characteristics of the provided stack and queues implementations.

| Concurrent  Object    |        Provided Implementations           | Memory Reclamation                        |
| --------------------- | ----------------------------------------- | ----------------------------------------- |
| Concurrent Queues     | CC-Queue, DSM-Queue and H-Queue [1]       | Supported                                 |
|                       | SimQueue [2,10]                           | Not supported                             |
|                       | OsciQueue [3]                             | Supported                                 |
|                       | CLH-Queue [5,6]                           | Supported                                 |
|                       | MS-Queue [7]                              | Hazard Pointers (not provided by Synch)   |
|                       | LCRQ [11,12]                              | Hazard Pointers (not provided by Synch)   |
|                       | FC-Queue [14]                             | Supported                                 |
| Concurrent Stacks     | CC-Stack, DSM-Stack and H-Stack [1]       | Supported                                 |
|                       | SimStack [2,10]                           | Supported (since v2.4.0)                  |
|                       | OsciStack [3]                             | Supported                                 |
|                       | CLH-Stack [5,6]                           | Supported                                 |
|                       | LF-Stack [8]                              | Hazard Pointers (not provided by Synch)   |
|                       | FC-Stack [14]                             | Supported                                 |


## Memory reclamation limitations

In the current design of the reclamation mechanism, each thread uses a single private pool for reclaiming memory. In a producer-consumer scenario where a set of threads performs only enqueue operations (or push operations in case of stacks) and all other threads perform dequeue operations (or pop operations in case of stacks), insufficient memory reclamation is performed since each memory pool is only accessible by the thread that owns it. We aim to improve this in future versions of the Synch framework.


# API documentation

A complete API documentation is provided in [https://nkallima.github.io/sim-universal-construction/index.html](https://nkallima.github.io/sim-universal-construction/index.html).

# Code example for a simple benchmark
We now describe a very simple example-benchmark that uses the Application Programming Interface (API) of the provided runtime. This simple benchmark measures the performance of Fetch&Add instructions in multi-core machines. The purpose of this simple benchmark is to measure the performance of Fetch&Add implementations (hardware or software).

```c
#include <stdio.h>
#include <stdint.h>

#include <primitives.h>
#include <threadtools.h>
#include <barrier.h>

#define N_THREADS 10
#define RUNS      1000000

volatile int64_t object CACHE_ALIGN;
int64_t d1 CACHE_ALIGN, d2;
SynchBarrier bar CACHE_ALIGN;

inline static void *Execute(void *Arg) {
    long i, id = (long)Arg;

    synchBarrierWait(&bar);
    if (id == 0) d1 = synchGetTimeMillis();

    for (i = 0; i < RUNS; i++)
        synchFAA64(&object, 1);

    synchBarrierWait(&bar);
    if (id == 0) d2 = synchGetTimeMillis();

    return NULL;
}

int main(int argc, char *argv[]) {
    object = 1;

    synchBarrierSet(&bar, N_THREADS);
    synchStartThreadsN(N_THREADS, Execute, SYNCH_DONT_USE_UTHREADS);
    synchJoinThreadsN(N_THREADS - 1);

    printf("time: %ld (ms)\tthroughput: %.2f (millions ops/sec)\n", 
           (d2 - d1), RUNS * N_THREADS / (1000.0 * (d2 - d1)));

    return 0;
}
```

This example-benchmark creates `N_THREADS`, where each of them executes `RUNS` Fetch&Add operations in a shared 64-bit integer. At the end of the benchmark the throughput (i.e. Fetch&Add operations per second) is calculated. By seting varous values for `N_THREADS`, this benchmark is able to measure strong scaling.

The `synchStartThreadsN` function (provided by the API defined in `threadtools.h`) in main, creates `N_THREADS` threads and each of the executes the `Execute` function declared in the same file. The `SYNCH_DONT_USE_UTHREADS_` argument imposes `synchStartThreadsN` to create only Posix threads; in case that the user sets the corresponding fibers argument to `M` > 0, then `synchStartThreadsN` will create `N_THREADS` Posix threads and each of them will create `M` user-level (i.e. fiber) threads. The `synchJoinThreadsN` function (also provided by `threadtools. h`) waits until all Posix and fiber (if any) threads finish the execution of the `Execute` function. The Fetch&Add instruction on 64-bit integers is performed by the `synchFAA64` function provided by the API of `primitives.h`.

The threads executing the `Execute` function use the `SynchBarrier` re-entrant barrier object for simultaneously starting to perform Fetch&Add instructions on the shared variable `object`. This barrier is also re-used before the end of the `Execute` function in order to allow thread with `id = 0` to measure the amount of time that the benchmark needed for completion. The `synchBarrierSet` function in `main` initializes the `SynchBarrier` object. The `synchBarrierSet` takes as an argument a pointer to the barrier object and the number of threads `N_THREADS` that are going to use it. Both `synchBarrierSet` and `synchBarrierWait` are provided by the API of `barrier.h`

At the end of the benchmark, `main` calculates and prints the average throughput of Fetch&Add operations per second achieved by the benchmark.


# If you want to cite us

```latex
@article{Kallimanis2021,
  doi = {10.21105/joss.03143},
  url = {https://doi.org/10.21105/joss.03143},
  year = {2021},
  publisher = {The Open Journal},
  volume = {6},
  number = {64},
  pages = {3143},
  author = {Nikolaos D. Kallimanis},
  title = {Synch: A framework for concurrent data-structures and benchmarks},
  journal = {Journal of Open Source Software}
}
```

# Releases

An extensive list of the recent releases and their features is provided at [https://github.com/nkallima/sim-universal-construction/releases](https://github.com/nkallima/sim-universal-construction/releases).

# License

The Synch framework is provided under the [LGPL-2.1 License](https://github.com/nkallima/sim-universal-construction/blob/main/LICENSE).

# Code of conduct

[Code of conduct](https://github.com/nkallima/sim-universal-construction/blob/main/.github/CODE_OF_CONDUCT.md).

# References

[1]. Panagiota Fatourou, and Nikolaos D. Kallimanis. "Revisiting the combining synchronization technique". ACM SIGPLAN Notices. Vol. 47. No. 8. ACM, PPoPP 2012.

[2]. Panagiota Fatourou, and Nikolaos D. Kallimanis. "A highly-efficient wait-free universal construction". Proceedings of the twenty-third annual ACM symposium on Parallelism in algorithms and architectures (SPAA), 2011.

[3]. Panagiota Fatourou, and Nikolaos D. Kallimanis. "Lock Oscillation: Boosting the Performance of Concurrent Data Structures" Proceedings of the 21st International Conference on Principles of Distributed Systems (Opodis), 2017.

[4]. Yoshihiro Oyama, Kenjiro Taura, and Akinori Yonezawa. "Executing parallel programs with synchronization bottlenecks efficiently". Proceedings of the International Workshop on Parallel and Distributed Computing for Symbolic and Irregular Applications. Vol. 16. 1999.

[5]. Travis S. Craig. "Building FIFO and priority-queueing spin locks from atomic swap". Technical Report TR 93-02-02, Department of Computer Science, University of Washington, February 1993.

[6]. Peter Magnusson, Anders Landin, and Erik Hagersten. "Queue locks on cache coherent multiprocessors". Parallel Processing Symposium, 1994. Proceedings., Eighth International. IEEE, 1994.
    
[7]. Maged M. Michael, and Michael L. Scott. "Simple, fast, and practical non-blocking and blocking concurrent queue algorithms". Proceedings of the fifteenth annual ACM symposium on Principles of distributed computing. ACM, 1996.
    
[8]. R. Kent Treiber. "Systems programming: Coping with parallelism". International Business Machines Incorporated, Thomas J. Watson Research Center, 1986.

[9]. John M. Mellor-Crummey, and Michael L. Scott. "Algorithms for scalable synchronization on shared-memory multiprocessors". ACM Transactions on Computer Systems (TOCS) 9.1 (1991): 21-65.

[10]. Panagiota Fatourou, and Nikolaos D. Kallimanis. "Highly-efficient wait-free synchronization". Theory of Computing Systems 55.3 (2014): 475-520.

[11]. Adam Morrison, and Yehuda Afek. "Fast concurrent queues for x86 processors". Proceedings of the 18th ACM SIGPLAN symposium on Principles and practice of parallel programming. 2013.

[12]. Adam Morrison, and Yehuda Afek. Source code for LCRQ. http://mcg.cs.tau.ac.il/projects/lcrq.

[13]. Guy E. Blelloch, and Yuanhao Wei. "Brief Announcement: Concurrent Fixed-Size Allocation and Free in Constant Time." 34th International Symposium on Distributed Computing (DISC 2020). Schloss Dagstuhl-Leibniz-Zentrum für Informatik, 2020.

[14]. Danny Hendler, Itai Incze, Nir Shavit, and Moran Tzafrir. Flat combining and the synchronization-parallelism tradeoff. In Proceedings of the twenty-second annual ACM symposium on Parallelism in algorithms and architectures (SPAA 2010), pp. 355-364.

[15]. Danny Hendler, Itai Incze, Nir Shavit, and Moran Tzafrir. Source code for flat-combing. https://github.com/mit-carbon/Flat-Combining


# Contact

For any further information, please do not hesitate to
send an email at nkallima (at) ics.forth.gr. Feedback is always valuable.
# Contribution

## Source code structure
The Synch framework consists of 3 main parts, i.e. the Runtime/Primitives, the library of concurrent data-structures and the benchmarks  (see the figure below) . The Runtime/Primitives part provides some basic functionality for creating and managing threads, functionality for basic atomic primitives (e.g. Compare&Swap, Fetch&Add, fences, simple synchronization barriers, etc.), mechanisms for basic memory allocation/management (e.g. memory pools, etc.), functionality for measuring time, reporting CPU counters, etc. Furthermore, the Runtime/Primitives provides a simple and lightweight library of user level-threads that can be used in order to evaluate the provided data-structures and algorithms. 

The Concurrent library utilizes the building blocks of the Runtime/Primitives layer in order to provide all the concurrent data-structures (e.g. combining objects, queues, stacks, etc.). For almost every concurrent data-structure or synchronization mechanism, Synch provides at least one benchmark for evaluating its performance.

<p align="center">
    <img src="resources/code_structure.png" alt="The code-structure of the Synch framework" width="50%">
</p>

## First, discuss your issue
Before contributing to this repository, you are encouraged to first discuss the issue, change or contribution that you wish to make. In order to do this, please open an issue or send an email. We are open to discuss new ideas, improvements or any other contributions. Please note that we have a [Code of conduct](https://github.com/nkallima/sim-universal-construction/blob/main/.github/CODE_OF_CONDUCT.md) that should be considered in all your interactions with the project.

## Ready to make a change
- Fork the repository.
- Create a new branch for your patch.
- Make your updates/changes.
- Open a pull request.
- Submit your pull request targeting the current development branch (i.e. the branch with the highest version number).
- Get it reviewed.

## Coding conventions
Before committing a new patch, follow the steps below for checking the coding style of your patch.
- Ensure that you have installed `clang-format` package (version >= 9).
- After cloning the repository you need to run `git config --local core.hooksPath .git_hooks` to enable recommended code style checks (placed in `.clang_format` file) during git commit. If the tool discovers inconsistencies, it will create a patch file. Please follow the instructions to apply the patch before opening a pull request.

In general, please conform your coding styling to the following conventions:
- ident style:
    * Please avoid using tabs, use spaces instead of tabs.
- comments:
    * Comments start with `//`.
    * For exposing documentation to Doxygen start comments with `///`.
- structs:
    * Self-explanatory names: a struct `DoubleEndedQueueStruct` will tell the developer what the struct is used for.
    * Use a capital-letter for the first character of the struct.
    * Please use `typedef` for simplifying struct-naming, e.g. `DoubleEndedQueueStruct` instead of `struct DoubleEndedQueueStruct`.
    * All the public structs provided by the runtime/primitives to the end-user should start with the `Synch` prefix (starting from v.3.0.0), e.g. `SynchMemoryStruct`, etc.
- functions: for functions the following rules hold:
    * Self-explanatory names: a function `getMemory()` will tell the developer what it returns as well as `getThreadId()`, etc.
    * Avoid using a capital-letter for the first character of the function.
    * All the public functions provided by the runtime/primitives to the end-user should start with the `synch` prefix (starting from v.3.0.0), i.e. `synchGetMemory`.
    * All the internal functions should NOT start with the `synch` prefix, i.e. `getMemory`.
- definitions:
    * Use capital letters for all the definitions.
    * It is strongly recommended all the public definitions to start with the `SYNCH_` prefix.
- memory allocation and alignment:
    * For memory management, please do not use `malloc`, `calloc`, etc. functions directly.
    * Please use the memory management functionality provided by `primitives.h`.
    * In case that the functionality of `primitives.h` does not meet the needs of your patch, please consider to expand it.
- atomic primitives and other processor primitives:
    * Do not use `asm` assembly code or compiler intrinsics directly in your code.
    * Please use the primitives provided by `primitives.h`. 
    * In case that the functionality of `primitives.h` does not meet the needs of your patch, please consider to expand it.
- file-naming conventions:
    * The source-code of each concurrent data-structure that is provided is placed under the `libconcurrent/concurrent` directory.
    * Each concurrent data-structure provides a public API that is a `.h` header-file placed under the `libconcurrent/includes` directory.
    * For each concurrent data-structure, at least one benchmark is provided under the `benchmarks` directory.
    * All benchmarks are placed under the benchmarks directory. The source-code of each of them has the `bench.c` suffix, while the prefix is the name of the concurrent data-structure that is benchmarked (i.e. the `benchmarks/ccsynchbench.c` file is a benchmark for the CC-Synch combining object).

## API & Code documentation
- It is strongly recommended to sufficiently comment your code.
- At least, ensure that all the public functionality of your patch provides adequate documentation through Doxygen.

## Basic correctness validation
Whenever the `DEBUG` macro is defined in `libconcurrent/config.h`, most of the provided benchmarks make some basic sanity-checks in order to ensure that the provided concurrent data-structures behave appropriately. For example, in concurrent stack benchmarks where each thread executes a specific amount of Push/Pop pairs of operations on a concurrent stack, the `DEBUG` macro ensures that after the end of the benchmark the stack is empty of elements (given that the amount of elements in the initialization phase of the benchmark is zero). 

The `validate.sh` script compiles the sources in `DEBUG` mode and runs a big set of benchmarks with various numbers of threads. After running each of the benchmarks, the script evaluates the `DEBUG` output and in case of success it prints `PASS`. In case of a failure, the script simply prints `FAIL`.

In case that you want to contribute a new concurrent data-structure, you should provide sanity-check code in a benchmark that should be activated whenever `DEBUG` macro is defined. Please, don't forget to update the `validate.sh` script appropriately.

## Final steps before pulling a request
- Update the README.md with details.
- Log your contribution/changes to the CHANGELOG.md file.
- Increase the version numbers (if needed) in any example-files and the README.md to the new version that this pull Request would represent. The versioning scheme that we use (for releases later than v2.2.0) is [SemVer](https://semver.org/). In short, given a version number MAJOR.MINOR.PATCH, increment the:
    1. MAJOR version when you make incompatible API changes,
    2. MINOR version when you add functionality in a backwards compatible manner, and
    3. PATCH version when you make backwards compatible bug fixes.

# Expected Performance

This document provides a brief discussion of the expected performance for some concurrent objects provided by the Synch framework. For a more detailed performance analysis of the provided objects, the reader is encouraged to look at [1,2,3,10].

## Evaluation testbed

The provided brief performance analysis is performed  in a 64-core AMD Epyc multiprocessor (abbr. Epyc) consisting of 2 Epyc 7501 processors that is able to handle 128 threads. Each of these processors consists of 4 silicon dies and each silicon die consist of a NUMA node. Each die contains 8 processing cores, where each core is handling 2 threads. Thus, the AMD multiprocessor consists of 8 NUMA nodes, where each NUMA node consists of 8 processing cores.

## Combining objects performance

At the first experiment, the CC-Synch [1], DSM-Synch [1], H-Synch [1] and PSim [2] combining objects are compared with a universal object based on MCS spin-locks [9] (provided by the `mcsbench.run` benchmark), a simple lock-free implementation (provided by the `lfuobjectbench.run` benchmark), and the flat-combining (abbr. FC) technique provided by [11, 12].

In the experiment of the following figure, the performance of the evaluated combining objects is presented while performing the Fetch&Multiply synthetic benchmark. In this benchmark, a simple concurrent object that supports an atomic Fetch&Multiply float object is implemented. This object is simple enough to exhibit any overheads that a synchronization technique may induce while simulating a very simple concurrent object under high contention. The figure shows the average throughput (Fetch&Multiply operations per second) that each synchronization technique achieves when it executes 10^7 Fetch&Multiply operations (i.e., by setting `-r=10000000`). Since for different values of threads t, each thread executes  10^7/t Fetch&Multiply operations, this experiment measures strong scaling. Specifically, the horizontal axis represents the number of threads t, while the vertical axis displays the throughput (in millions of operations per second) that each synchronization technique achieves. For each value of t, the experiment has been performed 10 times and averages have been calculated using the `bench.sh` script. A random number of dummy loop iterations (i.e, by setting `-w=512`) have been inserted between the execution of two consecutive Fetch&Multiply operations by each thread. In this way, a random work load large enough to avoid unrealistically low cache miss ratios and long runs is simulated. However, this workload is not big enough to substantially reduce the contention. In the experiments performed on the Epyc machine, H-Synch outperforms all other synchronization techniques since it exploits the hierarchical NUMA nature of the Epyc machine. More specifically, H-Synch achieves up to 3 times higher throughput than flat-combining and outperforms CC-Synch by a factor of up to 1.23. The simple lock-free implementation of Fetch&Multiply is slightly slower than PSim and faster that flat-combining. Also, CC-Synch is up to 6 times faster than the object based on MCS spin-locks. CC-Synch performs also very well; its performance is close to that of DSM-Synch and in par with the lock-free implementation.

<p align="center">
    <img src="resources/fam_epyc_w512.png" width="80%">
</p>

In the experiment illustrated in the following figure, the behavior of the evaluated objects for different amounts of random work (i.e. for different numbers of dummy loop iterations inserted between the executions of two Fetch&Multiply) is studied. We fix the number of threads to 128 and we perform the experiment for different random work values (0−32k). It is shown that for a wide range of values (64−1024), there are no big differences on the throughput exhibited by the evaluated algorithms. The reason for this is that for all these values the synchronization cost is the dominant performance factor. In cases that the random work is too high (greater than 16k), the throughput of all algorithms degrades and the performance differences among them become minimal since the amount of random work becomes then the dominant performance factor.

<p align="center">
    <img src="resources/fam_epyc_local_work.png" width="80%">
</p>

## Performance of stacks and queues

The figure below illustrates the performance of some concurrent stack implementations provided by the Synch framework. More specifically, this figure presents the performance of H-Stack [1], DSM-Stack [1], CC-Stack [1], the CLH-Stack that is a simple stack implementation based on CLH queue-locks [5,6], the wait-free SimStack [2], the lock-free stack of [7], and FC-Stack that is a stack implementation based on flat-combining [11,12]. Each of the t threads executes 10^7/t pairs of push and pop operations starting from an empty data structure. This experiment is performed for different values of t. Similarly to the first experiment, a random local work (up to 512 dummy loop iterations) is simulated between the execution of two consecutive operations by the same thread. The stack was initially empty. H-Stack leverages the hierarchical NUMA characteristics of the Epyc machine and it outperforms CC-Stack by a factor of up to 1.25, while being up to 6 times faster than FC-Stack. CC-Stack outperforms DSM-Queue by being 1.16 times faster that SimStack.

<p align="center">
    <img src="resources/stacks_epyc_w512.png" width="80%">
</p>

Similarly to the previous experiment for stacks, the following figure shows the performance of some concurrent queue implementations provided by the Synch framework. More specifically, this figure presents the performance of H-Queue [1], DSM-Queue [1], CC-Queue [1], the CLH-Queue that is the two lock queue implementation by Michael and Scott [5,6,7], the wait-free SimQueue [1], the lock free queue implementation presented by Michael and Scott in [7], and FC-Queue that is a queue implementation based on flat-combining [11,12]. Each of the t threads executes 10^7/t pairs of enqueue and dequeue operations, starting from an empty data structure. This experiment is performed for different values of t. Similarly to previous experiments, a random local work (up to 512 dummy loop iterations) is simulated between the execution of two consecutive operations by the same thread. The queue was initially empty. H-Queue, which is a queue implementation based on H-Synch, leverages the hierarchical nature of the Epyc machine and outperforms all the other queue implementations. More specifically, H-Queue outperforms CC-Queue by a factor of up to 1.24, while being up to 3 times faster than FC-Queue. CC-Queue and DSM-Queue achieve very similar performance while being faster than SimQueue.

<p align="center">
    <img src="resources/queues_epyc_w512.png" width="80%">
</p>

# References

[1]. Panagiota Fatourou, and Nikolaos D. Kallimanis. "Revisiting the combining synchronization technique". ACM SIGPLAN Notices. Vol. 47. No. 8. ACM, PPoPP 2012.

[2]. Panagiota Fatourou, and Nikolaos D. Kallimanis. "A highly-efficient wait-free universal construction". Proceedings of the twenty-third annual ACM symposium on Parallelism in algorithms and architectures (SPAA), 2011.

[3]. Panagiota Fatourou, and Nikolaos D. Kallimanis. "Lock Oscillation: Boosting the Performance of Concurrent Data Structures" Proceedings of the 21st International Conference on Principles of Distributed Systems (Opodis), 2017.

[4]. Yoshihiro Oyama, Kenjiro Taura, and Akinori Yonezawa. "Executing parallel programs with synchronization bottlenecks efficiently". Proceedings of the International Workshop on Parallel and Distributed Computing for Symbolic and Irregular Applications. Vol. 16. 1999.

[5]. Travis S. Craig. "Building FIFO and priority-queueing spin locks from atomic swap". Technical Report TR 93-02-02, Department of Computer Science, University of Washington, February 1993.

[6]. Peter Magnusson, Anders Landin, and Erik Hagersten. "Queue locks on cache coherent multiprocessors". Parallel Processing Symposium, 1994. Proceedings., Eighth International. IEEE, 1994.

[7]. Maged M. Michael, and Michael L. Scott. "Simple, fast, and practical non-blocking and blocking concurrent queue algorithms". Proceedings of the fifteenth annual ACM symposium on Principles of distributed computing. ACM, 1996.

[8]. R. Kent Treiber. "Systems programming: Coping with parallelism". International Business Machines Incorporated, Thomas J. Watson Research Center, 1986.

[9]. John M. Mellor-Crummey, and Michael L. Scott. "Algorithms for scalable synchronization on shared-memory multiprocessors". ACM Transactions on Computer Systems (TOCS) 9.1 (1991): 21-65.

[10]. Panagiota Fatourou, and Nikolaos D. Kallimanis. "Highly-efficient wait-free synchronization". Theory of Computing Systems 55.3 (2014): 475-520.

[11]. Danny Hendler, Itai Incze, Nir Shavit, and Moran Tzafrir. Flat combining and
the synchronization-parallelism tradeoff. In Proceedings of the 22nd Annual ACM
Symposium on Parallel Algorithms and Architectures, pages 355–364, 2010.

[12]. Danny Hendler, Itai Incze, Nir Shavit, and Moran Tzafrir. The Code for Flat-
Combining. http://github.com/mit-carbon/Flat-Combining.

# Contributor Covenant Code of Conduct

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
nkallima (at) gmail.com.
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
---
name: Bug report
about: Create a report to help us improve
title: "[BUG]"
labels: ''
assignees: ''

---

**Code version**
Version(s) of the code .

**Machine architecture/setup**
Describe your setup, i.e. machine architecture, number of cores/processors, etc.

**Compiler**
Report the compiler you use to compile the sources and its version.

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the bug.

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Additional context**
Add any other context about the problem here.
---
name: Feature request
about: Suggest an idea for this project
title: "[FEATURE]"
labels: enhancement
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
