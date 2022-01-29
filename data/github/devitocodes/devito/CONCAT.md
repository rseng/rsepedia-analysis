# Devito: Fast Stencil Computation from Symbolic Specification

[![Build Status for the Core backend](https://github.com/devitocodes/devito/workflows/CI-core/badge.svg)](https://github.com/devitocodes/devito/actions?query=workflow%3ACI-core)
[![Build Status with MPI](https://github.com/devitocodes/devito/workflows/CI-mpi/badge.svg)](https://github.com/devitocodes/devito/actions?query=workflow%3ACI-mpi)
[![Build Status on GPU](https://github.com/devitocodes/devito/workflows/CI-gpu/badge.svg)](https://github.com/devitocodes/devito/actions?query=workflow%3ACI-gpu)
[![Code Coverage](https://codecov.io/gh/devitocodes/devito/branch/master/graph/badge.svg)](https://codecov.io/gh/devitocodes/devito)
[![Slack Status](https://img.shields.io/badge/chat-on%20slack-%2336C5F0)](https://join.slack.com/t/devitocodes/shared_invite/zt-gtd2yxj9-Y31YKk_7lr9AwfXeL2iMFg)
[![asv](http://img.shields.io/badge/benchmarked%20by-asv-blue.svg?style=flat)](https://devitocodes.github.io/devito-performance)
[![PyPI version](https://badge.fury.io/py/devito.svg)](https://badge.fury.io/py/devito)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/devitocodes/devito/master)

[Devito](http://www.devitoproject.org) is a Python package to implement
optimized stencil computation (e.g., finite differences, image processing,
machine learning) from high-level symbolic problem definitions.  Devito builds
on [SymPy](http://www.sympy.org/en/index.html) and employs automated code
generation and just-in-time compilation to execute optimized computational
kernels on several computer platforms, including CPUs, GPUs, and clusters
thereof.

- [About Devito](#about-devito)
- [Installation](#installation)
- [Resources](#resources)
- [Performance](#performance)
- [Get in touch](#get-in-touch)
- [Interactive jupyter notebooks](#interactive-jupyter-notebooks)

## About Devito

Devito provides a functional language to implement sophisticated operators that
can be made up of multiple stencil computations, boundary conditions, sparse
operations (e.g., interpolation), and much more.  A typical use case is
explicit finite difference methods for approximating partial differential
equations. For example, a 2D diffusion operator may be implemented with Devito
as follows

```python
>>> grid = Grid(shape=(10, 10))
>>> f = TimeFunction(name='f', grid=grid, space_order=2)
>>> eqn = Eq(f.dt, 0.5 * f.laplace)
>>> op = Operator(Eq(f.forward, solve(eqn, f.forward)))
```

An `Operator` generates low-level code from an ordered collection of `Eq` (the
example above being for a single equation). This code may also be compiled and
executed

```python
>>> op(t=timesteps)
```

There is virtually no limit to the complexity of an `Operator` -- the Devito
compiler will automatically analyze the input, detect and apply optimizations
(including single- and multi-node parallelism), and eventually generate code
with suitable loops and expressions.

Key features include:

* A functional language to express finite difference operators.
* Straightforward mechanisms to adjust the discretization.
* Constructs to express sparse operators (e.g., interpolation), classic linear
  operators (e.g., convolutions), and tensor contractions.
* Seamless support for boundary conditions and adjoint operators.
* A flexible API to define custom stencils, sub-domains, sub-sampling,
  and staggered grids.
* Generation of highly optimized parallel code (SIMD vectorization, CPU and GPU
  parallelism via OpenMP, multi-node parallelism via MPI, blocking, aggressive
  symbolic transformations for FLOP reduction, etc.).
* Distributed NumPy arrays over multi-node (MPI) domain decompositions.
* Inspection and customization of the generated code.
* Autotuning framework to ease performance tuning.
* Smooth integration with popular Python packages such as NumPy, SymPy, Dask,
  and SciPy, as well as machine learning frameworks such as TensorFlow and
  PyTorch.

## Installation

The easiest way to try Devito is through Docker using the following commands:
```
# get the code
git clone https://github.com/devitocodes/devito.git
cd devito

# start a jupyter notebook server on port 8888
docker-compose up devito
```
After running the last command above, the terminal will display a URL such as
`https://127.0.0.1:8888/?token=XXX`. Copy-paste this URL into a browser window
to start a [Jupyter](https://jupyter.org/) notebook session where you can go
through the [tutorials](https://github.com/devitocodes/devito/tree/master/examples)
provided with Devito or create your own notebooks.

[See here](http://devitocodes.github.io/devito/download.html) for detailed installation
instructions and other options. If you encounter a problem during installation, please
see the
[installation issues](https://github.com/devitocodes/devito/wiki/Installation-Issues) we
have seen in the past. 

## Resources

To learn how to use Devito,
[here](https://github.com/devitocodes/devito/blob/master/examples) is a good
place to start, with lots of examples and tutorials.

The [website](https://www.devitoproject.org/) also provides access to other
information, including documentation and instructions for citing us.

Some FAQ are discussed [here](https://github.com/devitocodes/devito/wiki/FAQ).

## Performance

If you are interested in any of the following

* Generation of parallel code (CPU, GPU, multi-node via MPI);
* Performance tuning;
* Benchmarking operators;

then you should take a look at this
[README](https://github.com/devitocodes/devito/blob/master/benchmarks/user).

You may also be interested in
[TheMatrix](https://github.com/devitocodes/thematrix) -- a cross-architecture
benchmarking framework showing the performance of several production-grade
seismic operators implemented with Devito. This is now our flagship project
towards neat, open, and reproducible science.

## Get in touch

If you're using Devito, we would like to hear from you. Whether you
are facing issues or just trying it out, join the
[conversation](https://join.slack.com/t/devitocodes/shared_invite/zt-gtd2yxj9-Y31YKk_7lr9AwfXeL2iMFg).

## Interactive jupyter notebooks
The tutorial jupyter notebook are available interactively at the public [binder](https://mybinder.org/v2/gh/devitocodes/devito/master) jupyterhub. 
# Contributor Covenant Code of Conduct

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
reported by contacting the project team as a group
 [devito-maint](mailto:devito-maint@imperial.ac.uk)
or confidentially by contacting
 [Gerard Gorman](mailto:g.gorman@imperial.ac.uk),
 and/or [Fabio Luporini](mailto:f.luporini12@imperial.ac.uk).
All complaints will be reviewed and investigated and will result in a response
that is deemed necessary and appropriate to the circumstances. The project team
is obligated to maintain confidentiality with regard to the reporter of an
incident. Further details of specific enforcement policies may be posted
separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see
https://www.contributor-covenant.org/faq

# Contributing to Devito

We welcome third-party contributions, and we would love you to become an active contributor!

Software contributions are made via GitHub pull requests to https://github.com/devitocodes/devito. If you are planning a large contribution, we encourage you to engage with us frequently to ensure that your effort is well-directed. See below for more details.

Devito is distributed under the MIT License, https://github.com/devitocodes/devito/blob/master/LICENSE.md. The act of submitting a pull request or patch (with or without an explicit Signed-off-by tag) will be understood as an affirmation of the following:

 Developer's Certificate of Origin 1.1

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

### Reporting issues

There are several options:
* Talk to us. You can join our Slack team via this [link](https://join.slack.com/t/devitocodes/shared_invite/zt-gtd2yxj9-Y31YKk_7lr9AwfXeL2iMFg). Should you have installation issues, or should you bump into something that appears to be a Devito-related bug, do not hesitate to get in touch. We are always keen to help out.
* File an issue on [our GitHub page](https://github.com/devitocodes/devito/issues).

### Making changes

First of all, read of [code of conduct](https://github.com/devitocodes/devito/blob/master/CODE_OF_CONDUCT.md) and make sure you agree with it.

The protocol to propose a patch is:
* [Recommended, but not compulsory] Talk to us on Slack about what you're trying to do. There is a great chance we can support you.
* As soon as you know what you need to do, [fork](https://help.github.com/articles/fork-a-repo/) Devito.
* Create a branch with a suitable name.
* Write code following the guidelines below. Commit your changes as small logical units.
* Commit messages should adhere to the format `<tag>: <msg>`, where `<tag>` could be, for example, "ir" (if the commit impacts the intermediate representation), "operator", "tests", etc. We may ask you to rebase the commit history if it looks too messy.
* Write tests to convince us and yourself that what you've done works as expected. Commit them.
* Run **the entire test suite**, including the new tests, to make sure that you haven't accidentally broken anything else.
* Push everything to your Devito fork.
* Submit a Pull Request on our repository.
* Wait for us to provide feedback. This may require a few iterations.

Tip, especially for newcomers: prefer short, self-contained Pull Requests over lengthy, impenetrable, and thus difficult to review, ones.

### Coding guidelines

Some coding rules are "enforced" (and automatically checked by our Continuous Integration systems), some are "strongly recommended", others are "optional" but welcome.

* We _enforce_ [PEP8](https://www.python.org/dev/peps/pep-0008/), with a few exceptions, listed [here](https://github.com/devitocodes/devito/blob/master/setup.cfg#L3)
* We _enforce_ a maximum line length of 90 characters.
* We _enforce_ indentation via 4 spaces.
* We _suggest_ to use ``flake8`` to check the above points locally, before filing a Pull Request.
* We _strongly recommend_ to document any new module, class, routine, ... with [NumPy-like docstrings](https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html#example-numpy) ("numpydoc").
* We _strongly recommend_ imports to be at the top of a module, logically grouped and, within each group, to be alphabetically ordered. As an example, condider our [__init__.py](https://github.com/devitocodes/devito/blob/master/devito/__init__.py): the first group is imports from the standard library; then imports from third-party dependencies; finally, imports from devito modules.
* We _strongly recommend_ to follow standard Python coding guidelines:
  - Use camel caps for class names, e.g. ``class FooBar``.
  - Method names must start with a small letter; use underscores to separate words, e.g. ``def _my_meth_...``.
  - Private class attributes and methods must start with an underscore.
  - Variable names should be explicative (Devito prefers "long and clear" over "short but unclear").
  - Comment your code, and do not be afraid of being verbose. The first letter must be capitalized. Do not use punctuation (unless the comment consists of multiple sentences).
* We _like_ that blank lines are used to logically split blocks of code implementing different (possibly sequential) tasks.

### Adding tutorials or examples

We always look forward to extending our [suite of tutorials and examples](https://www.devitoproject.org/devito/tutorials.html) with new Jupyter Notebooks. Even something completely new, such as a new series of tutorials showing your work with Devito, would be a great addition.
# Binder directory

This directory contains a ```requirements.txt``` file with the added requirement of **matplotlib**, thus enabling jupyter notebooks hosted by mybinder.org to display graphics.

Use the binder badge below to start an interactive jupyterhub server.

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/devitocodes/devito/master)# Benchmarking Devito

There are two subdirectories:

* `benchmarks/user`, which provides a Python script, `benchmark.py`, to
  evaluate the performance of some relevant Operators defined in `/examples`.
  `python benchmark.py --help` explains how to configure a benchmark run.
  Users interested in benchmarking Devito may want to explore this approach.
* `benchmarks/regression`, which provides a series of performance regression
  tests. This is used by continuous integration, to ensure that new features
  to be merged into trunk do not cause slowdowns. The performance regression
  framework is based on [airspeed velocity](https://asv.readthedocs.io/en/stable/).
# Benchmarking a Devito Operator

## Running benchmarks

`benchmark.py` implements a minimalist framework to evaluate the performance of
a Devito Operator while varying:

* the problem size (e.g., shape of the computational grid);
* the discretization (e.g., space- and time-order of the input/output fields);
* the simulation time (in milliseconds);
* the performance optimization level;
* the autotuning level.

Running `python benchmark.py --help` will display a list of useful options.

## A step back: configuring your machine for reliable benchmarking

If you are tempted to use your laptop to run a benchmark, you may want to
reconsider: heat and power management may affect the results you get in an
unpredictable way.

It is important that *both* the Python process running Devito (process*es* if
running with MPI) and the OpenMP threads spawned while running an Operator are
pinned to specific CPU cores, to get reliable and deterministic results. There
are several ways to achieve this:

* Through environment variables. All MPI/OpenMP distributions provide a set of
  environment variables to control process/thread pinning.  Devito also
  supplies the `set_omp_pinning.sh` program (under `/scripts`), which helps
  with thread pinning (though, currently, only limited to Intel architectures).
* Through a program such as `numactl` or `taskset`.

If running on a NUMA system, where multiple nodes of CPU cores ("sockets") and
memory are available, pinning becomes even more important (and therefore
deserves more attention). On a NUMA system, a core can access both local and
remote memory nodes, but the latency is obviously smaller in the former case.
Thus, if a process/thread is pinned to a given core, it is important that as
much accessed data as possible is allocated in local memory (ideally, the entire
working set). There are multiple scenarios that are worth considering:

* Purely sequential run (no OpenMP, no MPI). Use `numactl` or `taskset` to pin
  the Python process. On a NUMA system, this also ensures that all data gets
  allocated in local memory.
* OpenMP-only run. When encountering a parallel region, an Operator will spawn a
  certain number of OpenMP threads. By default, as many threads as available
  logical cores are created. This can be changed by setting the OpenMP-standard
  `OMP_NUM_THREADS` environment variable to a different value. When might we
  want to do this?
  - Unless on a hyperthreads-centerd system, such as an Intel Knights Landing,
    spawning only as many threads as *physical* cores usually results in
    slightly better performance due to less contention for hardware resources.
  - Since, here, we are merely interested in benchmarking, when running on a
    NUMA system we should restrain an Operator to run on a single node of CPU
    cores (i.e., a single socket), as in practice one should always use MPI
    across multiple sockets (hence resorting to the MPI+OpenMP mode). There is
    also one caveat: we need to make sure that the Python process runs on the
    same socket as the OpenMP threads, otherwise some data might get allocated
    on remote memory. For this, `numactl` or `taskset` are our friends.
* MPI-only run. Process pinning should be implemented exploiting the proper MPI
  environment variables.
* MPI+OpenMP. The typical execution mode is: one MPI process per socket, and
  each MPI process spawns a set of OpenMP threads upon entering an OpenMP
  parallel region.  Pinning is typically enforced via environment variables.

Some more information about pinning is available
[here](https://hpc-wiki.info/hpc/Binding/Pinning).

There are many ways one can check that pinning is working as expected. A
recommended tool for rapid visual inspection is [htop](http://hisham.hm/htop/).

## Enabling OpenMP

To switch on multi-threaded execution of Operators via OpenMP, the following
environment variable must be set:
```
DEVITO_LANGUAGE=openmp
```
One has two options: either set it explicitly or prepend it to the Python
command. In the former case, assuming a bash shell:
```
export DEVITO_LANGUAGE=openmp
```
In the latter case:
```
DEVITO_LANGUAGE=openmp python benchmark.py ...
```

## Enabling MPI

To switch on MPI, one should set
```
DEVITO_MPI=1
```
and run with `mpirun -n number_of_processes python benchmark.py ...`

Devito supports multiple MPI schemes for halo exchange. See the `Tips` section
below.

## The optimization level

`benchmark.py` allows to set optimization mode, as well as several optimization
options, via the `--opt` argument. Please refer to
[this](https://github.com/devitocodes/devito/blob/master/examples/performance/00_overview.ipynb)
notebook for a comprehensive list of all optimization modes and options
available in Devito. You may also want to take a look at the example command
lines a few sections below.

## Auto-tuning

Auto-tuning can significantly improve the run-time performance of an Operator. It
can be enabled on an Operator basis:
```
op = Operator(...)
op.apply(autotune=True)
```
The auto-tuner will discover a suitable block shape for each blocked loop nest
in the generated code.

With `autotune=True`, the auto-tuner operates in `basic` mode, which only attempts
a small batch of block shapes. With `autotune='aggressive'`, the auto-tuning phase
will likely take up more time, but it will also evaluate more block shapes.

By default, `benchmark.py` runs Operators with auto-tuning in aggressive mode,
that is as `op.apply(autotune='aggressive')`. This can be changed with the
`-a/--autotune` flags. In particular, `benchmark.py` uses the so called
"pre-emptive" auto-tuning, which implies two things:

* The Operator's output fields are copied, and the Operator will write to these
  copies while auto-tuning. So the memory footprint is temporarily larger during
  this phase.
* The auto-tuning phase produces values that are eventually ditched;
  afterwards, the actual computation takes place. The execution time of the
  latter does not include auto-tuning.

Note that in production runs one should/would rather use the so called "runtime
auto-tuning":
```
op.apply(autotune=('aggressive', 'runtime'))
```
in which auto-tuning, as the name suggests, will be performed during the first N
timesteps of the actual computation, after which the best block shapes are
selected and used for all remaining timesteps.

## Choice of the backend compiler

The "backend compiler" takes as input the code generated by Devito and
translates it into a shared object. Supported backend compilers are `gcc`,
`icc`, `pgcc`, `clang`. For each of these compilers, Devito uses some preset compilation
flags (e.g., -O3, -march=native, etc).

The default backend compiler is `gcc`. To change it, one should set the
`DEVITO_ARCH` environment variable to a different value; run
```
from devito import print_defaults
print_defaults()
```
to get all possible `DEVITO_ARCH` values.

## Benchmark verbosity

Run with `DEVITO_LOGGING=DEBUG` to find out the specific performance
optimizations applied by an Operator, how auto-tuning is getting along, and to
emit more performance metrics.

## Tips

* The most powerful MPI mode is called "full", and is activated setting
  `DEVITO_MPI=full` instead of `DEVITO_MPI=1`. The "full" mode supports
  computation/communication overlap.
* When auto-tuning is enabled, one should always run in performance mode:
  ```
  from devito import mode_performance
  mode_perfomance()
  ```
  This is automatically turned on by `benchmark.py`

## Example commands

The isotropic acoustic wave forward Operator in a `512**3` grid, space order
12, and a simulation time of 100ms:
```
python benchmark.py run -P acoustic -d 512 512 512 -so 12 --tn 100
```
Like before, but with auto-tuning in `basic` mode:
```
python benchmark.py run -P acoustic -d 512 512 512 -so 12 -a basic --tn 100
```
It is also possible to run a TTI forward operator -- here in a 512x402x890
grid:
```
python benchmark.py run -P tti -d 512 402 890 -so 12 -a basic --tn 100
```
Same as before, but telling devito not to use temporaries to store the
intermediate values which stem from mixed derivatives:
```
python benchmark.py run -P tti -d 512 402 890 -so 12 -a basic --tn 100 --opt
"('advanced', {'cire-mingain: 1000000'})"
```
Do not forget to pin processes, especially on NUMA systems; below, we use
`numactl` to pin processes and threads to one specific NUMA domain.
```
numactl --cpubind=0 --membind=0 python benchmark.py ...
```
While a benchmark is running, you can have some useful programs running in
background in other shells. For example, to monitor pinning:
```
htop
```
or to keep the memory footprint under control:
```
watch numastat -m
```

## The run-jit-backdoor mode

As of Devito v3.5 it is possible to customize the code generated by Devito.
This is often referred to as the ["JIT backdoor"
mode](https://github.com/devitocodes/devito/wiki/FAQ#can-i-manually-modify-the-c-code-generated-by-devito-and-test-these-modifications).
With ``benchmark.py`` we can exploit this feature to manually hack and test the
code generated for a given benchmark. So, we first run a problem, for example
```
python benchmark.py run-jit-backdoor -P acoustic -d 512 512 512 -so 12 --tn 100
```
As you may expect, the ``run-jit-backdoor`` mode accepts exactly the same arguments
as the ``run`` mode. Eventually, you will see a message along the lines of
```
You may now edit the generated code in
`/tmp/devito-jitcache-uid1000/31e8d25408f369754e2b7a26f4439944dc7683e2.c`. Then
save the file, and re-run this benchmark.
```
At this point, just follow the instructions on screen. The next time you run
the benchmark, the modified C code will be re-compiled and executed. Thus,
you will see the performance impact of your changes.

## Running on HPC clusters

`benchmark.py` can be used to evaluate MPI on multi-node systems:
```
mpiexec python benchmark.py ...
```
In `bench` mode, each MPI rank will produce a different `.json` file
summarizing the achieved performance in a structured format.

Further, we provide `make-pbs.py`, a simple program to generate PBS files
to submit jobs on HPC clusters. Take a look at `python make-pbs.py --help`
for more information, and in particular `python make-pbs.py generate --help`.
`make-pbs.py` is especially indicated if interested in running strong scaling
experiments.

## Benchmark output

The GFlops/s and GPoints/s performance, Operational Intensity (OI) and
execution time are emitted to standard output at the end of each run.  You may
find this
[FAQ](https://github.com/devitocodes/devito/wiki/FAQ#how-does-devito-compute-the-performance-of-an-operator)
useful.
Example runs:

* `python3 run_advisor.py --name isotropic --path <path-to-devito>/examples/seismic/acoustic/acoustic_example.py`
* `python3 run_advisor.py --name tti_so8 --path <path-to-devito>/examples/seismic/tti/tti_example.py --exec-args "-so 8"`
* `python3 run_advisor.py --name iso_ac_so6 --path <path-to-devito>/benchmarks/user/benchmark.py --exec-args "bench -P acoustic -so 6 --tn 200 -d 100 100 100 --autotune off -x 1"`

After the run has finished you should be able to plot a roofline with the results and export roofline data to JSON using:
* `python3 roofline.py --name Roofline --project <advisor-project-name>`

To create a read-only snapshot for use with Intel Advisor GUI, use:
* `advixe-cl --snapshot --project-dir=<advisor-project-name> pack -- /<new-snapshot-name>`

Prerequisites:
* Support guaranteed only for Intel Advisor as installed with Intel Parallel Studio v 2020 Update 2
  and Intel oneAPI 2021; earlier years may not work; other 2020/2021 versions, as well as later years,
  may or may not work.
* In Linux systems you may need to enable system-wide profiling by setting:
  - `/proc/sys/kernel/yama/ptrace_scope` to `0`
  - `/proc/sys/kernel/perf_event_paranoid` to `1`

* `numactl` must be available on the system. If not available, install with:
	`sudo apt-get install numactl`
* Install `pandas` and `matplotlib`. They are not included in the core Devito installation.

Limitations:

* Untested with more complicated examples.
* Untested on Intel KNL (we might need to ask `numactl` to bind to MCDRAM).
* Running the `tripcounts` analysis takes a lot, despite starting in paused
  mode. This analysis, together with the `survey` analysis, is necessary to
  generate a roofline. Both are run by `run_advisor.py`.
* Requires python3, untested in earlier versions of python and conda environments
* Currently requires download of repository and running `pip3 install .`, the scripts
  are currently not included as a package with the user installation of Devito

TODO:

* Give a name to the points in the roofline, otherwise it's challenging to
  relate loops (code sections) to data.
* Emit a report summarizing the configuration used to run the analysis
  (threading, socket binding, ...).

Useful links:
* [ Memory-Level Roofline Analysis in Intel® Advisor ](https://software.intel.com/content/www/us/en/develop/articles/memory-level-roofline-model-with-advisor.html " Memory-Level Roofline Analysis in Intel® Advisor ")
* [CPU / Memory Roofline Insights
Perspective](https://software.intel.com/content/www/us/en/develop/documentation/advisor-user-guide/top/optimize-cpu-usage/cpu-roofline-perspective.html "CPU / Memory Roofline Insights
Perspective")
* [ Roofline Resources for Intel® Advisor Users ](https://software.intel.com/content/www/us/en/develop/articles/advisor-roofline-resources.html " Roofline Resources for Intel® Advisor Users ")# Performance regression

Based on [airspeed velocity](https://asv.readthedocs.io/en/stable/).
## How to navigate this directory

Examples and tutorials are provided in the form of single Python files, as Jupyter
notebooks, or as mini-apps built on top of Devito.

Jupyter notebooks are files with extension `.ipynb`. To execute these, run
`jupyter notebook`, and then click on the desired notebook in the window that
pops up in your browser. In alternative, you may explore the pre-rendered
notebooks directly on GitHub or, for a potentially smoother experience, [with
nbviewer](https://nbviewer.jupyter.org/github/devitocodes/devito/tree/master/examples/).

We recommend newcomers to start with the following sets of tutorials:

* `userapi`: Gentle introduction to symbolic computation with Devito.
* `cfd`: A series of introductory notebooks showing how to use Devito to
  implement finite difference operators typical of computational fluid
  dynamics. These are based on the excellent blog ["CFD Python:12 steps to
  Navier-Stokes"](http://lorenabarba.com/blog/cfd-python-12-steps-to-navier-stokes/)
  by the Lorena A. Barba group.

A set of more advanced examples are available in `seismic`:

* `seismic/tutorials`: A series of Jupyter notebooks of incremental complexity,
  showing a variety of Devito features in the context of seismic inversion
  operators. Among the discussed features are custom stencils, staggered
  grids, tensor notation, and time blocking.
* `seismic/acoustic`: Example implementations of isotropic acoustic forward,
  adjoint, gradient and born operators, suitable for full-waveform inversion
  methods (FWI).
* `seismic/tti`: Example implementations of several anisotropic acoustic
  forward operators (TTI).
* `seismic/elastic`: Example implementation of an isotropic elastic forward
  operator. `elastic`, unlike `acoustic` and `tti`, fully exploits the
  tensorial nature of the Devito symbolic language.
* `seismic/viscoelastic`: Example implementation of an isotropic viscoelastic
  forward operator. Like `elastic`, `viscoelastic` exploits tensor functions
  for a neat and compact representation of the discretized partial differential
  equations.
* `seismic/self-adjoint`: Self-adjoint energy conserving pseudo-acoustic
  operators, including notebooks for implementation of the nonlinear forward,
  the forward and adjoint linearized Jacobian, and tests proving accuracy and 
  correctness.

Further:

* `mpi`: Jupyter notebooks explaining how MPI works in Devito.
* `finance`: Jupyter notebooks with examples of applying Devito to partial differential equations with financial applications.
* `misc`: Example operators outside the context of finite differences and
* `performance`: Jupyter notebooks explaining the optimizations applied by Devito, the options available to steer the optimization process, how to run on GPUs, and much more.

For developers:

* `compiler`: A set of notebooks exploring the architecture of the Devito
  compiler. This is still in its infancy.

## More resources

* Articles, presentations, posters and much more concerning Devito is available
  [here](https://www.devitoproject.org/publications). The entries are ordered
  chronologically -- those at the top being the most recent ones, for each
  section.
* The user documentation is available [here](http://devitocodes.github.io/devito/).
### Summary of optimization parameters

* Levels
  * `noop`: disable optimizations
  * `advanced`: all optimizations
  * `advanced-fsg`: alternative optimization pipeline

* Options (type, default)
  * Parallelism:
    * `openmp` (boolean, False): enable/disable OpenMP parallelism
    * `par-collapse-ncores` (int, 4): control loop collapsing
    * `par-collapse-work` (int, 100): control loop collapsing
    * `par-chunk-nonaffine` (int, 3): control chunk size in nonaffine loops
    * `par-dynamic-work` (int, 10): switch between dynamic and static scheduling
    * `par-nested` (int, 2): control nested parallelism
  * Blocking:
    * `blockinner` (boolean, False): enable/disable loop blocking along innermost loop
    * `blocklevels` (int, 1): 1 => classic loop blocking; 2 for two-level hierarchical blocking; etc.
  * CIRE:
    * `min-storage` (boolean, False): smaller working set size, less loop fusion
    * `cire-rotate` (boolean, False): smaller working set size, fewer parallel dimensions
    * `cire-maxpar` (boolean, False): bigger working set size, more parallelism
    * `cire-ftemps` (boolean, False): give user control over the allocated temporaries
    * `cire-mingain` (int, 10): minimum gain to optimize away a redundant expression
    * `cire-schedule` ((str, int), 'automatic'): scheduling strategy for derivatives
  * Device-specific:
    * `gpu-fit` (boolean, False): list of saved TimeFunctions that fit in the device memory
    * `par-disabled` (boolean, True): enable/disable parallelism on the host
  * Misc:
    * `linearize` (boolean, False): linearize array accesses


### Optimization parameters by platform

* Parallelism

|                     |        CPU          |         GPU        |
|---------------------|---------------------|--------------------|
| openmp              | :heavy_check_mark:  | :heavy_check_mark: |
| par-collapse-ncores | :heavy_check_mark:  |         :x:        |
| par-collapse-work   | :heavy_check_mark:  |         :x:        |
| par-chunk-nonaffine | :heavy_check_mark:  | :heavy_check_mark: |
| par-dynamic-work    | :heavy_check_mark:  |         :x:        |
| par-nested          | :heavy_check_mark:  |         :x:        |

* Blocking

|                     |        CPU          |         GPU        |
|---------------------|---------------------|--------------------|
| blockinner          | :heavy_check_mark:  |         :x:        |
| blocklevels         | :heavy_check_mark:  |         :x:        |

* CIRE

|                     |        CPU          |         GPU        |
|---------------------|---------------------|--------------------|
| min-storage         | :heavy_check_mark:  |         :x:        |
| cire-rotate         | :heavy_check_mark:  |         :x:        |
| cire-maxpar         | :heavy_check_mark:  | :heavy_check_mark: |
| cire-ftemps         | :heavy_check_mark:  | :heavy_check_mark: |
| cire-mingain        | :heavy_check_mark:  | :heavy_check_mark: |
| cire-schedule       | :heavy_check_mark:  | :heavy_check_mark: |

* Device-specific

|                     |        CPU          |         GPU        |
|---------------------|---------------------|--------------------|
| gpu-fit             |        :x:          | :heavy_check_mark: |
| par-disabled        |        :x:          | :heavy_check_mark: |

* Misc

|                     |        CPU          |         GPU        |
|---------------------|---------------------|--------------------|
| linearize           | :heavy_check_mark:  | :heavy_check_mark: |
# Absorbing boundary conditions for the acoustic wave equation  - Applications in seismic problems.

## Authors: Felipe A. G. Silva, Saulo R. M. Barros and Pedro S. Peixoto
    Institute of Mathematics and Statistics - Applied Mathematics Department
   (felipe.augusto.guedes@gmail.com, saulo@ime.usp.br, pedrosp@ime.usp.br)

**Important Informations:** These notebooks are part of the Project Software Technologies for Modeling and Inversion (STMI) at RCGI in the  University of Sao Paulo. 

The objective of these notebooks is to present several schemes which are designed to reduce artificial reflections on boundaries in the numerical solution of the acoustic wave equation with finite differences. We consider several methods, covering absorbing boundary conditions and absorbing boundary layers. Among the schemes, we have implemented:

- Clayton's A1 and A2 boundary conditions;
- Higdon second order boundary conditions;
- Sochaki's type of Damping boundary layer;
- Perfectly Matched Layer (PML);
- Hybrid Absorbing Boundary Conditions (HABC);

The computational implementation of the methods above is done within the framework of <a href="https://www.devitoproject.org/">Devito</a>, which is aimed to produce highly optimized code for finite differences discretizations, generated from high level symbolic problem definitions. Devito presents a work structure in Python and generates code in C ++, which can be taylored for high performance on different computational platforms. The notebooks are organized as follows:

- <a href="01_introduction.ipynb">1. Introduction and description of the acoustic problem;</a>
- <a href="02_damping.ipynb">2. Implementation of Sochaki's damping;</a>
- <a href="03_pml.ipynb">3. PML implementation;</a>
- <a href="04_habc.ipynb">4. HABC (Hybrid absorbing boundary conditions. These encompass also the absorbing boundary conditions A1, A2 and Higdon).;</a>

The notebooks bring a theoretical description of the methods together with the Devito implementation, which can be used for  the simulations of interest. We choose a reference problem, described in the notebook <a href="1_introduction.ipynb">Introduction to Acoustic Problem</a>. The spatial and temporal discretizations used throughout the notebooks are also presented in this introductory notebook, together with other relevant concepts to be used overall. Therefore, one should first assimilate the contents of this notebook. 

In the remaining notebooks, we incrementally describe several numerical techniques to reduce artificial boundary reflections. It is better to follow the order of the notebooks, since concepts are used afterward. We include simulations demonstrating the use of the methods. By changing some parameters, the user would be able to carry out several tests.
# Devito Self Adjoint modeling operators

## These operators are contributed by Chevron Energy Technology Company (2020)

These operators are based on simplfications of the systems presented in:
<br>**Self-adjoint, energy-conserving second-order pseudoacoustic systems for VTI and TTI media for reverse migration and full-waveform inversion** (2016)
<br>Kenneth Bube, John Washbourne, Raymond Ergas, and Tamas Nemeth
<br>SEG Technical Program Expanded Abstracts
<br>https://library.seg.org/doi/10.1190/segam2016-13878451.1

## Tutorial goal

The goal of this series of tutorials is to generate -- and then test for correctness -- the modeling and inversion capability in Devito for variable density visco- acoustics. We use an energy conserving form of the wave equation that is *self adjoint*, which allows the same modeling system to be used for all for all phases of finite difference evolution required for quasi-Newton optimization:
- **nonlinear forward**, nonlinear with respect to the model parameters
- **Jacobian forward**, linearized with respect to the model parameters 
- **Jacobian adjoint**, linearized with respect to the model parameters

These notebooks first implement and then test for correctness for three types of modeling physics.

| Physics         | Implementation            | Notebook                          |
|:----------------|:--------------------------|:----------------------------------|
| Isotropic       | Nonlinear ops             | [sa_01_iso_implementation1.ipynb] |
| Isotropic       | Linearized ops            | [sa_02_iso_implementation2.ipynb] |
| Isotropic       | Correctness tests         | [sa_03_iso_correctness.ipynb]     |
|-----------------|---------------------------|-----------------------------------|
| VTI Anisotropic | Nonlinear/linearized ops  | [sa_11_vti_implementation.ipynb]  |
| VTI Anisotropic | Correctness tests         | [sa_12_vti_correctness.ipynb]     |
|-----------------|---------------------------|-----------------------------------|
| TTI Anisotropic | Nonlinear/linearized ops  | [sa_21_tti_implementation.ipynb]  |
| TTI Anisotropic | Correctness tests         | [sa_22_tti_correctness.ipynb]     |
|:----------------|:--------------------------|:----------------------------------|

[sa_01_iso_implementation1.ipynb]: sa_01_iso_implementation1.ipynb
[sa_02_iso_implementation2.ipynb]: sa_02_iso_implementation2.ipynb
[sa_03_iso_correctness.ipynb]:     sa_03_iso_correctness.ipynb
[sa_11_vti_implementation.ipynb]: sa_11_vti_implementation.ipynb
[sa_12_vti_correctness.ipynb]:     sa_12_vti_correctness.ipynb
[sa_21_tti_implementation.ipynb]: sa_21_tti_implementation.ipynb
[sa_22_tti_correctness.ipynb]:     sa_22_tti_correctness.ipynb

## Running unit tests
- if you would like to see stdout when running the tests, use
```py.test -c testUtils.py```

## TODO
- [X] Devito-esque equation version of setup_w_over_q
- [ ] figure out weird test failure depending on the order of equations in operator

**Equation order 1**
```
    return Operator([dm_update] + eqn + rec_term, subs=spacing_map,
                    name='IsoJacobianAdjOperator', **kwargs)
```
**Equation order 2**
```
    return Operator(eqn + rec_term + [dm_update], subs=spacing_map,
                    name='IsoJacobianAdjOperator', **kwargs)
```
    - With Equation order 1, all tests pass
    - With Equation order 2, there are different outcomes for tests 
    - Possibly there is a different path chosen through the AST, and different c code is generated?

- [ ] replace the conditional logic in the stencil with comprehension
```
    space_fd = sum([getattr(b * getattr(field, 'd%s'%d.name)(x0=d+d.spacing/2)),
        'd%s'%d.name)(x0=d-d.spacing/2)) for d in field.dimensions[1:]])
```
- [ ] Add memoized methods back to wavesolver.py
- [ ] Add ensureSanityOfFields methods for iso, vti, tti
- [ ] Add timing info via logging for the w_over_q setup, as in initialize_damp
- [ ] Add smoother back to setup_w_over_q method
```
     eqn1 = Eq(wOverQ, val)
     Operator([eqn1], name='WOverQ_Operator_init')()
     # If we apply the smoother, we must renormalize output to [qmin,qmax]
     if sigma > 0:
         smooth = gaussian_smooth(wOverQ.data, sigma=sigma)
         smin, smax = np.min(smooth), np.max(smooth)
         smooth[:] = qmin + (qmax - qmin) * (smooth - smin) / (smax - smin)
         wOverQ.data[:] = smooth
     eqn2 = Eq(wOverQ, w / wOverQ)
     Operator([eqn2], name='WOverQ_Operator_recip')()
```
- [X] Correctness tests
  - [X] Analytic response in the far field
  - [X] Modeling operator linearity test, with respect to source
  - [X] Modeling operator adjoint test, with respect to source
  - [X] Nonlinear operator linearization test, with respect to model/data
  - [X] Jacobian operator linearity test, with respect to model/data
  - [X] Jacobian operator adjoint test, with respect to model/data
  - [X] Skew symmetry test for shifted derivatives

## To save generated code 

```
f = open("operator.c", "w")
print(op, file=f)
f.close()
```