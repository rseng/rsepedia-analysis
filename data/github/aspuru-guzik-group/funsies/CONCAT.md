# funsies
is a python library and execution engine to build reproducible,
fault-tolerant, distributed and composable computational workflows.

- 🐍 Workflows are specified in pure python.
- 🐦 Lightweight with few dependencies.
- 🚀 Easy to deploy to compute clusters and distributed systems.
- 🔧 Can be embedded in your own apps.
- 📏 First-class support for static analysis. Use
  [mypy](http://mypy-lang.org/) to check your workflows!

Workflows are encoded in a [redis server](https://redis.io/) and executed
using the distributed job queue library [RQ](https://python-rq.org/). A hash
tree data structure enables automatic and transparent caching and incremental
computing.

[Source docs can be found
here.](https://aspuru-guzik-group.github.io/funsies/) Some example funsies
scripts can be found in the [recipes folder.](./recipes)

## Installation
Using `pip`, 

```bash
pip install funsies
```

This will enable the `funsies` CLI tool as well as the `funsies` python
module. Python 3.7, 3.8 and 3.9 are supported. To run workflows, you'll need a
Redis server, version 4.x or higher. On Linux Redis can be installed using conda,

```bash
conda install redis
```

pip,

```bash
pip install redis-server
```

or your system package manager. On Mac OSX, Redis can be downloaded using
Homebrew,

```bash
brew install redis
```

(Windows is not supported by Redis, but a third-party package can be obtained
from [this repository](https://github.com/tporadowski/redis). This has **not**
been tested, however.)

## Hello, funsies!
To run workflows, three components need to be connected:

- 📜 a python script describing the workflow
- 💻 a redis server that holds workflows and data
- 👷 worker processes that execute the workflow

funsies is distributed: all three components can be on different computers or
even be connected at different time. Redis is started using `redis-server`,
workers are started using `funsies worker` and the workflow is run using
python.

For running on a single machine, the `start-funsies` script takes care of starting the database and workers,

```bash
start-funsies \
    --no-pw \
    --workers 2
```

Here is an example workflow script,

```python
from funsies import Fun, reduce, shell
with Fun():
    # you can run shell commands
    cmd = shell('sleep 2; echo 👋 🪐')
    # and python ones
    python = reduce(sum, [3, 2])
    # outputs are saved at hash addresses
    print(f"my outputs are saved to {cmd.stdout.hash[:5]} and {python.hash[:5]}")
```

The workflow is just python, and is run using the python interpreter,

```bash
$ python hello-world.py
my outputs are saved to 4138b and 80aa3
```

The `Fun()` context manager takes care of connecting to the database. The
script should execute immediately; no work is done just yet because workflows
are lazily executed.

To execute the workflow, we trigger using the hashes above using the CLI,

```bash
$ funsies execute 4138b 80aa3
```

Once the workers are finished, results can be printed directly to stdout using
their hashes,

```bash
$ funsies cat 4138b
👋 🪐
$ funsies cat 80aa3
5
```

They can also be accessed from within python, from other steps in the
workflows etc. Shutting down the database and workers can also be performed
using the CLI,

```bash
$ funsies shutdown --all
```

## How does it work?

The design of **funsies** is inspired by
[git](https://git-scm.com/book/en/v2/Git-Internals-Git-Objects) and
[ccache](https://ccache.dev/). All files and variable values are abstracted
into a provenance-tracking DAG structure. Basically, "files" are identified
entirely based on what operations lead to their creation. This (somewhat
opinionated) design produces interesting properties that are not common in
workflow engines:

#### Incremental computation

funsies automatically and transparently saves all input and output "files".
This produces automatic and transparent checkpointing and [incremental
computing](https://en.wikipedia.org/wiki/Incremental_computing). Re-running
the same funsies script, **even on a different machine**, will not perform any
computations (beyond database lookups). Modifying the script and re-running it
will only recompute changed results. 

In contrast with e.g. Make, this is not based on modification date but
directly on the data history, which is more robust to changes in the workflow.

#### Decentralized workflows

Workflows and their elements are not identified based on any global indexing
scheme. This makes it possible to generate workflows fully dynamically from
any connected computer node, to merge or compose DAGs from different databases
and to dynamically re-parametrize them, etc.

#### No local file operations

All "files" are encoded in a redis instance or to a data directory, with no
local filesystem management required. funsies workers can even operate without
any permanent data storage, as is often the case in file-driven workflows
using only a container's [tmpfs](https://docs.docker.com/storage/tmpfs/).

## Recovering from failures

Raised exceptions in python codes, worker failures, missing output files and
other error conditions are automatically caught by funsies workers, providing
fault tolerance to workflows. Errors are logged on `stderr` with full
traceback and can be recovered from the database.

Steps that depend on failed ones propagate those
errors and their provenance. Errors can then be dealt with wherever it is most
appropriate to do so [using techniques from functional
programming.](https://fsharpforfunandprofit.com/rop/) 

As an example, consider a workflow that first runs a CLI program `simulate`
that ought to produce a `results.csv` file, which is subsequently analyzed
using a python function `analyze_data()`,

```python
import funsies as f

sim = f.shell("simulate data.inp", inp={"data.inp":"some input"}, out=["results.csv"])
final = f.reduce(analyze_data, sim.out["results.csv"])
```

In a normal python program, `analyze_data()` would need to guard against the
possibility that `results.csv` is absent, or risk a fatal exception. In the
above funsies script, if `results.csv` is not produced, then it is replaced by
an instance of `Error` which tracks the failing step. The workflow engine
automatically shortcircuit the execution of `analyze_data` and insteads
forward the `Error` to `final`. In this way, the value of `final` provides
direct error tracing to the failed step. Furthermore, it means that
`analyze_data` does not need it's own error handling code if its output is
optional or if the error is better dealt with in a later step.

This error-handling approach is heavily influenced by the `Result<T,E>` type
from [the Rust programming language.](https://doc.rust-lang.org/std/result/)


## Is it production-ready?

🧪 warning: funsies is research-grade code ! 🧪

At this time, the funsies API is fairly stable. However, users should know
that database dumps are not yet fully forward- or backward-compatible, and
breaking changes are likely to be introduced on new releases.

## Related projects
funsies is intended as a lightweight alternative to industrial workflow
engines, such as [Apache Airflow](https://airflow.apache.org/) or
[Luigi](https://github.com/spotify/luigi). We rely heavily on awesome python
libraries: [RQ library](https://github.com/rq/rq),
[loguru](https://github.com/Delgan/loguru),
[Click](https://click.palletsprojects.com/) and
[chevron](https://github.com/noahmorrison/chevron). We are inspired by
[git](https://git-scm.com/book/en/v2/Git-Internals-Git-Objects),
[ccache](https://ccache.dev/),
[snakemake](https://snakemake.readthedocs.io/en/stable/)
[targets](https://github.com/ropensci/targets),
[rain](https://github.com/substantic/rain) and others. A comprehensive list of
other worfklow engine can be found
[here.](https://github.com/pditommaso/awesome-pipeline)


## License

funsies is provided under the MIT license.

## Contributing

All contributions are welcome! Consult [the CONTRIBUTING](./CONTRIBUTING.md)
file for help. Please file issues for any bugs and documentation problems.
# How to contribute to funsies

funsies is a free software project, and we welcome every kind of contributions
(documentations, bug reports, fixes, etc.)

If you encounter an issue, please fill a report [using the github Issues
page.](https://github.com/aspuru-guzik-group/funsies/issues)

If you have a fix, use the [pull request functionality on
Github.](https://github.com/aspuru-guzik-group/funsies/pulls) Make sure your
commits are given concise, explanatory names and limit your file edits to the
minimal relevant parts. Always run formatters before comitting (as described
below).

## Instructions for contributors

This repo has a CI workflow that will run the test suite, mypy, linting etc.,
and all tests and lints should pass. If you want to run CI locally before
pushing, all of it is [automated using Nox for local
development.](https://nox.thea.codes/en/stable/) To install nox in a python
env, use

```bash
pip install nox
```

To automatically format code so that it passes linting, you can use

```bash
nox -rs fmt
```

To lint the code and run the mypy type checker, use

```bash
nox -rs lint
nox -rs mypy
```

Tests can be run using `nox -rs tests` which will run tests for three python
versions (if they are installed). Some of the tests are rather slow and so are
skipped by default, you can run the whole set using `nox -rs tests -- --cov`.
You'll need a redis server installed and available on `$PATH`. 

Note that nox will take care of installing packages from PyPI for each of the
above steps, so you shouldn't need to do anything besides installing nox
itself to get a working dev environment.

funsies is formatted using [black](https://github.com/psf/black) and
[isort](https://pypi.org/project/isort/), which are run automatically using
nox, as described above.

Importantly, **funsies is a statically typed program,** and every function
boundaries needs annotations. 

Note also that the code tries to maintain (as much as possible) a functionally
pure paradigm for operations encoded in the database. This means that adding
new functionality can be somewhat tricky. If you find yourself meddling with
the internals described below, please talk to me (@clavigne) first!


## Internals: how to read this code

funsies is not a particularly complex codebase (yet) and should be fairly
readable. Here is a basic summary of the internal architecture.

The low-level workflow generation and encoding can be found in
[_graph.py](https://github.com/aspuru-guzik-group/funsies/blob/master/src/funsies/_graph.py)
 and
[_funsies.py](https://github.com/aspuru-guzik-group/funsies/blob/master/src/funsies/_funsies.py). 
The former contains hash calculation, getters and setters for data and
dependencies, while the latter contains data structure for encoding
operations. Data is serialized in
[_serdes.py](https://github.com/aspuru-guzik-group/funsies/blob/master/src/funsies/_serdes.py)

At the lowest level, workflow execution uses runner functions registered in
[_run.py](https://github.com/aspuru-guzik-group/funsies/blob/master/src/funsies/_run.py)
and called from `_run.run_op()`. The correct ordering of executions (as well
as other graph traversal functions) are computed in
[_dag.py](https://github.com/aspuru-guzik-group/funsies/blob/master/src/funsies/_dag.py).

User functions are found in
[ui.py](https://github.com/aspuru-guzik-group/funsies/blob/master/src/funsies/ui.py),
[fp.py](https://github.com/aspuru-guzik-group/funsies/blob/master/src/funsies/fp.py)
and others for shell operations, python operations etc. Console entry points
(the `funsies` cli command) are grouped in
[_cli.py](https://github.com/aspuru-guzik-group/funsies/blob/master/src/funsies/_cli.py).
 

---
title: 'funsies: A minimalist, distributed and dynamic workflow engine'
tags:
  - workflow
  - Python
  - redis
  - decentralized
  - computational chemistry
authors:
  - name: Cyrille Lavigne^[Corresponding author.]
    orcid: 0000-0003-2778-1866
    affiliation: 1
  - name: Alán Aspuru-Guzik
    orcid: 0000-0002-8277-4434
    affiliation: "1, 2, 3, 4"
affiliations:
  - name: Department of Computer Science, University of Toronto, 40 St. George St, Toronto, Ontario M5S 2E4, Canada
    index: 1
  - name: Chemical Physics Theory Group, Department of Chemistry, University of Toronto, 80 St. George St, Toronto, Ontario M5S 3H6, Canada
    index: 2
  - name: Vector Institute for Artificial Intelligence, 661 University Ave Suite 710, Toronto, Ontario M5G 1M1, Canada
    index: 3
  - name: Lebovic Fellow, Canadian Institute for Advanced Research (CIFAR), 661 University Ave, Toronto, Ontario M5G, Canada
    index: 4
date: 27 April 2021
bibliography: paper.bib
---

# Summary

Large-scale, high-throughput computational investigations are increasingly
common in chemistry and physics. Until recently, computational chemistry was
primarily performed using all-in-one monolithic software
packages [@smith:2020; @aquilante:2020; @kuhne:2020;
@apra:2020; @barca:2020; @romero:2020]. However, the
limits of individual programs become evident when tackling complex
multifaceted problems. As such, it is increasingly common to use multiple
disparate software packages in a single computational pipeline, 
 often stitched together using shell scripts in
languages such as Bash, or using Python and other interpreted languages.

These complex computational pipelines are difficult to scale and automate, as
they often include manual steps and significant “human-in-the-loop” tuning.
Shell scripting errors are often undetected, which can compromise
scientific results. Conversely, exception-based error handling, the standard
approach in Python, can readily bring a computational workflow to a halt when
exceptions are not properly caught [@weimer:2008].

`funsies` is a set of Python programs and modules to describe, execute and
analyze computational workflows, with first-class support for shell scripting.
It includes a lightweight, decentralized workflow engine backed by a NoSQL
store. Using `funsies`, external programs and Python-based computations
are easily mixed together. Errors are detected and propagated throughout
computations. Automatic, transparent incremental computing (based on a hash
tree data structure) provides a convenient environment for iterative
prototyping of computationally expensive workflows.

# Statement of need


Modern workflow management programs used in the private sector, such as Apache
Airflow and Uber's Cadence, are robust and extremely scalable, but are
difficult to deploy. Scientific workflow management systems, many of which are
compiled in [@awesome_pipeline] and systematically reviewed in
[@molder_sustainable_2021], are easier to set up on high-performance computing
clusters, but are tuned to the needs of specific disciplines, such as
bioinformatics or machine learning. This includes, for example, the use of
configuration file formats (YAML, JSON, etc.), packaging tools (for example,
conda or Docker), locked-in compute providers (Amazon Web Services, Google
Cloud) and storage formats that may be common in specific scientific fields
but not throughout the greater community.

For our own group's research program, we wanted to have available a
lightweight workflow management system that could be readily deployed to new
and varied computational facilities and local workstations with minimal
effort. This system had to support our existing shell-based and Python-based
scripts, and be flexible enough for rapid prototyping all the way to
large-scale computational campaigns, and provide an embeddable solution that
can be bundled within other software [@lavigne_automatic_2020]. Finally, we
were looking for a tool that could integrate data generation and storage, to
avoid the common practice of transforming the filesystem into what is
effectively a schema-less database. We developed `funsies` to address those
needs.


# Features and Implementation

`funsies` is a Python library and a set of associated command-line tools.
Using the `funsies` library, general computational workflows are described in
lazily evaluated Python code. Operations in `funsies` are taken to be pure,
that is, all operation outputs are entirely and solely determined by their
inputs. Workflows are orchestrated using Python by manipulating pointers to
yet-to-be-calculated data. Workflow instructions are transparently translated
and saved as graph elements in a Redis database.

Computational evaluation is initiated by the user asking for specific output
value. The task graph from these final outputs is walked back all the way to
those operations with no dependencies. These initial operations are then
queued for execution. Lightweight worker processes, instantiated from the
command line on local or remote machines, connect to the Redis database and
start executing the workflow. For each operation, the worker checks if outputs
are already cached, and if not, executes the associated function and saves its
outputs. It then enqueues any dependents for execution, by itself or by other
workers. In this way, the entire computational graph is evaluated in a
distributed, decentralized fashion without any scheduler or manager program.
Errors in workflows are handled using a functional approach inspired by Rust
[@klabnik_rust_2019]. Specifically, exceptions are propagated through workflow
steps, canceling dependent tasks, without interrupting valid workflow
branches. This provides both easy error tracing and a high degree of fault
tolerance.


The main distinguishing feature of `funsies` is the hash tree structure that
is used to encode all operations and their inputs. The causal hashing approach
used in `funsies` can also be found in Snakemake [@molder_sustainable_2021] as
an optional component and the (now defunct) Koji workflow
system [@maymounkov_koji_2018], as part of the Nix package
manager [@dolstra_nix_2004] and in the Git version control
system [@chacon_pro_2014]. In `funsies`, we replace all filesystem operations
with hash addressed operations; that is all I/O operations and dependencies
are tracked.

Every operation has a hash address that is computed from the hash values of
its dependencies and a hashed identifier for the associated operation on data.
In this way, the consistency of data dependencies is strongly enforced.
Changes to data and operations are automatically and transparently propagated,
as changing a single dependency will cause a rehash of all its dependents,
effectively producing a new workflow with no associated data that needs to be
recomputed. Alternatively, if data already exists at a specific hash address,
then it was generated from the same operations that produced that hash. In
this way, the hash tree structure enables transparent and automatic
incremental recomputing. 

Using hash addresses also enables decentralization, as we can rely on the
unlikeliness of hash collisions [@stevens_first_2017] to eliminate centralized
locks. An important advantage of this approach is that it allows worker
processes to generate their own workflows of tasks dynamically. Results from
these dynamic workflows can be collected and used further in the workflow
description, provided they can be reduced to a number of outputs known at
compile time, a technique similar to MapReduce [@dean_mapreduce_2004].

As of now, we have published one project [@pollice:2021] that used an earlier
iteration of `funsies`, and are using it in multiple ongoing inquiries. We
provide several sample workflows on GitHub, with a focus on computational
chemistry, quantum computing, and high-performance computing infrastructure.

We intend to maintain `funsies` and of course welcome [collaborations from
contributors around the world.](https://github.com/aspuru-guzik-group/funsies/blob/master/CONTRIBUTING.md)


# Acknowledgements

We acknowledge testing by early users Cher-Tian Ser, Kjell
Jorner and Gabriel dos Passos Gomes. CL also
thanks Chris Crebolder for help setting up documentation pages.
We acknowledge the Defense Advanced Research Projects Agency (DARPA) under the
Accelerated Molecular Discovery Program under Cooperative Agreement No.
HR00111920027 dated August 1, 2019. The content of the information presented
in this work does not necessarily reflect the position or the policy of the
Government. A. A.-G. thanks Dr. Anders G. Frøseth for his generous support. A.
A.-G. also acknowledges the generous support of Natural Resources Canada and
the Canada 150 Research Chairs program. We thank Compute Canada for providing
computational resources.

# References
In `funsies`, workflows written in pure python are saved to
[redis](https://redis.io/), a distributed, in-memory data store. Workflow
execution is performed using the minimal distributed queuing library
[RQ](https://python-rq.org/). Workflows are automatically parallelized and
computed incrementally. Command-line tools are provided that allow funsies to
be easily integrated in pre-existing shell workflows.

Incremental computation and caching are generated by hashing the steps
required to generate all file objects. Workflows are encoded using a [Merkle
tree](https://en.wikipedia.org/wiki/Merkle_tree) data structure, similar to
what is used by [ccache](https://ccache.dev/manual/4.2.html#_how_ccache_works)
for incremental computation or by distributed version control systems such as [mercurial](https://ericsink.com/vcbe/html/repository_structure.html).

Workflows are written in pure python using a set of primitives such as shell
commands (using the `shell()` function) and python glue code (encoded with
`py()`). Errors are handled using a functional programming (a `errors.Result`
monad) inspired by [Rust](https://doc.rust-lang.org/std/result/). In practice,
python code can raise exceptions and shell commands can fail in specifc
branches without compromising the execution of the whole workflow, and the
origin of errors is readily traced.



# Funsies by example

Here we show some example funsies script:

#### [Using funsies with SLURM for computational chemistry](slurm-conformers/README.md)
#### [Implementing mergesort as a funsies dynamic workflow](silly/README.md)
## funsies + SLURM = 💖
A major reason why funsies was created is to orchestrate computationally
expensive complicated workflows at high-performance computing centers. In
effect, funsies tries to fill the niche between full-blown workflow engine
like [airflow](https://airflow.apache.org/) and hastily stitched together
pieces of bash script, [GNU parallel](https://www.gnu.org/software/parallel/)
and FORTRAN 77.

In this recipe, I describe how to integrate funsies with
[slurm](https://slurm.schedmd.com/documentation.html) the de facto standard
resource manager for scientific computing clusters.

### Problem statement
Molecules are usually represented as simple stick diagrams. Obviously, they
don't move while on paper, which gives the impression that they are static.
That's not the case in real life.

Molecules are dynamic entities and they freely flop around in vacuum. The
stable geometries of a molecule (its potential energy minima) that can rapidly
interconvert at room temperature are called conformers. The floppier the
molecule, the more conformers it has. 

In this problem, we will use funsies to compute conformer energies for a
floppy molecule, and verify that the conformers can interconvert at room
temperature.

### Computational workflow
For our computational workflow, we start from the SMILES representation of an
alkane diol. Then,
1. We use [OpenBabel](http://openbabel.org/wiki/Main_Page) to systematically
   find its conformers.
2. We optimize each of those conformers with
   [xtb](https://github.com/grimme-lab/xtb) and sort them.
3. We output the result as a JSON file.

The entire workflow is in [workflow.py](./workflow.py). This workflow contains
dynamic workflow generation (to account for the fact that the number of
conformers is initially not known), shell commands and some python functions.

### Deploy and execute
To deploy, we first setup a conda environment,
```bash
conda create -n funsies
conda install -c conda-forge xtb openbabel
conda install redis-server
pip install funsies
```
This will create an environment with xtb, openbabel, funsies and their
dependencies. Now that all this is setup, we submit using [a standard SLURM submission
script](./slurm-submit.sh),
```bash
sbatch slurm-submit.sh
```
The header of the submission script should be modified to match that of your
SLURM account.

What does the submission script do? It instantiate a redis server, starts a
number of workers on individual compute nodes, tells worker processes how to
find the server, runs the python workflow then shutdown workers and server.

All this in only 9 lines of shell!

The final step of the job dumps an image of the redis database in
`results.rdb`. This image includes all temporary results. To run a "mock"
version of the computation, simply rename or copy the file to `dump.rdb` and
start a server. For example, doing
```bash
cp results.rdb dump.rdb
redis-server &
python workflow.py
```
will run through the workflow again locally, but without recomputing anything.
To look at the computational graph using graphviz one can simply do,
```bash
cp results.rdb dump.rdb
redis-server &
funsies graph > graph.dot
dot -Tpdf graph.dot > graph.pdf
```
This will generate [the graph shown here.](./graph.pdf) 
## Mergesort

Here, we have a [funsies-based implementation](./mergesort.py) of the
mergesort algorithm, using recursion. This is quite possibly the least
efficient way to sort a list of integers (parallel though! 😁), but it does
demonstrate quite effectively dynamic DAG generation. The attached script will
sort a random list of 120 integers, which requires 7 nested workflows, all
generated recursively and dynamically. The final graph [is rather
interesting.](./graph.pdf)

The main challenge in this toy example is that we have to terminate our
recursions without explicitly extracting the funsies data (for euhm
performance reasons?). We use error propagation for this: basically in each
nested sub-workflow, we raise the `StopRecursion` exception if less than 2
elements are present and stop recursing deeper. What we get is basically
recursive, nested
[MapReduce.](https://hadoop.apache.org/docs/current/hadoop-mapreduce-client/hadoop-mapreduce-client-core/MapReduceTutorial.html)

Although this is rather silly, a similar approach could conceivably be used
for large-scale search problems using [a divide-and-conquer
algorithm](https://en.wikipedia.org/wiki/Bisection_method). It also
demonstrates how funsies can be used for pure python problems.



