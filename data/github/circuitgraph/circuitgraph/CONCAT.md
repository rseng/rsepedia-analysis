# Changelog
All notable changes to this project will be documented in this file.

## [0.1.3] - 2021-09-24
### Added
- More robust checks for incorrect circuit construction
- More robust parsing, including faster parsing using regex
- Simple circuit visualization
- `X` type for nodes (similar to `0`, and `1`)
- Sequential unroll transform

## [0.1.2] - 2021-01-24
### FIXED
- Synthesis now works with python3.6 again

### Added
- DesignCompiler synthesis

## [0.1.1] - 2021-01-22
### FIXED
- Parsing is now being included correctly

### Added
- Lark based verilog parsing to vastly speed up reading verilog netlists

### Changed
- Replaced memory elements with a more general blackbox-based scheme

## [0.0.3] - 2020-09-08
### Fixed
- Image link in README is external so that it will appear in pypi

## [0.0.2] - 2020-09-08
### Fixed
- Constant valued nodes are now written to verilog
- README installation documentation to point to the pypi repo

### Removed
- pyverilator support, because it breaks installation on python distributions built without tkinter
- pyeda dependency because it is no longer needed for parsing

## [0.0.1] - 2020-08-30
Initial release
<img src="https://raw.githubusercontent.com/circuitgraph/circuitgraph/master/docs/circuitgraph.png" width="300">

# CircuitGraph

[![Build Status](https://app.travis-ci.com/circuitgraph/circuitgraph.svg?branch=master)](https://app.travis-ci.com/github/circuitgraph/circuitgraph)
[![codecov](https://codecov.io/gh/circuitgraph/circuitgraph/branch/master/graph/badge.svg)](https://codecov.io/gh/circuitgraph/circuitgraph)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/pre-commit/pre-commit)

[**CircuitGraph**](https://circuitgraph.github.io/circuitgraph/) is a library for working with hardware designs as graphs. CircuitGraph provides an interface to do this built on [NetworkX](https://networkx.github.io), along with integrations with other useful tools such as sat solvers and the [Yosys](http://www.clifford.at/yosys/) synthesis tool, and input/output to verilog.

## Overview

The `Circuit` class is at the core of the library and it is essentially a wrapper around a [NetworkX](https://networkx.github.io) graph object. This graph is accessable through the `graph` member variable of `Circuit` and can be used as an entrypoint to the robust NetworkX API.

Here's a simple example of reading in a verilog file, adding a node to the graph, and writing back to a new file.

```python
import circuitgraph as cg

c = cg.from_file('/path/to/circuit.v')
# Add an AND gate to the circuit that takes as input nets o0, o1, o2, o3
c.add('g', 'and', fanin=[f'o{i}' for i in range(4)])
cg.to_file(c, '/path/to/output/circuit.v')
```

The documentation can be found [here](https://circuitgraph.github.io/circuitgraph/).

## Installation

CircuitGraph requires Python3. 
The easiest way to install is via PyPi:
```shell
pip install circuitgraph
```
To install from the release, download and:
```shell
pip install circuitgraph-<release>.tar.gz
```

Finally, to install in-place with the source, use:
```shell
cd <install location>
git clone https://github.com/circuitgraph/circuitgraph.git
cd circuitgraph
pip install -r requirements.txt
pip install -e .
```
### Optional Packages

In addition to the packages enumerated in `requirements.txt`, there are a few tools you can install to enable additional functionality.

If you would like to use the satisfiability functionality, install [PySAT](https://pysathq.github.io).

If you would like to perform synthesis you can install either Cadence Genus or [Yosys](http://www.clifford.at/yosys/). If you're going to use Genus, you must provide the path to a synthesis library to use by setting the `CIRCUITGRAPH_GENUS_LIBRARY_PATH` variable. 

## Contributing

If you have ideas on how to improve this library we'd love to hear your suggestions. Please open an issue. 
If you want to develop the improvement yourself, please consider the information below.

CI Testing and coverage is setup using [Travis CI](https://travis-ci.org/) and [Codecov](https://codecov.io). 
 If you would like to generate coverage information locally, install coverage and codecov.
```shell
pip install coverage codecov 
make coverage
```

Documentation is built using pdoc3.
```shell
pip install pdoc3
make doc
```

Tests are run using the builtin unittest framework.
```shell
make test
```

Code should be formatted using [black](https://black.readthedocs.io/en/stable/). 
[Pre-commit](https://pre-commit.com) is used to automatically run black on commit. 
```shell
pip install black pre-commit
pre-commit install
```

## Citation

If you use this software for your research, we ask you cite this publication:
https://joss.theoj.org/papers/10.21105/joss.02646

```
@article{sweeney2020circuitgraph,
  title={CircuitGraph: A Python package for Boolean circuits},
  author={Sweeney, Joseph and Purdy, Ruben and Blanton, Ronald D and Pileggi, Lawrence},
  journal={Journal of Open Source Software},
  volume={5},
  number={56},
  pages={2646},
  year={2020}
}
```

## Acknowledgements

Circuitgraph icon designed by [ncasti](https://github.com/ncasti).
---
title: 'CircuitGraph: A Python package for Boolean circuits'
tags:
  - Python
  - Boolean circuits
  - satisfiability
  - graph
  - electronic design automation
authors:
- name: Joseph Sweeney
  affiliation: 1
- name: Ruben Purdy
  affiliation: 1
- name: Ronald D Blanton
  affiliation: 1
- name: Lawrence Pileggi
  affiliation: 1
affiliations:
- name: Department of Electrical and Computer Engineering, Carnegie Mellon University, Pittsburgh, PA 15213
  index: 1
date: 13 August 2020
bibliography: paper.bib

---

# Summary

A Boolean circuit is a fundamental mathematical model ubiquitous in the 
design of modern computers. The model consists of a directed graph wherein 
nodes are logic gates with corresponding Boolean functions and edges are wires 
which determine the composition of the gates. `CircuitGraph` is a open-source Python library
for manipulating and analyzing Boolean circuits. 

# Statement of need 

Analysis, manipulation, and generation of Boolean circuits is fundamental to many aspects of
digital hardware design, cryptography, constraint solving, and other areas. 
Highly optimized software for processing Boolean circuits exists. Unfortunately
it generally is proprietary, with expensive license fees. Furthermore, these
options suffer from poor documentation, are closed source, and typically 
rely on Tool control language (Tcl). While simple, Tcl is slow, has limited
libraries and supporting community, and is unnecessarily verbose. These reasons
motivate the development of our open source solution. While this software will 
directly benefit our lab as a research platform, it certainly has application 
in other environments such as the classroom. 

# Functionality

The functionality of `CircuitGraph` has been tailored to our research needs, however,
the library is easily extensible to many other applications of Boolean circuits. In the 
following sub sections, we highlight some of the library's key features.

The core of the library is the `Circuit` class, which internally uses a `networkx.DiGraph` 
data structure from @SciPyProceedings_11. The class implements key Boolean circuit functionalities 
on top of the graph as we describe below.

## Interfaces

Compatibility with existing systems is a primary goal for our library. Towards this end, 
we have built interfaces for a subset of Verilog, the most commonly used Boolean circuit format.
This library supports generic stuctural Verilog which is the typical output of synthesis tools. 
Specifically, the library can parse combinational gates in the following forms. We also provide an interface 
to parse sequential elements. 

```verilog
assign a = b|(c^d);
xor(e,f,g);
```
Additionally, we have provided a library of generic and benchmark circuits that can be quickly instantiated.

```python
import circuitgraph as cg
c0 = cg.from_file('path/circuit.v')
c1 = cg.from_file('path/circuit.bench')
c2 = cg.from_lib('c17')
```

## Composition

A common issue found in similar tools is the poor expressivity of circuit construction 
primitives. We aim to build a simple, but powerful syntax for creating and connecting nodes
in a circuit. The ease of our syntax is enabled by the Python language. 
An example of this syntax is below.

```python
# add an OR gate named 'a'
c0.add('a','or')

# create an AND gate with circuit inputs in a single line. Input connections to the gate can be spcified with the fanin argument, output connections with fanout. 
c0.add('g','and',fanin=[c.add(f'in_{i}','input') for i in range(4)])
```

## Synthesis
We provide an interface to common synthesis tools including `yosys` from @wolf2019yosys and `Cadence Genus`. This allows 
the user to run basic synthesis routines on circuits from within Python. Specifically, we support the generic multi-level synthesis routines of both tools. 
```python
# synthesize circuit with yosys
c_syn = cg.syn(c0, "Yosys")
```

## Satisfiability

Satisfiability is an essential problem related to Boolean circuits. Surprisingly, commercial 
synthesis tools do not directly support its use (although the open source tools yosys does). 
We add satisfiability to our library which in turn enables a wide array of analysis including
sensitization, sensitivity, and influence. Our implementation utilizes `pysat` from @imms-sat18. 
The main interface is simple allowing the user to determine
satisfiability of a circuit under a set of variable assignments. To develop more complex routines, the user can also access the underlying `pysat.solver` instance. 
In conjunction with satisfiability, we provide interfaces to approximate and exact model count algorithms. 

```python
# check satisfiability assuming 'a' is False
cg.sat(c0,{'a':False})

# get number of solutions to circuit with 'a' False
cg.model_count(c0,{'a':False})
```

# Future Work
We plan on adding support for the BLIF and Bench formats. Additionally, we may expand the compatibility with Verilog standards if a need is shown. Support for timing-based synthesis may be useful in some scenarios. Other improvements could interfaces to open source simulation and Automatic Test Pattern Generation (ATPG) tools. 


# Requirements

As previously mentioned, `CircuitGraph` relies on the `networkx` and `pysat` libraries. Additionally, it uses `pyverilog`
to parse verilog netlists. 

# References
