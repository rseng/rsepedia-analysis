<p align="center">
  <a href="https://github.com/fAndreuzzi/BisPy" target="_blank" >
    <img alt="BisPy" src="logo.png" width="400" />
  </a>
</p>

![Python package](https://github.com/fAndreuzzi/BisPy/workflows/Python%20package/badge.svg?branch=master)
<a href='https://coveralls.io/github/fAndreuzzi/BisPy'><img src='https://coveralls.io/repos/github/fAndreuzzi/BisPy/badge.svg' alt='Coverage Status' /></a>
[![GitHub license](https://img.shields.io/github/license/Naereen/StrapDown.js.svg)](https://github.com/Naereen/StrapDown.js/blob/master/LICENSE)
<img src='https://img.shields.io/badge/Code%20style-Black-%23000000'/>
[![Documentation Status](https://readthedocs.org/projects/bispy-bisimulation-in-python/badge/?version=latest)](https://bispy-bisimulation-in-python.readthedocs.io/en/latest/?badge=latest)
[![status](https://joss.theoj.org/papers/9d9c3ca0715d482938b5a450525cefa0/status.svg)](https://joss.theoj.org/papers/9d9c3ca0715d482938b5a450525cefa0)
[![PyPI version](https://badge.fury.io/py/BisPy.svg)](https://badge.fury.io/py/BisPy)

## Description

**BisPy** is a Python package for the computation of the maximum bisimulation
of directed graphs. At the moment it supports the following algorithms:

- Paige-Tarjan
- Dovier-Piazza-Policriti
- Saha

A brief introduction to the problem can be found
[here](https://bispy-bisimulation-in-python.readthedocs.io/en/latest/?badge=latest#a-brief-introduction-to-bisimulation).

## Usage

### Paige-Tarjan, Dovier-Piazza-Policriti

To compute the maximum bisimulation of a graph, first of all we import
`paige_tarjan` and `dovier_piazza_policriti` from **BisPy**, as well as the
package _NetworkX_, which we use to represent graphs:

```python
>>> import networkx as nx
>>> from bispy import compute_maximum_bisimulation, Algorithms
```

We then create a simple graph:

```python
>>> graph = nx.balanced_tree(2,3, create_using=nx.DiGraph)
```

It's important to set `create_using=nx.DiGraph` since **BisPy** works only with
_directed_ graphs. Now we can compute the _maximum bisimulation_ using
_Paige-Tarjan_'s algorithm, which is the default for the function
`compute_maximum_bisimulation`:

```python
>>> compute_maximum_bisimulation(graph)
[(7, 8, 9, 10, 11, 12, 13, 14), (3, 4, 5, 6), (1, 2), (0,)]
```

We can use _Dovier-Piazza-Policriti_'s algorithm as well:

```python
>>> compute_maximum_bisimulation(graph, algorithm=Algorithms.DovierPiazzaPolicriti)
[(7, 8, 9, 10, 11, 12, 13, 14), (3, 4, 5, 6),  (1, 2), (0,)]

```

We may also introduce a _labeling set_ (or _initial partition_):

```python
>>> compute_maximum_bisimulation(graph, initial_partition=[(0,7,10), (1,2,3,4,5,6,8,9,11,12,13,14)])
[(7, 10), (5, 6), (8, 9, 11, 12, 13, 14), (3, 4), (2,), (1,), (0,)]

```

### Saha

In order to use _Saha_'s algorithm we only need to import the following
function:

```python
>>> from bispy import saha
```

We call that function to obtain an object of type `SahaPartition`, which has a
method called `add_edge`. This method adds a new edge to the graph and
recomputes the maximum bisimulation incrementally:

```python
saha_partition = saha(graph)
```

(We reused the `graph` object which we defined in the previous paragraph). We
can now use the aforementioned method `add_edge` (note that when you call this
method the instance of `graph` which you passed is **not** modified):

```python
>>> for edge in [(1,0), (4,0)]:
...    maximum_bisimulation = saha_partition.add_edge(edge)
...    print(maximum_bisimulation)
[(3, 4, 5, 6), (7, 8, 9, 10, 11, 12, 13, 14), (0,), (2,), (1,)]
[(3, 5, 6), (7, 8, 9, 10, 11, 12, 13, 14), (0,), (2,), (1,), (4,)]
```

## Documentation

You can read the documentation (hosted on ReadTheDocs) at this
[link](https://bispy-bisimulation-in-python.readthedocs.io/en/latest/?badge=latest).

To build the HTML version of the docs locally use:

```bash
> cd docs
> make html
```

The generated html can be found in `docs/build/html`.

## Dependencies and installation

**BisPy** requires the modules `llist, networkx`. The code is tested
for _Python 3_, while compatibility with _Python 2_ is not guaranteed. It can
be installed using `pip` or directly from the source code.

### Installing via _pip_

To install the package:

```bash
> pip install bispy
```

To uninstall the package:

```bash
> pip uninstall bispy
```

### Installing from source

You can clone this repository on your local machine using:

```bash
> git clone https://github.com/fAndreuzzi/BisPy
```

To install the package:

```bash
> cd BisPy
> python setup.py install
```


## Testing

We are using **GitHub actions** for continuous intergration testing. To run
tests locally (`pytest` is required) use the following command from the root
folder of **BisPy**:

```bash
> pytest tests
```

## Authors and acknowledgements

**BisPy** is currently developed and mantained by **Francesco Andreuzzi**. You
can contact me at:

- andreuzzi.francesco at gmail.com
- fandreuz at sissa.it

The project has been developed under the supervision of professor **Alberto
Casagrande** (_University of Trieste_), which was my advisor for my _bachelor
thesis_.

## Reporting a bug

The best way to report a bug is using the
[Issues](https://github.com/fAndreuzzi/BisPy/issues) section. Please, be clear,
and give detailed examples on how to reproduce the bug (the best option would
be the graph which triggered the error you are reporting).

## How to contribute

We are more than happy to receive contributions on tests, documentation and
new features. Our [Issues](https://github.com/fAndreuzzi/BisPy/issues)
section is always full of things to do.

Here are the guidelines to submit a patch:

1. Start by opening a new [issue](https://github.com/fAndreuzzi/BisPy/issues)
   describing the bug you want to fix, or the feature you want to introduce.
   This lets us keep track of what is being done at the moment, and possibly
   avoid writing different solutions for the same problem.

2. Fork the project, and setup a **new** branch to work in (_fix-issue-22_, for
   instance). If you do not separate your work in different branches you may
   have a bad time when trying to push a pull request to fix a particular
   issue.

3. Run [black](https://github.com/psf/black) before pushing
   your code for review.

4. Any significant changes should almost always be accompanied by tests. The
   project already has good test coverage, so look at some of the existing
   tests if you're unsure how to go about it.

5. Provide menaningful **commit messages** to help us keeping a good _git_
   history.

6. Finally you can submbit your _pull request_!

## References

We consulted the following resources during the development of **BisPy**:

- Saha, Diptikalyan. "An incremental bisimulation algorithm." International
  Conference on Foundations of Software Technology and Theoretical Computer
  Science. Springer, Berlin, Heidelberg, 2007.
  [DOI](https://doi.org/10.1007/978-3-540-77050-3_17)
- Dovier, Agostino, Carla Piazza, and Alberto Policriti. "A fast bisimulation
  algorithm." International Conference on Computer Aided Verification.
  Springer, Berlin, Heidelberg, 2001.
  [DOI](https://doi.org/10.1007/3-540-44585-4_8)
- Gentilini, Raffaella, Carla Piazza, and Alberto Policriti. "From bisimulation
  to simulation: Coarsest partition problems." Journal of Automated Reasoning
  31.1 (2003): 73-103. [DOI](https://doi.org/10.1023/A:1027328830731)
- Paige, Robert, and Robert E. Tarjan. "Three partition refinement algorithms."
  SIAM Journal on Computing 16.6 (1987): 973-989.
  [DOI](https://doi.org/10.1137/0216062)
- Hopcroft, John. "An n log n algorithm for minimizing states in a finite
  automaton." Theory of machines and computations. Academic Press, 1971.
  189-196.
- Aczel, Peter. "Non-well-founded sets." (1988).
- Kanellakis, Paris C., and Scott A. Smolka. "CCS expressions, finite state
  processes, and three problems of equivalence." Information and computation
  86.1 (1990): 43-68. [DOI](<https://doi.org/10.1016/0890-5401(90)90025-D>)
- Sharir, Micha. "A strong-connectivity algorithm and its applications in data
  flow analysis." Computers & Mathematics with Applications 7.1 (1981): 67-72.
  [DOI](<https://doi.org/10.1016/0898-1221(81)90008-0>)
- Cormen, Thomas H., et al. Introduction to algorithms. MIT press, 2009.
  (ISBN: 9780262533058)

## License

See the [LICENSE](LICENSE) file for license rights and limitations (MIT).
---
title: "BisPy: Bisimulation in Python"
tags:
  - Bisimulation
  - Graph theory
  - Graph algorithms
authors:
  - name: Francesco Andreuzzi
    orcid: 0000-0002-9508-7801
    affiliation: "1,2"
affiliations:
  - name: Internation School of Advanced Studies, SISSA, Trieste, Italy
    index: 1
  - name: Università degli Studi di Trieste
    index: 2
date: 11 July 2021
bibliography: paper.bib
---

# Summary

A binary relation $\mathcal{B}$ on the set $V$ of the nodes of a directed graph
is a bisimulation if the following condition is satisfied [@gentilini]:

\begin{gather} (a,b) \in \mathcal{B} \implies \begin{cases} a \to a' &\implies
\exists b' \in V \mid (a',b') \in \mathcal{B} \land b \to b'\\ b \to b'
&\implies \exists a' \in V \mid (a',b') \in \mathcal{B} \land a \to a'
\end{cases} \end{gather}

A _labeling function_ $\ell : V \to L$ may be introduced, in which case the
graph becomes a _Kripke structure_ and the additional condition
$(a,b) \in \mathcal{B} \implies \ell(a) = \ell(b)$ must be satisfied.

<p style="text-align: center;">

![On the left, a balanced tree paired with a labeling function, which induces a partition on $V$ of cardinality 2. We visually represent the corresponding maximum bisimulation on the right, computed using \texttt{BisPy}.\label{fig:example}](example.png)

</p>

The notion of _bisimulation_ and in particular of _maximum bisimulation_ —
namely the bisimulation which contains all the other bisimulations on the graph
— has applications in modal logic, formal verification, and concurrency theory
[@kanellakis], and is used for graph reduction as well [@gentilini]. The fact
that _graphs_ may be used to create digital models of a wide span of complex
systems makes bisimulation a useful tool in many different cases. For this
reason several algorithms for the computation of maximum bisimulation have been
studied throughout the years, and it is now known that the problem has an
$O(|E| \log |V|)$ algorithmic solution [@paigetarjan], where $V$ is the set of
nodes in the graph, and $E$ is the set of edges of the graph.

$\texttt{BisPy}$ is a Python package for the computation of maximum
bisimulation.

# Statement of need

To the best of our knowledge, \texttt{BisPy} is the first Python project to
address the problem presented above, and to meet the objectives of healthy
open source software, namely extensive testing, documentation, and intuitive
code commenting.

We think that our project may be a useful tool to study practical cases for
students approaching the field — since the notion of bisimulation may be
somewhat counterintuitive at first glance — as well as for established
researchers, who may use \texttt{BisPy} to study improvements on particular
types of graphs and to compare new algorithms with the state of the art.

It is interesting to observe that the package \texttt{BisPy}, briefly presented
below, contains the implementation of more than one algorithm for the
computation of maximum bisimulation, and every algorithm uses a peculiar
strategy to obtain the result. For this reason, we think that our package may be
useful to assess the performance of different approaches on a particular
problem.

# \texttt{BisPy}

Our package contains the implementation of the following algorithms:

- Paige-Tarjan [-@paigetarjan], which employs an insight from the famous
  algorithm for the minimization of finite states automata [@hopcroft];
- Dovier-Piazza-Policriti [-@dovier], which uses the notion of _rank_ to
  optimize the overhead of splitting the initial partition, and can be computed
  — prior the execution of the algorithm — using an $O(|V|+|E|)$ procedure
  [@sharir; @tarjan];
- Saha [-@saha], which can be used to update the maximum bisimulation of a
  graph after the addition of a new edge, and is more efficient than the
  computation _from scratch_ in some cases (the computational complexity
  depends on how much the maximum bisimulation changes due to the
  modification).

Our implementations have been tested and documented deeply; moreover we split
the algorithms into smaller functions, which we prefer to having a
monolithic block of code in order to improve readability and testability. This
kind of modularity allows us to reuse functions across multiple algorithms,
since several procedures are shared (e.g., $\texttt{split}$ is used in all
three of the algorithms that we mentioned above, while the computation of rank is
carried out only in the last two), and for the same reason we think that the
addition of new functionalities would be straightforward since we have already
implemented a significant set of common functions.

# Example

We present the code that we used to generate the example shown in
\autoref{fig:example}. First of all we import the modules needed to generate
the graph (\texttt{BisPy} takes \texttt{NetworkX} directed graphs in input) and
to compute the maximum bisimulation.

```python
>>> import networkx as nx
>>> from bispy import compute_maximum_bisimulation
```

After that we generate the graph, which as we mentioned before is a balanced
tree with _branching-factor_=2 and _depth_=3. We also create a list of tuples
that represents the labeling function which we employed in the example.

```python
>>> graph = nx.balanced_tree(2,3, create_using=nx.DiGraph)
>>> labels = [(0,1,2,3,4,5,6,7,9,10,11,12,13),(8,14)]
```

We can now compute the maximum bisimulation of the Kripke structure taken into
account as follows:

```python
>>> compute_maximum_bisimulation(graph, labels)
[(4,5),(7,9,10,11,12,13),(8,14),(3,6),(0,),(1,2)]
```

The visualization shown above has been drawn using the library
\texttt{PyGraphviz}. \texttt{BisPy} provides the requested output in the form
of a list of tuples, each of which contains the labels of all the nodes that
are members of an equivalence class of the maximum bisimulation.

# Performance

We briefly examine some performance results on two different kinds of graphs:

- _Balanced trees_ [@clrs] with variable branching factor $r$ and height $h$,
  for which we are going to use the notation $B_T(r,h)$;
- _Erdős-Rényi graphs_ [-@clrs], also called _binomial graphs_, whose set $E$
  of edges is generated randomly (the cardinality $|E|$ is roughly $p|V|$).

The first experiment involves balanced trees, and consists of the computation
of the maximum bisimulation of trees with variable dimensions. The labeling set
is the trivial partition of the set $V$. The results are shown in the left side
of \autoref{fig:performance}. The quantity that varies along the x-axis is
$|E| \log |V|$, since this allows the presentation of data in a more natural
way.

The performance complies with the expected complexity $|E| \log |V|$: for
instance our implementation of Dovier-Piazza-Policriti takes about 1.425
seconds to compute the maximum bisimulation on $B_T(3,10)$, and 12.596 seconds
on $B_T(3,12)$. The value of the ratio
$\frac{|E_{B_T(3,12)}| \log |V_{B_T(3,12)}|}{|E_{B_T(3,10)}| \log |V_{B_T(3,10)}|}$
is approximately 10.7, therefore the growth of the time function respects
approximately the predicted behavior.

Concerning binomial graphs, we fixed $p=0.0005$ in order to obtain a graph of
some practical interest (as $p \to 1$ the graph becomes complete, as $p \to 0$
also $|E| \to 0$). This time we also consider Saha's incremental algorithm, and
we conduct the experiment as follows:

1. Generate a binomial graph with the aforementioned features;
2. Compute the maximum bisimulation using Paige-Tarjan's algorithm;
3. Add a random edge to the graph;
4. Compute the updated maximum bisimulation three times, using the three
   algorithms taken into account, and verify the time taken by each one.

Since the experiment is not deterministic (the graph and the new edge are
generated randomly) we evaluate and visualize the mean time taken by step 4 on
a sample of 1000 iterations of steps 1-4.

The knowledge of the old maximum bisimulation is of no interest for
non-incremental algorithms. However Saha's algorithm uses this input to reduce
the number of steps: the goal of the second experiment is in fact to illustrate
this improvement. The results are shown in the right side of
\autoref{fig:performance}.

<p style="text-align: center;">

![On the left side of the figure, the time taken by our implementations of Paige-Tarjan and Dovier-Piazza-Policriti to compute the maximum bisimulation of balanced trees with variable branching factor and height. On the right side, the time needed to update the maximum bisimulation of a binomial graph after the addition of a random edge (for this experiment we also consider Saha's incremental algorithm).\label{fig:performance}](performance.png)

</p>

We ran the experiments on a workstation with operating system _CentOS Linux_,
(x86*64), processor Intel(R) Core(TM) i7-4790 CPU (4 cores, 3.60GHz), and 16 GB
RAM. Graphs have been generated using functions from the Python package
\_NetworkX* [@networkx]. We measured time using the Python module _timeit_
[@pythondocs].

# Acknowledgements

We acknowledge the support received from Alberto Casagrande during the
preliminar theoretical study of the topic, as well as SISSA mathLab for
providing the hardware used to perform experiments on large graphs.

# References
