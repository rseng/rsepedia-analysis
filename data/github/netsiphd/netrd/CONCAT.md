[![DOI](https://joss.theoj.org/papers/10.21105/joss.02990/status.svg)](https://doi.org/10.21105/joss.02990)
[![PyPI version](https://badge.fury.io/py/netrd.svg)](https://badge.fury.io/py/netrd)
[![ReadTheDocs](https://img.shields.io/readthedocs/netrd.svg)](
    https://netrd.readthedocs.io)
![CI](https://github.com/netsiphd/netrd/workflows/build/badge.svg)

# netrd: A library for network {reconstruction, distances, dynamics}

This library provides a consistent, NetworkX-based interface to various
utilities for graph distances, graph reconstruction from time series data, and
simulated dynamics on networks. 

Some resources that maybe of interest:

* A [tutorial](https://netrd.readthedocs.io/en/latest/tutorial.html) on how to use the library
* The API [reference](https://netrd.readthedocs.io/en/latest/) 
* A [notebook](https://nbviewer.jupyter.org/github/netsiphd/netrd/blob/master/notebooks/example.ipynb) showing advanced usage

# Installation

`netrd` is easy to install through pip:

```
pip install netrd
```

If you are thinking about contributing to `netrd`, you can install a
development version by executing

```
git clone https://github.com/netsiphd/netrd
cd netrd
pip install .
```

# Usage

## Reconstructing a graph

<p align="center">
<img src="netrd_reconstruction_example.png" alt="example reconstruction" width="95%"/>
</p>

The basic usage of a graph reconstruction algorithm is as follows:

```python
from netrd.reconstruction import CorrelationMatrix
import numpy as np
# 100 nodes, 1000 observations
TS = np.random.random((100, 1000))

reconstructor = CorrelationMatrix()
G = reconstructor.fit(TS, threshold_type='degree', avg_k=15)
# or alternately, G = reconstructor.results['graph']
```

Here, `TS` is an N x L numpy array consisting of L
observations for each of N sensors. This constrains the graphs
to have integer-valued nodes.

The `results` dict object, in addition to containing the graph
object, may also contain objects created as a side effect of
reconstructing the network, which may be useful for debugging or
considering goodness of fit. What is returned will vary between
reconstruction algorithms.

Many reconstruction algorithms create a dense matrix of weights and
use additional parameters to describe how to create a sparse graph; the
[tutorial](https://netrd.readthedocs.io/en/latest/tutorial.html) has more
details on these parameters.


## Distances between graphs

<p align="center">
<img src="netrd_distance_example.png" alt="example distance" width="95%"/>
</p>

The basic usage of a distance algorithm is as follows:

```python
from netrd.distance import QuantumJSD
import networkx as nx
G1 = nx.fast_gnp_random_graph(1000, .1)
G2 = nx.fast_gnp_random_graph(1000, .1)

dist_obj = QuantumJSD()
distance = dist_obj.dist(G1, G2)
# or alternatively: distance = dist_obj.results['dist']
```

Here, `G1` and `G2` are `nx.Graph` objects (or subclasses such as
`nx.DiGraph`). The results dictionary holds the distance value, as
well as any other values that were computed as a side effect.

## Dynamics on graphs

<p align="center">
<img src="netrd_dynamics_example.png" alt="example distance" width="95%"/>
</p>

The basic usage of a dynamics algorithm is as follows:

```python
from netrd.dynamics import VoterModel
import networkx as nx
ground_truth = nx.karate_club_graph()

dynamics_model = VoterModel()
synthetic_TS = dynamics_model.simulate(ground_truth, 1000)
# this is the same structure as the input data to a reconstructor
# G = CorrelationMatrix().fit(synthetic_TS)
```

This produces a numpy array of time series data.


# Contributing

Contributing guidelines can be found in [CONTRIBUTING.md](CONTRIBUTING.md).


# Publications

* McCabe, S., Torres, L., LaRock, T., Haque, S. A., Yang, C.-H., Hartle, H., and
Klein, B. (2021). netrd: A library for network reconstruction and graph
distances. *Journal of Open Source Software* 6(62): 2990.
doi:&nbsp;[10.21105/joss.02990](https://doi.org/10.21105/joss.02990).
arXiv:&nbsp;[2010.16019](https://arxiv.org/abs/2010.16019).
    + paper detailing the methods used in this package

* Hartle H., Klein B., McCabe S., Daniels A., St-Onge G., Murphy C., and
Hébert-Dufresne L. (2020). Network comparison and the within-ensemble graph
distance. *Proceedings of the Royal Society A* 476: 20190744.
doi:&nbsp;[10.1098/rspa.2019.0744](http://doi.org/10.1098/rspa.2019.0744).
arXiv:&nbsp;[2008.02415](https://arxiv.org/abs/2008.02415).
    + recent work introducing a baseline measure for comparing graph distances
# netrd

Welcome to `netrd` and thanks for your interest in contributing! During development please
make sure to keep the following checklists handy. They contain a summary of all
the important steps you need to take to contribute to the package. As a general statement,
the more familiar you already are with git(hub), the less relevant the detailed instructions
below will be for you.

For an introduction and overview of the project, check out these
[slides](https://docs.google.com/presentation/d/1nnGAttVH5sjzqzHJBIirBSyhbK9t2BdaU6kHaTGdgtM/edit?usp=sharing).


## Types of Contribution

There are multiple ways to contribute to `netrd` (borrowed below list from [here](https://github.com/uzhdag/pathpy/blob/master/CONTRIBUTING.rst)):

#### Report Bugs

To report a bug in the package, open an issue at https://github.com/netsiphd/netrd/issues.

Please include in your bug report:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

#### Fix Bugs

Look through the GitHub issues for bugs. Anything tagged with "bug" and "help
wanted" is open to whoever wants to implement it.

#### Implement Features or New Methods

Look through the GitHub issues for features. Anything tagged with "enhancement"
and "help wanted" is open to whomever wants to implement it. If you know of a 
method that is implemented in another programming language, feel free to 
translate it into python here. If you don't want to translate it yourself, feel 
free to add an issue at https://github.com/netsiphd/netrd/issues. If you have 
read through this document and still have questions, also open an issue. When 
in doubt, open an issue.

#### Improve Documentation

Documentation is just as important as the code it documents. Please feel
free to submit PRs that are focused on fixing, improving, correcting, or
refactoring documentation. Documentation lives [here](https://netrd.readthedocs.io/en/latest/).

#### Submit Feedback

The best way to send feedback is to open an issue at https://github.com/netsiphd/netrd/issues.

If you are proposing to implement a distance metric, reconstruction algorithm, or dynamical process, 
see more details below. 

If you are proposing a feature not directly related to implementing a new method:

* Explain in detail why the feature is desirable and how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that your contributions
  are welcome!

##### A Brief Note On Licensing
Often, python code for an algorithm of interest already exists. In the interest of avoiding repeated reinvention of the wheel, we welcome code from other sources being integrated into `netrd`. If you are doing this, we ask that you be explicit and transparent about where the code came from and which license it is released under. The safest thing to do is copy the license from the original code into the header documentation of your file. For reference, this software is [licensed under MIT](https://github.com/tlarock/netrd/blob/master/LICENSE).

## Setup
Before starting your contribution, you need to complete the following instructions once.
The goal of this process is to fork, download and install the latest version of `netrd`.

1. Log in to GitHub.

2. Fork this repository by pressing 'Fork' at the top right of this
   page. This will lead you to 'github.com/<your_account>/netrd'. We refer
   to this as your personal fork (or just 'your fork'), as opposed to this repository
   (github.com/netsiphd/netrd), which we refer to as the 'upstream repository'.

3. Clone your fork to your machine by opening a console and doing

   ```
   git clone https://github.com/<your_account>/netrd.git
   ```

   Make sure to clone your fork, not the upstream repo. This will create a
   directory called 'netrd/'. Navigate to it and execute

   ```
   git remote add upstream https://github.com/netsiphd/netrd.git
   ```

   In this way, your machine will know of both your fork (which git calls
   `origin`) and the upstream repository (`upstream`).

4. During development, you will probably want to play around with your
   code. For this, you need to install the `netrd` package and have it
   reflect your changes as you go along. For this, open the console and
   navigate to the `netrd/` directory, and execute

	```
	pip install -e .
	```

	From now on, you can open a Jupyter notebook, ipython console, or your
    favorite IDE from anywhere in your computer and type `import netrd`.


These steps need to be taken only once. Now anything you do in the `netrd/`
directory in your machine can be `push`ed into your fork. Once it is in
your fork you can then request one of the organizers to `pull` from your
fork into the upstream repository (by submitting a 'pull request'). More on this later!


## Before you start coding

Once you have completed the above steps, you are ready to choose an algorithm to implement and begin coding.

1. Choose which algorithm you are interested in working on.

2. Open an issue at https://github.com/netsiphd/netrd/issues by clicking the "New Issue" button. 

	* Title the issue "Implement XYZ method", where XYZ method is a shorthand name for the distance, reconstruction or dynamics method you plan to implement.
	* Leave a comment that includes a brief motivation for why you want to see this method in `netrd`, as well as any key citations.
	* If such an issue already exists for the method you are going to write, it is not necessary to open another. However, it is a good idea to leave a comment letting others know you are going to work on it.

2. In your machine, create the file where your algorithm is going to
   live. If you chose a distance algorithm, copy
   an existing distance, such as `/netrd/distance/nbd.py`, into 
   `netrd/distance/<algorithm_name>.py`. Similarly, if you chose a reconstruction
   algorithm, copy an existing reconstruction method into 
   `netrd/reconstruction/<algorithm_name>.py`. Please keep in mind that
   <algorithm_name> will be used inside the code, so try to choose
   something that looks "pythonic". In particular, <algorithm_name> cannot
   include spaces, should not include upper case letters, and should use underscores 
   rather than hyphens.

3. Open the newly created file and edit as follows. At the very top you
   will find a string describing the algorithm. Edit this to describe the algorithm you
   are about to code, and preferably include a citation and link to any relevant papers. 
   Also add your name and email address (optional). Do not delete the line 
   `from .base import BaseDistance` or `from .base import BaseReconstructor`. 
   In the next line, change the class name to something appropriate. Guidelines here
   are to use `CamelCaseLikeThis` and not `snake_case_like_this`. (These
   are python guidelines, not ours!)

2. If you are implementing a distance method, you need to edit
   `netrd/distance/__init__.py`. Open it and add the following line:

	```
	from .<your_file_name> import <YourAlgorithmName>
	```

	* Note: there is one dot (.) before <your_file_name>. This is important!
	* Note: this line must go BEFORE the line with `__all__ = []`.

   If you are implementing a reconstruction method, you need to edit
    `netrd/reconstruction/__init__.py` instead, with the same line.

   This line tells the `netrd` package where to find your code.

3. In order for your contribution to have automated documentation
   generated, you need to edit the file
   `netrd/doc/source/distance.rst`. (Or
   `netrd/doc/source/reconstruction.rst`). Under the 'Available distances'
   heading you will find a list of currently implemented methods inside a
   `.. autosummary` command. Add the bane of your algorithm to this list,
   in alphabetical order.

4. Go back to editing <your_file_name>. After the line that starts with
   `class`, there is a function called `dist` for distances
   and `fit` for reconstructors. This is where the magic happens! There are
   some comments inside of these functions to guide the development, so
   please make sure to read them! You can also use the already existing code
   to get a feel for how you might design your implementation. Feel free to add 
   or remove anything and everything you feel you need. For example, if you 
   need auxiliary functions, you can add those as standalone functions. 
   Please try to base your style on already existing code as much as possible. 
   However, if you really need to do something differently, go ahead and we will 
   discuss how to make it fit with the rest of the package. In particular, 
   the `dist` or `fit` functions _must_ have those names and _must_ receive 
   those parameters (`dist` receives two graphs, while `fit` receives one time series). 
   You can add more parameters if your method needs them, but only _after_ the
   ones already in place. Every time you add a new parameter, it _must_ 
   have a default value.

5. If you need other auxiliary files, we got you covered! Place your
   Jupyter notebooks in the `netrd/notebooks` folder, and any data files
   you may need in the `netrd/data` folder. If you are willing to write a
   short documentation file (this may be plain text, a notebook, latex,
   etc) place that inside `netrd/docs`.
	NOTE: If you want to contribute these files to the git repository,
	please be mindful of the size and format of the files (especially
	when considering uploading data).


## After you finish coding

1. This project enforces a consistent coding style through the use of the
   [Black](https://black.readthedocs.io/en/stable/) autoformatter. Before
   committing your changes, please run Black to make sure your code passes our
   automated tests.

   ```
   # if you don't have black installed on your system, just pip install:
   # pip install black

   black --skip-string-normalization netrd/
   black --skip-string-normalization tests/
   ```

   If you do not do this, our tests will likely fail and your code will not be
   merged without further changes.

2. We also want to make sure there are no unused imports; for that we'll use
   [flake8](https://flake8.pycqa.org/en/latest/). This is a general-purpose
   linter that we're using to catch unused imports and wildcard imports. You
   can probably check these manually, but our test suite will run this command:

   ```
   # if you don't have flake8 installed on your system, just pip install:
   # pip install flake8

   flake8 --select F401,F403 netrd/
   flake8 --select F401,F403 tests/
   ```

   And if there are any import warnings, it cause the tests to fail and will
   have to be fixed before merging.

3. After updating your local code, the first thing to do is tell git which files
   you have been working on. (This is called staging.) If you worked on a
   distance algorithm, for example, do

   ```
   git add netrd/distance/<your_file> netrd/distance/__init__.py
   ```

4. Next tell git to commit (or save) your changes:

	```
	git commit -m 'Write a commit message here. This will be public and
	should be descriptive of the work you have done. Please be as explicit
	as possible, but at least make sure to include the name of the method
	you implemented. For example, the commit message may be: add
	implementation of SomeMethod, based on SomeAuthor and/or SomeCode.'
	```

5. Now you have to tell git to do two things. First, `pull` the latest changes from
   the upstream repository (in case someone made changes while you were coding), 
   then `push` your changes and the updated code from your machine to your fork:

	```
	git pull upstream master
	git push origin master
	```

	NOTE: If you edited already existing files, the `pull` may result in
	conflicts that must be merged. If you run in to trouble here, ask
	for help!

6. Finally, you need to tell this (the upstream) repository to include your
   contributions. For this, we use the GitHub web interface. At the top of
   this page, there is a 'New Pull Request' button. Click on it, and it
   will take you to a page titled 'Compare Changes'. Right below the title,
   click on the blue text that reads 'compare across forks'. This will show
   four buttons. Make sure that the first button reads 'base fork:
   netsiphd/netrd', the second button reads 'base: master', the third
   button reads 'head fork: <your_username>/netrd', and the fourth button
   reads 'compare: master'. (If everything has gone according to plan, the
   only button you should have to change is the third one - make sure you
   find your username, not someone else's.) After you find your username,
   GitHub will show a rundown of the differences that you are adding to the
   upstream repository, so you will be able to see what changes you are
   contributing. If everything looks correct, press 'Create Pull
   Request'.
	NOTE: Advanced git users may want to develop on branches other
	than master on their fork. That is totally fine, we won't know
	the difference in the end anyway.


That's it! After you've completed these steps, maintainers will be notified 
and will review your code and changes to make sure that everything is in place. 
Some automated tests will also run in the background to make sure that your 
code can be imported correctly and other sanity checks. Once that is all done, 
one of us will either accept your Pull Request, or leave a message requesting some
changes (you will receive an email either way).
---
title: 'netrd: A library for network reconstruction and graph distances'
tags:
  - Python
  - network science
  - network reconstruction
  - graph distance
  - network dynamics
authors:
  - name: Stefan McCabe
    orcid: 0000-0002-7180-145X
    affiliation: 1
  - name: Leo Torres
    orcid: 0000-0002-2675-2775
    affiliation: 1
  - name: Timothy LaRock
    orcid: 0000-0003-0801-3917
    affiliation: 1
  - name: Syed Arefinul Haque
    orcid: 0000-0002-8371-2366
    affiliation: 1
  - name: Chia-Hung Yang
    orcid: 0000-0002-4936-808X
    affiliation: 1
  - name: Harrison Hartle
    orcid: 0000-0002-0917-6112
    affiliation: 1
  - name: Brennan Klein
    orcid: 0000-0001-8326-5044
    affiliation: "1, 2"
affiliations:
 - name: Network Science Institute, Northeastern University, Boston, MA, USA
   index: 1
 - name: Laboratory for the Modeling of Biological and Socio-Technical Systems, Northeastern University, Boston, USA
   index: 2
date: 29 October 2020
bibliography: paper.bib
---

# Statement of need

Complex systems throughout nature and society are often best represented as *networks*. Over the last two decades, alongside the increased availability of large network datasets, we have witnessed the rapid rise of network science [@Amaral2004; @Vespignani2008; @Newman2010; @Barabasi2016]. This field is built around the idea that an increased understanding of the complex structural properties of a variety systems will allow us to better observe, predict, and even control the behavior of these systems.

However, for many systems, the data we have access to is not a direct description of the underlying network. More and more, we see the drive to study networks that have been inferred or reconstructed from non-network data---in particular, using *time series* data from the nodes in a system to infer likely connections between them [@Brugere2018; @Runge2018]. Selecting the most appropriate technique for this task is a challenging problem in network science. Different reconstruction techniques usually have different assumptions, and their performance varies from system to system in the real world. One way around this problem could be to use several different reconstruction techniques and compare the resulting networks. However, network comparison is also not an easy problem, as it is not obvious how best to quantify the differences between two networks, in part because of the diversity of tools for doing so.

The `netrd` Python package seeks to address these two parallel problems in network science.

# Summary

`netrd` offers, to our knowledge, the most extensive collection of both network reconstruction techniques and network comparison techniques (often referred to as *graph distances*) in a single library. Below, we expand on these two main functionalities of the `netrd` package.

The first core use of `netrd` is to reconstruct networks from time series data. Given time series data, $TS$, of the behavior of $N$ nodes / components / sensors of a system over the course of $L$ timesteps, and given the assumption that the behavior of every node, $v_i$, may have been influenced by the past behavior of other nodes, $v_j$, there are dozens of techniques that can be used to infer which connections, $e_{ij}$, are likely to exist between the nodes. That is, we can use one of many *network reconstruction* techniques to create a network representation, $G_r$, that attempts to best capture the relationships between the time series of every node in $TS$. `netrd` lets users perform this network reconstruction task using 17 different techniques. This means that up to 17 different networks can formed created from a single time series dataset. For example, in \autoref{fig:ground} we show the outputs of 15 different reconstruction techniques applied to time series data generated from an example network [@Sugihara2012; @Mishchenko2011; @Hoang2019; @Sheikhattar2018; @Friedman2008; @Edelman2005; @Zeng2013; @Donges2009; @Barucca2014; @Ledoit2003; @Stetter2012; @Peixoto2019].

Practitioners often apply these network reconstruction algorithms to real time series data. For example, in neuroscience, researchers often try to reconstruct functional networks from time series readouts of neural activity [@Mishchenko2011]. In economics, researchers can infer networks of influence between companies based on time series of changes in companies' stock prices [@Squartini2018]. At the same time, it is often quite helpful having the freedom to *simulate* arbitrary time series dynamics on randomly generated networks. This provides a controlled setting to assess the performance of network reconstruction algorithms. For this reason, the `netrd` package also includes a number of different techniques for simulating dynamics on networks.

The second core use of `netrd` addresses a common goal when studying networks: describing and quantifying the differences between two networks. This is a challenging problem, as there are countless axes upon which two networks can differ; as such, a number of *graph distance* measures have emerged over the years attempting to address this problem. As is the case for many hard problems in network science, it can be difficult to know which (of many) measures are suited for a given setting. In `netrd`, we consolidate over 20 different graph distance measures into a single package [@Jaccard1901; @Hamming1950; @Jurman2015; @Golub2013; @Donnat2018; @Carpi2011; @Bagrow2019; @DeDomenico2016; @Chen2018; @Hammond2013; @Monnig2018; @Tsitsulin2018; @Jurman2011; @Ipsen2002; @Torres2019; @Mellor2019; @Schieber2017; @Koutra2016; @Berlingerio2012]. \autoref{fig:dists} shows an example of just how different these measures can be when comparing two networks, $G_1$ and $G_2$. This submodule in `netrd` has already been used in recent work with a novel characterization of the graph distance literature [@Hartle2020].

This package builds on commonly used Python packages (e.g. `networkx` [@SciPyProceedings11], `numpy` [@Harris2020], `scipy` [@SciPy2020]) and is already a widely used resource for network scientists and other multidisciplinary researchers. With ongoing open-source development, we see this as a tool that will continue to be used by all sorts of researchers to come.

# Related software packages

In the network reconstruction literature, there are often software repositories that detail a single technique or a few related ones. For example Lizier (2014) implemented a Java package (portable to Python, octave, R, Julia, Clojure, MATLAB) that uses information-theoretic approaches for inferring network structure from time-series data [@Lizier2014]; Runge et al. (2019) created a Python package that combines linear or nonlinear conditional independence tests with a causal discovery algorithm to reconstruct causal networks from large-scale time series datasets [@Runge2019]. These are two examples of powerful and widely used packages though neither includes as wide-ranging techniques as `netrd` (nor were they explicitly designed to). In the graph distance literature, the same trend is broadly true: many one-off software repositories exist for specific measures. However, there are some packages that do include multiple graph distances; for example, Wills (2017) created a `NetComp` package that includes several variants of a few distance measures included here [@Wills2017].


# Figures
![**Example of the network reconstruction pipeline.** (Top row) A sample network, its adjacency matrix, and an example time series, $TS$, of node-level activity simulated on the network. (Bottom rows) The outputs of 15 different network reconstruction algorithms, each using $TS$ to create a new adjacency matrix that captures key structural properties of the original network.\label{fig:ground}](allRecons_withGroundtruth_SherringtonKirkpatrick.pdf)

![**Example of the graph distance measures in `netrd`.** Here, we measure the graph distance between two networks using 20 different distance measures from `netrd`.\label{fig:dists}](netrd_distance_example.pdf)

# Acknowledgements

The authors thank Kathryn Coronges, Mark Giannini, and Alessandro Vespignani for contributing to the coordination of the 2019 Network Science Institute "Collabathon", where much of the development of this package began. The authors acknowledge the support of ten other contributors to this package: Guillaume St-Onge, Andrew Mellor, Charles Murphy, David Saffo, Carolina Mattsson, Ryan Gallagher, Matteo Chinazzi, Jessica Davis, Alexander J. Gates, and Anton Tsitulin. **Funding:** This research was supported by the Network Science Institute at Northeastern University.

# References