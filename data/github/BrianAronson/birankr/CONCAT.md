---
title: 'BiRank: Fast and Flexible Ranking on Bipartite Networks with R and Python'
tags:
  - bipartite network
  - PageRank
  - ranking
  - centrality
  - R
  - Python
authors:
  - name: Kai-Cheng Yang
    affiliation: 1
  - name: Brian Aronson
    affiliation: 2
  - name: Yong-Yeol Ahn
    affiliation: 1
affiliations:
 - name: Luddy School of Informatics, Computing, and Engineering, Indiana University, Bloomington, IN
   index: 1
 - name: Department of Sociology, Indiana University, Bloomington, IN
   index: 2
date: 18 February 2020
bibliography: paper.bib
---

# Summary

Bipartite (two-mode) networks are ubiquitous.
Common examples include networks of collaboration between scientists and their shared papers, networks of affiliation between corporate directors and board members, networks of patients and their doctors, and networks of competition between companies and their shared consumers.
Bipartite networks are commonly reduced to unipartite networks for further analysis, such as calculating node centrality (e.g. PageRank, see Figure \ref{illustration}(c)).
However, one-mode projections often destroy important structural information [@lehmann2008biclique] and can lead to imprecise network measurements. Moreover, there are numerous ways to obtain unipartite networks from a bipartite network, each of which has different characteristics and idiosyncrasies [@bass2013using].

To overcome the issues of one-mode projection, we present BiRank, an R and Python package that performs PageRank on bipartite networks directly.
The BiRank package contains several ranking algorithms that generalize PageRank to bipartite networks by propagating the probability mass (or importance scores) across two sides of the networks repeatedly using the following equations:

$$ \mathbf{T} = \alpha S_T \mathbf{B} + (1-\alpha)\mathbf{T}^0 $$
$$ \mathbf{B} = \beta S_B \mathbf{T} + (1-\beta)\mathbf{B}^0 $$

until they converge (see Figure \ref{illustration}(a)), where $\mathbf{T},\mathbf{B}$ are the ranking values for the top and bottom nodes, elements in $\mathbf{T}^0$ and $\mathbf{B}^0$ are set to $1/|\mathbf{T}|$ and $1/|\mathbf{B}|$ by default, $\alpha$ and $\beta$ are damping factors and set to 0.85 by default, $S_T, S_B$ are the transition matrices.
Unlike the one-mode projected PageRank, BiRank algorithms generate the ranking values for nodes from both sides simultaneously and take account of the full network topology without any information loss.

![(a) BiRank algorithms perform the ranking process on the bipartite networks directly and generate the ranking values for the top and bottom nodes simultaneously. (b) A bipartite network with three top nodes and four bottom nodes. (c) After the one-mode projection, a unipartite network of the bottom nodes is generated. PageRank can be performed to generate the ranking values of the bottom nodes. \label{illustration}](illustration.png) 

Our package implements the most notable and straightforward operationalizations of biparitite PageRanks including HITS [@kleinberg1999authoritative; @liao2014network], CoHITS [@deng2009generalized], BGRM [@rui2007bipartite], and Birank [@he2016birank].
The algorithms mainly differ in the way they normalize node ranks in the iterations (see Table \ref{normalizers}).

: A summary of transition matrices used in different BiRank algorithms.
$K_T$ and $K_B$ are diagonal matrices with generalized degrees (sum of the edge weights) on the diagonal, i.e.  $(K_T)_{ii} = \sum_j w_{ij}$ and $(K_B)_{jj} = \sum_i w_{ij}$.
$w_{ij}$ is the element on row $i$ and column $j$ of the bipartite network adjacency matrix $W^{|T|\times |B|}$.  \label{normalizers}

+-----------------------+--------------------------------+---------------------------+
| **Transition matrix** | $S_B$                          | $S_T$                     |
+-----------------------+--------------------------------+---------------------------+
| HITS                  | $W^\top$                       | $W$                       |
+-----------------------+--------------------------------+---------------------------+
| Co-HITS               | $W^\top K_T^{-1}$              | $W K_B^{-1}$              |
+-----------------------+--------------------------------+---------------------------+
| BGRM                  | $K_B^{-1} W^\top K_T^{-1}$     | $K_T^{-1} W K_B^{-1}$     |
+-----------------------+--------------------------------+---------------------------+
| Birank                | $K_B^{-1/2} W^\top K_T^{-1/2}$ | $K_T^{-1/2} W K_B^{-1/2}$ |
+-----------------------+--------------------------------+---------------------------+

Our guiding philosophy is to make the package as flexible as possible, given the diverse array of problems and data formats that are used in network analysis, while achieving good performance.
We therefore provide a number of convenience options for incorporating edge weights into rank estimations, estimating ranks on different types of input (edge lists, dense matrices, and sparse matrices), multiple file formats (as vectors, lists, or data frames), and for estimating PageRank on the one-mode projection of a network.
Moreover, this implementation uses efficient data storage and algorithms to ensure good performance and scalability.
For example, regardless of the algorithm of choice, it takes less than 10 seconds and less than 1GB of RAM to estimate ranks on a bipartite network containing half a million top nodes, more than two million bottom nodes, and about three million edges on a machine with 16 AMD EPYC 7000 series 2.5 GHz processors.

As a demonstration, we apply HITS, CoHITS, and one-mode projected PageRank to the Marvel Universe collaboration network  [@alberich2002marvel].
The Marvel Universe collaboration network comprises a network of affiliation with ties between every Marvel comic book (n = 12,849) and every character (n = 6,444) who appeared in those books. To give a sense of this network's structure, Figure \ref{marvel_network} illustrates a small sociogram of characters within ten comic books of this dataset.

![Sociogram of character-book ties within 10 comic books of the Marvel Universe collaboration network. \label{marvel_network}](marvel_network.png)

Table \ref{marvel_rank} presents the five characters with the highest ranking values from each algorithm.
Results are similar, with Captain America and Iron Man occurring in all three ranking algorithms. 
However, discrepancies arise from differences in the underlying ranking algorithms and how they interact with the network's structure. 
PageRank on the one mode projection first converts comic-character ties to character-character ties. Without information about the structure of characters-comic ties, PageRank mainly prioritizes nodes with a large number of transitive ties in the original network. For example, Wolverine has a higher PageRank than the Thing but Wolverine appears in much fewer comic books than the Thing. Instead, Wolverine's high PageRank is a result of his co-presence in comic books with large numbers of other characters. In contrast, the Thing tends to repeatedly appear in central comic books with other central characters in the Marvel universe, hence the Thing has a high CoHITS rank but a lower PageRank than Wolverine.


: Top five characters in the Marvel Universe collaboration network ranked by HITS, CoHITS and PageRank with one-mode projection. \label{marvel_rank}

+---------+-----------------+-----------------+------------------------+
| **Rank**| **HITS**        | **CoHITS**      |**Projection+PageRank** |
+=========+=================+=================+========================+
| 1st     | Captain America | Spider-man      | Captain America        |
+---------+-----------------+-----------------+------------------------+
| 2nd     | Iron man        | Captain America | Spider-man             |
+---------+-----------------+-----------------+------------------------+
| 3rd     | Thing           | Iron man        | Iron man               |
+---------+-----------------+-----------------+------------------------+
| 4th     | Human torch     | Hulk            | Wolverine              |
+---------+-----------------+-----------------+------------------------+
| 5th     | Mr. fantastic   | Thing           | Thor                   |
+---------+-----------------+-----------------+------------------------+

Differences between how HITS and CoHITS estimate ranks on the Marvel Universe collaboration network are more complicated. CoHITS normalizes the transition matrix by the outdegree of the source nodes, and therefore places somewhat less value on connections from highly connected characters and from highly connected comic books than HITS. As a result, CoHITS tends to assign higher ranks to characters who are connected to a more diverse array of comic books than does HITS. This difference is best illustrated by the inclusion of Mr. Fantastic in HITS' top-ranked characters and the inclusion of Spider Man in CoHITS' top-ranked characters: Spider Man appears in nearly twice as many comic books as Mr. Fantastic and collaborates with a significantly wider cast of characters than Mr. Fantastic; however, Mr. Fantastic tends to appear in highly central comic books with large character casts. It is open to interpretation as to which measure of centrality is better, but in many applications, we tend to prefer CoHITS over HITS as CoHITS ranks are less influenced by the presence of outliers with extreme degrees [@aronson_yang_odabas_ahn_perry_2020].

It is also worth mentioning that assigning different edge weights to the network can significantly affect ranking results.
Our package offers flexibility by allowing different combinations of algorithms and edge weights.
We leave the choice to the users' discretion.

Despite the ubiquity of bitpartite networks, bipartite PageRank algorithms are missing from the popular network packages, and our package serves to close this gap. Our target audience includes researchers and data scientists who deal with bipartite networks. To improve the accessibility, both R (birankr) and Python (birankpy) versions of the package are available. The documentation of BiRank consists of manual pages for its method functions, example usages, and unit tests.

# Acknowledgement

The authors acknowledge support from National Institute on Drug Abuse (grant R01 DA039928).

# References
# BiRank R and Python package

JOSS paper: [![DOI](https://joss.theoj.org/papers/10.21105/joss.02315/status.svg)](https://doi.org/10.21105/joss.02315)

Python package:
[![PyPI version](https://badge.fury.io/py/birankpy.svg)](https://badge.fury.io/py/birankpy)
[![Downloads](https://pepy.tech/badge/birankpy)](https://pepy.tech/project/birankpy)

R package:
[![Travis build status](https://travis-ci.org/BrianAronson/birankr.svg?branch=master)](https://travis-ci.org/BrianAronson/birankr)
[![R Downloads](https://cranlogs.r-pkg.org/badges/grand-total/birankr)](https://cranlogs.r-pkg.org/badges/grand-total/birankr)

Bipartite (two-mode) networks are ubiquitous.
When calculating node centrality measures in bipartite networks, a common approach is to apply PageRank on the one-mode projection of the network.
However, the projection can cause information loss and distort the network topology.
For better node ranking on bipartite networks, it is preferable to use a ranking algorithm that fully accounts for the topology of both modes of the network.

We present the BiRank package, which implements bipartite ranking algorithms HITS, CoHITS, BGRM, and BiRank.
BiRank provides convenience options for incorporating node-level weights into rank estimations, allowing maximum flexibility for different purpose.
It can efficiently handle networks with millions of nodes on a single midrange server.
Both R and Python versions are available.

## R version: birankr

### Overview

[CRAN package](https://cran.r-project.org/package=birankr) with highly efficient functions for estimating various rank (centrality) measures of nodes in bipartite graphs (two-mode networks) including HITS, CoHITS, BGRM, and BiRank. Also provides easy-to-use tools for incorporating or removing edge-weights during rank estimation, projecting two-mode graphs to one-mode, efficiently estimating PageRank in one-mode graphs, and for converting edgelists and matrices to sparseMatrix format. Best of all, the package's rank estimators can work directly with common formats of network data including edgelists (class `data.frame`, `data.table`, or `tbl_df`) and adjacency matrices (class `matrix` or `dgCMatrix`).

### Installation

This package can be directly installed via CRAN with `install.packages("birankr")`. Alternatively, newest versions of this package can be installed with `devtools::install_github("BrianAronson/birankr")`

### Example

Let's pretend we have a dataset (`df`) containing patient-provider ties (`patient_id` and `provider_id`) among providers that have ever prescribed an opioid:

```r
df <- data.frame(
    patient_id = sample(x = 1:10000, size = 10000, replace = T),
    provider_id = sample(x = 1:5000, size = 10000, replace = T)
)
```

We are interested in identifying patients who are likely doctor shopping. We assume that a highly central patient in the patient-doctor network is likely to be a person who is deliberately identifying more "generous" opioid prescribers. We therefore estimate a patients' rank in this network with the CoHITS algorithm:

```r
df.rank <- br_cohits(data = df)
```

Note that rank estimates are scaled according to the size of the network, with more nodes tending to result in smaller ranks. Due to this, it is often advisable to rescale rank estimates more interpretable numbers. For example, we could rescale such that the mean rank = 1 with the following data.table syntax:

```r
df.rank <- data.table(df.rank)
df.rank[, rank := rank/mean(rank)]
```

Finally, we decide to identify the IDs and ranks of the highest ranking patients in `df`:

```r
head(df.rank[order(rank, decreasing = T), ], 10)
```

For a more detailed example, check out [examples/Marvel_social_network.md](https://github.com/BrianAronson/birankr/blob/master/examples/Marvel_social_network.md), where we use the ranking algorithm to analyze the Marvel comic book social network.

### Function overview

Below is a brief outline of each function in this package:

- **bipartite_rank**
  - Estimates any type of bipartite rank.
- **br_bgrm**
  - Estimates ranks with BGRM algorithm
- **br_birank**
  - Estimates ranks with BiRank algorithm
- **br_cohits**
  - Estimates ranks with CoHITS algorithm
- **br_hits**
  - Estimates ranks with HITS algorithm
- **pagerank**
  - Estimates ranks with PageRank algorithm
- **project_to_one_mode**
  - Creates a one mode projection of a sparse matrix
- **sparsematrix_from_edgelist**
  - Creates a sparsematrix from an edgelist
- **sparsematrix_from_matrix**
  - Creates a sparsematrix from a matrix
- **sparsematrix_rm_weights**
  - Removes edge weights from a sparsematrix

### Documentation

Full documentation of `birankr` can be found in [birankr.pdf](https://github.com/BrianAronson/birankr/blob/master/birankr.pdf).

### Tests

To run the unit tests, install the birankr and devtools packages and run:

```
devtools::test("birankr")
```

## Python version: birankpy

### History

- Nov.10, 2021 (v1.0.1): drop support for python3.5; add support for python3.9

### Overview

`birankpy` provides functions for estimating various rank measures of nodes in bipartite networks including HITS, CoHITS, BGRM, and BiRank.
It can also project two-mode networks to one-mode, and estimate PageRank on it.
`birankpy` allows user-defined edge weights.
Implemented with sparse matrix, it's highly efficient.

### Dependencies

- `networkx`
- `pandas`
- `numpy`
- `scipy`

### Installation

Install with `pip`:

```bash
pip install birankpy
```

### Example

Let's pretend we have an edge list `edgelist_df` containing ties between top nodes and bottom nodes:

| top_node | bottom_node |
| -------- | ----------- |
| 1        | a           |
| 1        | b           |
| 2        | a           |
| ...      | ..          |
| 123      | z           |

To performing BiRank on this bipartite network, just:

```python
bn = birankpy.BipartiteNetwork()

bn.set_edgelist(edgelist_df,  top_col='top_node', bottom_col='bottom_node')

top_birank_df, bottom_birank_df = bn.generate_birank()
```

For a more detailed example, check out [examples/Marvel_social_network.ipynb](https://github.com/BrianAronson/birankr/blob/master/examples/Marvel_social_network.ipynb), where we use the ranking algorithm to analyze the Marvel comic book social network.

### Documentation

See documentation for `birankpy` at [birankpy doc](https://github.com/BrianAronson/birankr/blob/master/birankpy/README.md).

### Tests

To run the unit tests, first go to the `tests` directory and then run:

```bash
python test_birankpy.py
```

# Community Guidelines

## How to Contribute

In general, you can contribute to this project by creating [issues](https://github.com/BrianAronson/birankr/issues).
You are also welcome to contribute to the source code directly by forking the project, modifying the code, and creating [pull requests](https://github.com/BrianAronson/birankr/pulls).
If you are not familiar with pull requests, check out [this post](https://guides.github.com/activities/forking/).
Please use clear and organized descriptions when creating issues and pull requests.

## Bug Report and Support Request

You can use [issues](https://github.com/BrianAronson/birankr/issues) to report bugs and seek support.
Before creating any new issues, please check for similar ones in the issue list first.
This resubmission addresses the following concerns:

>    Found the following (possibly) invalid URLs:
>      URL: https://cran.r-project.org/web/packages/birankr/index.html
>        From: README.md
>        Status: 200
>        Message: OK
>        CRAN URL not in canonical form
>      The canonical URL of the CRAN page for a package is
>        https://CRAN.R-project.org/package=pkgname
> 
>    Found the following (possibly) invalid file URIs:
>      URI: examples/Marvel_social_network.md
>        From: README.md
>      URI: birankr.pdf
>        From: README.md
>      URI: examples/Marvel_social_network.ipynb
>        From: README.md
>      URI: birankpy/README.md
>        From: README.md

Done. I have changed links in the readme to use absolute path names and I provided the canonical URL for the CRAN page. 


## Test environments
* local Windows 10 install, R 3.4.1
* x86_64-redhat-linux-gnu, R 3.5.2

## R CMD check results
There were no ERRORs, WARNINGs, or NOTEs. 

## Downstream dependencies
There are currently no downstream dependencies for this package# BiRankpy

Bipartite (two-mode) networks are ubiquitous.
When calculating node centrality measures in bipartite networks, a common approach is to apply PageRank on the one-mode projection of the network.
However, the projection can cause information loss and distort the network topology.
For better node ranking on bipartite networks, it is preferable to use a ranking algorithm that fully accounts for the topology of both modes of the network.

We present the BiRank package, which implements bipartite ranking algorithms HITS, CoHITS, BGRM, and Birank.
BiRank provides convenience options for incorporating node-level weights into rank estimations, allowing maximum flexibility for different purpose.
It can efficiently handle networks with millions of nodes on a single midrange server.
Both R and Python versions.


# Overview

`birankpy` provides functions for estimating various rank measures of nodes in bipartite networks including HITS, CoHITS, BGRM, and Birank.
It can also project two-mode networks to one-mode, and estimat PageRank on it.
`birankpy` allows user-defined edge weights.
Implemented with sparse matrix, it's highly efficient.

### Example
Let's pretend we have a edge list `edgelist_df` containing ties between top nodes and bottom nodes:

top_node | bottom_node
------------ | -------------
1 | a
1 | b
2 | a


To performing BiRank on, just:

```python
bn = birankpy.BipartiteNetwork()

bn.set_edgelist(edgelist_df,  top_col='top_node', bottom_col='bottom_node')

top_birank_df, bottom_birank_df = bn.generate_birank()
```
﻿# birankr 1.0.1
- Fixed a bug preventing rank algorithms from working properly on dense matrices.
- Fixed a bug sometimes causing wrong ID names to be returned in dataframes
- Fixed a bug preventing PageRank from working on some types of sparse matrices.
- Fixed a bug preventing proper conversion of unipartite dense matrices to sparse matrices.
- Minor syntax changes within bipartite ranking functions to make them more error-proof.
- Added unit tests to prevent bugs from occur following future updates.# API documentation of birankpy

## PageRank

```python
pagerank(adj, d=0.85, max_iter=200, tol=0.0001, verbose=False)
```

Return the PageRank of the nodes in a graph using power iteration.
This function takes the sparse matrix as input directly, avoiding the overheads of converting the network to a networkx Graph object and back.

This function takes the sparse matrix as input directly, avoiding the overheads
of converting the network to a networkx Graph object and back.

Input:

parameter | type | note
----------|------|-----
adj | scipy.sparsematrix | Adjacency matrix of the graph
d   | float              | Dumping factor
max_iter | int | Maximum iteration times
tol | float | Error tolerance to check convergence
verbose | boolean | If print iteration information

Output:

type | note
----------|------
numpy.ndarray|The PageRank values

## BiRank
```python
birank(W, normalizer='HITS', alpha=0.85, beta=0.85, max_iter=200, tol=0.0001, verbose=False)
```

Calculate the PageRank of bipartite networks directly.
See paper https://ieeexplore.ieee.org/abstract/document/7572089/
for details.
Different normalizer yields very different results.
More studies are needed for deciding the right one.

Input:

parameter | type | note
----------|------|-----
W | scipy.sparsematrix | Adjacency matrix of the bipartite network D\*P
normalizer | string | Choose which normalizer to use, see the paper for details
alpha | float | Damping factors for the rows
beta | float | Damping factors for the columns
max_iter | int | Maximum iteration times
tol | float | Error tolerance to check convergence
verbose | boolean | If print iteration information

Output:

variable | type | note
----------|------|-----
d | numpy.ndarray | The BiRank for rows and columns
p | numpy.ndarray | The BiRank for rows and columns

## UnipartiteNetwork
```python
UnipartiteNetwork(self)
```

Class for handling unipartite networks using scipy's sparse matrix
Designed to for large networkx, but functionalities are limited

###  

```python
set_adj_matrix(self, id_df, W, index_col=None)
```

Set the adjacency matrix of the network.

Input:

parameter | type | note
----------|------|-----
id_df | pandas.DataFrame | The mapping between node and index
W | scipy.sparsematrix | The adjacency matrix of the network; the node order in id_df has to match W
index_col | string | column name of the index


```python
generate_pagerank(self, **kwargs)
```

Generate the PageRank values for the network using `pagerank()`.
The parameters are the same with `pagerank()`.


## BipartiteNetwork
```python
BipartiteNetwork(self)
```

Class for handling bipartite networks using scipy's sparse matrix
Design to for large networkx, but functionalities are limited

```python
load_edgelist(self, edgelist_path, top_col, bottom_col, weight_col='None', sep=',')
```

Method to load an edgelist.

Inputs:

parameter | type | note
----------|------|-----
edge_list_path | string | the path to the edgelist file
top_col | string | column of the top nodes
bottom_col | string | column of the bottom nodes
weight_col | string | column of the edge weights
sep | string | the seperators of the edgelist file

Suppose the bipartite network has D top nodes and P bottom nodes.
The edgelist file should have the format similar to the example:

top | bottom | weight
----|--------|-------
t1 | b1 | 1
t1 | b2 | 1
t2 | b1 | 2
...|...|...
tD | bP | 1

The edgelist file needs at least two columns for the top nodes and bottom nodes. An optional column can carry the edge weight.
You need to specify the columns in the method parameters.
The network is represented by a D*P dimensional matrix.

```python
set_edgelist(self, df, top_col, bottom_col, weight_col=None)
```

Method to set the edgelist.

Inputs:

parameter | type | note
----------|------|-----
df | pandas.DataFrame | the edgelist with at least two columns
top_col | string | column of the edgelist dataframe for top nodes
bottom_col | string | column of the edgelist dataframe for bottom nodes
weight_col | string | column of the edgelist dataframe for edge weights


The edgelist should be represented by a dataframe.
The dataframe needs at least two columns for the top nodes and bottom nodes. An optional column can carry the edge weight.
You need to specify the columns in the method parameters.

```python
unipartite_projection(self, on)
```
Project the bipartite network to one side of the nodes to generate a unipartite network.

Input:

parameter | type | note
----------|------|-----
on | string | Name of the column to project the network on

Output:

| type | note
------|-----
UnipartiteNetwork | The projected unipartite network

If projected on top nodes, the resulting adjacency matrix has dimension: D\*D.
If projected on bottom nodes, the resulting adjacency matrix has dimension: P\*P.

```python
generate_birank(self, **kwargs)
```

Output:

variable | type | note
---------|------|-----
top_df   | pandas.DataFrame | BiRank values for the top nodes
bottom_df   | pandas.DataFrame | BiRank values for bottom nodes


Generate the BiRank values for the top and bottom nodes simultaneously using `birank()`.
The parameters are the same with `birank()`.
This directory contains examples of using the Python and R packages to analyze the Marvel social network.

`Marvel_social_network.ipynb` is for the Python version; `Marvel_social_network.Rmd` and `Marvel_social_network.md` are for the R version.
Marvel Social Network Example
================

### Overview

This notebook contains an example of how to use the R version of this
package to analyze the [Marvel character social
network](https://arxiv.org/abs/cond-mat/0202174), a bipartite network of
ties between characters and comic books. An edge between a character and
a comic book indicates that the character appears in that comic book.

### Installation

Let’s first load the package. Uncomment the first line if you have not
yet installed the birankr package.

``` r
#install.packages("birankr")
library(birankr)
```

    ## Loading required package: Matrix

    ## Loading required package: data.table

We will also need to download the data and format the column names. We
will use the cleaned edge list provided at
<http://syntagmatic.github.io/exposedata/marvel/data/source.csv>.

``` r
marvel_df <- 
  fread("http://syntagmatic.github.io/exposedata/marvel/data/source.csv")
names(marvel_df) <- c('character', 'comic_book')
```

### Browse data

To get a sense of the data, let’s print the number of unique characters
and comic books, and print the first few lines of the data.

``` r
marvel_df[, lapply(.SD, function(x) length(unique(x)))]
```

    ##    character comic_book
    ## 1:      6444      12849

``` r
head(marvel_df)
```

    ##               character comic_book
    ## 1: KILLRAVEN/JONATHAN R     AA2 35
    ## 2:             M'SHULLA     AA2 35
    ## 3: 24-HOUR MAN/EMMANUEL     AA2 35
    ## 4:            OLD SKULL     AA2 35
    ## 5:               G'RATH     AA2 35
    ## 6: 3-D MAN/CHARLES CHAN   M/PRM 35

### Estimate Ranks

Now we can try the bipartite\_rank algorithm with two different
normalizers: HITS and CoHITS. Because the first column of this data
contains the `character` column, the `character` column will be treated
as the senders (or top nodes), and the `comic_book` column will be
treated as the receivers (or bottom nodes). Since the algorithm defaults
to returning only the rankings of senders, the syntax below will provide
us with rank scores for the `character` column.

``` r
HITS_ranks <- bipartite_rank(data = marvel_df, normalizer='HITS')
CoHITS_ranks <- bipartite_rank(data = marvel_df, normalizer='CoHITS')
```

Notice that the results are slightly different, with the HITS normalizer
returning Captain America as the highest ranked/most central comic book
character, and the CoHITS normalizer returning Spider-Man as the highest
ranked comic book character.

``` r
head(HITS_ranks[order(HITS_ranks$rank, decreasing = T), ])
```

    ##                character       rank
    ## 14       CAPTAIN AMERICA 0.02703070
    ## 63  IRON MAN/TONY STARK  0.01993319
    ## 34  THING/BENJAMIN J. GR 0.01990385
    ## 115 HUMAN TORCH/JOHNNY S 0.01924820
    ## 113 MR. FANTASTIC/REED R 0.01882251
    ## 26  INVISIBLE WOMAN/SUE  0.01772762

``` r
head(CoHITS_ranks[order(CoHITS_ranks$rank, decreasing = T), ])
```

    ##                character        rank
    ## 80  SPIDER-MAN/PETER PAR 0.014299452
    ## 14       CAPTAIN AMERICA 0.011231753
    ## 63  IRON MAN/TONY STARK  0.009827441
    ## 23  HULK/DR. ROBERT BRUC 0.007843040
    ## 34  THING/BENJAMIN J. GR 0.007838295
    ## 105 THOR/DR. DONALD BLAK 0.007142875
