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
---
title: "Marvel Social Network Example"
output: github_document
---
### Overview
This notebook contains an example of how to use the R version of this package to analyze the [Marvel character social network](https://arxiv.org/abs/cond-mat/0202174), a bipartite network of ties between characters and comic books. An edge between a character and a comic book indicates that the character appears in that comic book.

### Installation
Let's first load the package. Uncomment the first line if you have not yet installed the birankr package.

```{r}
#install.packages("birankr")
library(birankr)
```

We will also need to download the data and format the column names. We will use the cleaned edge list provided at http://syntagmatic.github.io/exposedata/marvel/data/source.csv. 

```{r}
marvel_df <- 
  fread("http://syntagmatic.github.io/exposedata/marvel/data/source.csv")
names(marvel_df) <- c('character', 'comic_book')

```

### Browse data
To get a sense of the data, let's print the number of unique characters and comic books, and print the first few lines of the data.


```{r}
marvel_df[, lapply(.SD, function(x) length(unique(x)))]
head(marvel_df)
```

### Estimate Ranks
Now we can try the bipartite_rank algorithm with two different normalizers: HITS and CoHITS. Because the first column of this data contains the `character` column, the `character` column will be treated as the senders (or top nodes), and the `comic_book` column will be treated as the receivers (or bottom nodes). Since the algorithm defaults to returning only the rankings of senders, the syntax below will provide us with rank scores for the `character` column.

```{r}
HITS_ranks <- bipartite_rank(data = marvel_df, normalizer='HITS')
CoHITS_ranks <- bipartite_rank(data = marvel_df, normalizer='CoHITS')
```

Notice that the results are slightly different, with the HITS normalizer returning Captain America as the highest ranked/most central comic book character, and the CoHITS normalizer returning Spider-Man as the highest ranked comic book character.

```{r}
head(HITS_ranks[order(HITS_ranks$rank, decreasing = T), ])
head(CoHITS_ranks[order(CoHITS_ranks$rank, decreasing = T), ])
```% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/06_PageRank.R
\name{pagerank}
\alias{pagerank}
\title{Estimate PageRank}
\usage{
pagerank(
  data,
  is_bipartite = TRUE,
  project_mode = c("rows", "columns"),
  sender_name = NULL,
  receiver_name = NULL,
  weight_name = NULL,
  rm_weights = FALSE,
  duplicates = c("add", "remove"),
  return_data_frame = TRUE,
  alpha = 0.85,
  max_iter = 200,
  tol = 1e-04,
  verbose = FALSE
)
}
\arguments{
\item{data}{Data to use for estimating PageRank. Can contain unipartite or bipartite graph data, either formatted as an edge list (class data.frame, data.table, or tibble (tbl_df)) or as an adjacency matrix (class matrix or dgCMatrix).}

\item{is_bipartite}{Indicate whether input data is bipartite (rather than unipartite/one-mode). Defaults to TRUE.}

\item{project_mode}{Mode for which to return PageRank estimates. Parameter ignored if is_bipartite = FALSE. Defaults to "rows" (the first column of an edge list).}

\item{sender_name}{Name of sender column. Parameter ignored if data is an adjacency matrix. Defaults to first column of edge list.}

\item{receiver_name}{Name of sender column. Parameter ignored if data is an adjacency matrix. Defaults to the second column of edge list.}

\item{weight_name}{Name of edge weights. Parameter ignored if data is an adjacency matrix. Defaults to edge weights = 1.}

\item{rm_weights}{Removes edge weights from graph object before estimating PageRank. Defaults to FALSE.}

\item{duplicates}{How to treat duplicate edges if any in data. Parameter ignored if data is an adjacency matrix. If option "add" is selected, duplicated edges and corresponding edge weights are collapsed via addition. Otherwise, duplicated edges are removed and only the first instance of a duplicated edge is used. Defaults to "add".}

\item{return_data_frame}{Return results as a data frame with node names in the first column and ranks in the second column. If set to FALSE, the function just returns a named vector of ranks. Defaults to TRUE.}

\item{alpha}{Dampening factor. Defaults to 0.85.}

\item{max_iter}{Maximum number of iterations to run before model fails to converge. Defaults to 200.}

\item{tol}{Maximum tolerance of model convergence. Defaults to 1.0e-4.}

\item{verbose}{Show the progress of this function. Defaults to FALSE.}
}
\value{
A dataframe containing each node name and node rank. If return_data_frame changed to FALSE or input data is classed as an adjacency matrix, returns a vector of node ranks. Does not return node ranks for isolates.
}
\description{
Estimate PageRank (centrality scores) of nodes from an edge list or adjacency matrix. If data is a bipartite graph, estimates PageRank based on a one-mode projection of the input. If the data is an edge list, returns ranks ordered by the unique values in the supplied edge list (first by unique senders, then by unique receivers).
}
\details{
The default optional arguments are likely well-suited for most users. However, it is critical to change the is.bipartite function to FALSE when working with one mode data. In addition, when estimating PageRank in unipartite edge lists that contain nodes with outdegrees or indegrees equal to 0, it is recommended that users append self-ties to the edge list to ensure that the returned PageRank estimates are ordered intuitively.
}
\examples{
#Prepare one-mode data
    df_one_mode <- data.frame(
      sender = sample(x = 1:10000, size = 10000, replace = TRUE),
      receiver = sample(x = 1:10000, size = 10000, replace = TRUE)
    )

#Add self-loops for all nodes
    unique_ids <- unique(c(df_one_mode$sender, df_one_mode$receiver))
    df_one_mode <- rbind(df_one_mode, data.frame(sender = unique_ids,
    receiver = unique_ids))

#Estimate PageRank in one-mode data
    PageRank <- pagerank(data = df_one_mode, is_bipartite = FALSE)

#Estimate PageRank in two-mode data
    df_two_mode <- data.frame(
      patient_id = sample(x = 1:10000, size = 10000, replace = TRUE),
      provider_id = sample(x = 1:5000, size = 10000, replace = TRUE)
    )
    PageRank <- pagerank(data = df_two_mode)
}
\references{
Lawrence Page, Sergey Brin, Rajeev Motwani, and Terry Winograd. "The pagerank citation ranking: Bringing order to the web". Technical report, Stanford InfoLab, 1999
}
\keyword{Bipartite}
\keyword{PageRank}
\keyword{centrality}
\keyword{rank}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/10_BiRank.R
\name{br_birank}
\alias{br_birank}
\title{BiRanks}
\usage{
br_birank(
  data,
  sender_name = NULL,
  receiver_name = NULL,
  weight_name = NULL,
  rm_weights = FALSE,
  duplicates = c("add", "remove"),
  return_mode = c("rows", "columns", "both"),
  return_data_frame = TRUE,
  alpha = 0.85,
  beta = 0.85,
  max_iter = 200,
  tol = 1e-04,
  verbose = FALSE
)
}
\arguments{
\item{data}{Data to use for estimating BiRank. Must contain bipartite graph data, either formatted as an edge list (class data.frame, data.table, or tibble (tbl_df)) or as an adjacency matrix (class matrix or dgCMatrix).}

\item{sender_name}{Name of sender column. Parameter ignored if data is an adjacency matrix. Defaults to first column of edge list.}

\item{receiver_name}{Name of sender column. Parameter ignored if data is an adjacency matrix. Defaults to the second column of edge list.}

\item{weight_name}{Name of edge weights. Parameter ignored if data is an adjacency matrix. Defaults to edge weights = 1.}

\item{rm_weights}{Removes edge weights from graph object before estimating BiRank. Parameter ignored if data is an edge list. Defaults to FALSE.}

\item{duplicates}{How to treat duplicate edges if any in data. Parameter ignored if data is an adjacency matrix. If option "add" is selected, duplicated edges and corresponding edge weights are collapsed via addition. Otherwise, duplicated edges are removed and only the first instance of a duplicated edge is used. Defaults to "add".}

\item{return_mode}{Mode for which to return BiRank ranks. Defaults to "rows" (the first column of an edge list).}

\item{return_data_frame}{Return results as a data frame with node names in first column and ranks in the second column. If set to FALSE, the function just returns a named vector of ranks. Defaults to TRUE.}

\item{alpha}{Dampening factor for first mode of data. Defaults to 0.85.}

\item{beta}{Dampening factor for second mode of data. Defaults to 0.85.}

\item{max_iter}{Maximum number of iterations to run before model fails to converge. Defaults to 200.}

\item{tol}{Maximum tolerance of model convergence. Defaults to 1.0e-4.}

\item{verbose}{Show the progress of this function. Defaults to FALSE.}
}
\value{
A dataframe containing each node name and node rank. If return_data_frame changed to FALSE or input data is classed as an adjacency matrix, returns a vector of node ranks. Does not return node ranks for isolates.
}
\description{
Estimate BiRanks of nodes from an edge list or adjacency matrix. Returns a vector of ranks or (optionally) a list containing a vector for each mode.
}
\details{
If input data is an edge list, this function returns ranks ordered by the unique values in the supplied edge list. Data inputted as an edge list are always assumed to contain named vertex IDs rather than to reflect an index of vertex positions in a network matrix. Users who wish for their edge lists to reflect vertex indices are recommended to input their data as a matrix or as a sparse matrix. \cr \cr
Network isolates are assigned a value of \eqn{(1 - alpha) / (n\_columns)} or \eqn{(1 - beta) / (n\_rows)} depending on their mode in the network. These values will always be smaller than the minimum value assigned to non-isolated nodes in the given mode. However, estimates on network isolates are non-meaningful. Users are advised to treat isolates with caution. \cr \cr
Created by He et al. (2017) \doi{10.1109/TKDE.2016.2611584}, BiRank is a highly generalizable algorithm that was developed explicitly for use in bipartite graphs. In fact, He et al.'s implementation of BiRank forms the basis of this package's implementation of all other bipartite ranking algorithms. Like every other bipartite ranking algorithm, BiRank simultaneously estimates ranks across each mode of the input data. BiRank's implementation is also highly similar to BGRM in that it symmetrically normalizes the transition matrix. BiRank differs from BGRM only in that it normalizes the transition matrix by the square-root outdegree of the source node and the square-root indegree of the target node.
}
\examples{
#create edge list between patients and providers
    df <- data.table(
      patient_id = sample(x = 1:10000, size = 10000, replace = TRUE),
      provider_id = sample(x = 1:5000, size = 10000, replace = TRUE)
    )

#estimate BiRank ranks
    BiRank <- br_birank(data = df)
}
\references{
Xiangnan He, Ming Gao, Min-Yen Kan, and Dingxian Wang. "Birank: Towards ranking on bipartite graphs". \emph{IEEE Transactions on Knowledge and Data Engineering}, 29(1):57-71, 2016
}
\keyword{BiRank}
\keyword{Bipartite}
\keyword{centrality}
\keyword{rank}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/05_bipartite_rank.R
\name{bipartite_rank}
\alias{bipartite_rank}
\title{Bipartite Ranks}
\usage{
bipartite_rank(
  data,
  sender_name = NULL,
  receiver_name = NULL,
  weight_name = NULL,
  rm_weights = FALSE,
  duplicates = c("add", "remove"),
  normalizer = c("HITS", "CoHITS", "BGRM", "BiRank"),
  return_mode = c("rows", "columns", "both"),
  return_data_frame = TRUE,
  alpha = 0.85,
  beta = 0.85,
  max_iter = 200,
  tol = 1e-04,
  verbose = FALSE
)
}
\arguments{
\item{data}{Data to use for estimating rank. Must contain bipartite graph data, either formatted as an edge list (class data.frame, data.table, or tibble (tbl_df)) or as an adjacency matrix (class matrix or dgCMatrix).}

\item{sender_name}{Name of sender column. Parameter ignored if data is an adjacency matrix. Defaults to first column of edge list.}

\item{receiver_name}{Name of sender column. Parameter ignored if data is an adjacency matrix. Defaults to the second column of edge list.}

\item{weight_name}{Name of edge weights. Parameter ignored if data is an adjacency matrix. Defaults to edge weights = 1.}

\item{rm_weights}{Removes edge weights from graph object before estimating rank. Parameter ignored if data is an edge list. Defaults to FALSE.}

\item{duplicates}{How to treat duplicate edges if any in data. Parameter ignored if data is an adjacency matrix. If option "add" is selected, duplicated edges and corresponding edge weights are collapsed via addition. Otherwise, duplicated edges are removed and only the first instance of a duplicated edge is used. Defaults to "add".}

\item{normalizer}{Normalizer (algorithm) used for estimating node ranks (centrality scores). Options include HITS, CoHITS, BGRM, and BiRank. Defaults to HITS.}

\item{return_mode}{Mode for which to return ranks. Defaults to "rows" (the first column of an edge list).}

\item{return_data_frame}{Return results as a data frame with node names in the first column and ranks in the second column. If set to FALSE, the function just returns a named vector of ranks. Defaults to TRUE.}

\item{alpha}{Dampening factor for first mode of data. Defaults to 0.85.}

\item{beta}{Dampening factor for second mode of data. Defaults to 0.85.}

\item{max_iter}{Maximum number of iterations to run before model fails to converge. Defaults to 200.}

\item{tol}{Maximum tolerance of model convergence. Defaults to 1.0e-4.}

\item{verbose}{Show the progress of this function. Defaults to FALSE.}
}
\value{
A dataframe containing each node name and node rank. If return_data_frame changed to FALSE or input data is classed as an adjacency matrix, returns a vector of node ranks. Does not return node ranks for isolates.
}
\description{
Estimate bipartite ranks (centrality scores) of nodes from an edge list or adjacency matrix. Functions as a wrapper for estimating rank based on a number of normalizers (algorithms) including HITS, CoHITS, BGRM, and BiRank. Returns a vector of ranks or (optionally) a list containing a vector for each mode.
}
\details{
If input data is an edge list, this function returns ranks ordered by the unique values in the supplied edge list. Data inputted as an edge list are always assumed to contain named vertex IDs rather than to reflect an index of vertex positions in a network matrix. Users who wish for their edge lists to reflect vertex indices are recommended to input their data as a matrix or as a sparse matrix. \cr \cr
Network isolates are assigned a value of \eqn{(1 - alpha) / (n\_columns)} or \eqn{(1 - beta) / (n\_rows)} depending on their mode in the network. These values will always be smaller than the minimum value assigned to non-isolated nodes in the given mode. However, estimates on network isolates are non-meaningful. Users are advised to treat isolates with caution. \cr \cr
For information about the different normalizers available in this function, see the descriptions for the HITS, CoHITS, BGRM, and BiRank functions. However, below outlines the key differences between the normalizers, with \eqn{K_d} and \eqn{K_p} representing diagonal matrices with generalized degrees (sum of the edge weights) on the diagonal (e.g. \eqn{(K_d)_{ii} = \sum_j w_{ij}} and \eqn{(K_p)_{jj} = \sum_i w_{ij}}).
\tabular{lll}{
\strong{Transition matrix} \tab \strong{\eqn{S_p}} \tab \strong{\eqn{S_d}} \cr
--------------------- \tab --------------------- \tab --------------------- \cr
HITS \tab \eqn{W^T} \tab \eqn{W} \cr
Co-HITS \tab \eqn{W^T K_d^{-1}} \tab \eqn{W K_p^{-1}} \cr
BGRM \tab \eqn{K_p^{-1} W^T K_d^{-1}} \tab \eqn{K_d^{-1} W K_p^{-1}} \cr
BiRank \tab \eqn{K_p^{-1/2} W^T K_d^{-1/2}} \tab \eqn{K_d^{-1/2} W K_p^{-1/2}}
}
}
\examples{
#create edge list between patients and providers
    df <- data.table(
      patient_id = sample(x = 1:10000, size = 10000, replace = TRUE),
      provider_id = sample(x = 1:5000, size = 10000, replace = TRUE)
    )

#estimate CoHITS ranks
    CoHITS <- bipartite_rank(data = df, normalizer = "CoHITS")
}
\keyword{BGRM}
\keyword{BiRank}
\keyword{Bipartite}
\keyword{CoHITS}
\keyword{HITS}
\keyword{centrality}
\keyword{rank}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/08_CoHITS.R
\name{br_cohits}
\alias{br_cohits}
\title{CoHITS Ranks}
\usage{
br_cohits(
  data,
  sender_name = NULL,
  receiver_name = NULL,
  weight_name = NULL,
  rm_weights = FALSE,
  duplicates = c("add", "remove"),
  return_mode = c("rows", "columns", "both"),
  return_data_frame = TRUE,
  alpha = 0.85,
  beta = 0.85,
  max_iter = 200,
  tol = 1e-04,
  verbose = FALSE
)
}
\arguments{
\item{data}{Data to use for estimating CoHITS. Must contain bipartite graph data, either formatted as an edge list (class data.frame, data.table, or tibble (tbl_df)) or as an adjacency matrix (class matrix or dgCMatrix).}

\item{sender_name}{Name of sender column. Parameter ignored if data is an adjacency matrix. Defaults to first column of edge list.}

\item{receiver_name}{Name of sender column. Parameter ignored if data is an adjacency matrix. Defaults to the second column of edge list.}

\item{weight_name}{Name of edge weights. Parameter ignored if data is an adjacency matrix. Defaults to edge weights = 1.}

\item{rm_weights}{Removes edge weights from graph object before estimating CoHITS. Parameter ignored if data is an edge list. Defaults to FALSE.}

\item{duplicates}{How to treat duplicate edges if any in data. Parameter ignored if data is an adjacency matrix. If option "add" is selected, duplicated edges and corresponding edge weights are collapsed via addition. Otherwise, duplicated edges are removed and only the first instance of a duplicated edge is used. Defaults to "add".}

\item{return_mode}{Mode for which to return CoHITS ranks. Defaults to "rows" (the first column of an edge list).}

\item{return_data_frame}{Return results as a data frame with node names in the first column and ranks in the second column. If set to FALSE, the function just returns a named vector of ranks. Defaults to TRUE.}

\item{alpha}{Dampening factor for first mode of data. Defaults to 0.85.}

\item{beta}{Dampening factor for second mode of data. Defaults to 0.85.}

\item{max_iter}{Maximum number of iterations to run before model fails to converge. Defaults to 200.}

\item{tol}{Maximum tolerance of model convergence. Defaults to 1.0e-4.}

\item{verbose}{Show the progress of this function. Defaults to FALSE.}
}
\value{
A dataframe containing each node name and node rank. If return_data_frame changed to FALSE or input data is classed as an adjacency matrix, returns a vector of node ranks. Does not return node ranks for isolates.
}
\description{
Estimate CoHITS ranks of nodes from an edge list or adjacency matrix. Returns a vector of ranks or (optionally) a list containing a vector for each mode.
}
\details{
If input data is an edge list, this function returns ranks ordered by the unique values in the supplied edge list. Data inputted as an edge list are always assumed to contain named vertex IDs rather than to reflect an index of vertex positions in a network matrix. Users who wish for their edge lists to reflect vertex indices are recommended to input their data as a matrix or as a sparse matrix. \cr \cr
Network isolates are assigned a value of \eqn{(1 - alpha) / (n\_columns)} or \eqn{(1 - beta) / (n\_rows)} depending on their mode in the network. These values will always be smaller than the minimum value assigned to non-isolated nodes in the given mode. However, estimates on network isolates are non-meaningful. Users are advised to treat isolates with caution. \cr \cr
Created by Deng, Lyo, and Kind (2009) \doi{10.1145/1557019.1557051}, CoHITS was developed explicitly for use in bipartite graphs as a way to better-incorporate content information (the "Co" in CoHITS) in HITS ranks. Like HITS, CoHITS is based on a markov process for simultaneously estimating ranks across each mode of the input data. CoHITS primarily differs from HITS in that it normalizes the transition matrix by the out-degree of the source nodes, leading to an interpretation more similar to that of a random walk.
}
\examples{
#create edge list between patients and providers
    df <- data.table(
      patient_id = sample(x = 1:10000, size = 10000, replace = TRUE),
      provider_id = sample(x = 1:5000, size = 10000, replace = TRUE)
    )

#estimate CoHITS ranks
    CoHITS <- br_cohits(data = df)
}
\references{
Hongbo Deng, Michael R. Lyu, and Irwin King. "A generalized co-hits algorithm and its application to bipartite graphs". In \emph{Proceedings of the 15th ACM SIGKDD International Conference on Knowledge Discovery and Data Mining},  KDD '09,  pages 239-248,  New York,  NY, USA, 2009. ACM.
}
\keyword{Bipartite}
\keyword{CoHITS}
\keyword{centrality}
\keyword{rank}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/04_sparsematrix_from_edgelist.R
\name{sparsematrix_from_edgelist}
\alias{sparsematrix_from_edgelist}
\title{Convert edge list to sparse matrix}
\usage{
sparsematrix_from_edgelist(
  data,
  sender_name = NULL,
  receiver_name = NULL,
  weight_name = NULL,
  duplicates = c("add", "remove"),
  is_bipartite = T
)
}
\arguments{
\item{data}{Edge list to convert to sparse matrix. Must be in edge list format and of class data.frame, data.table, or tbl_df.}

\item{sender_name}{Name of sender column. Defaults to the first column of an edge list.}

\item{receiver_name}{Name of sender column. Defaults to the second column of an edge list.}

\item{weight_name}{Name of edge weights. Defaults to edge weight = 1.}

\item{duplicates}{How to treat duplicate edges from edge list. If option "add" is selected, duplicated edges and corresponding edge weights are collapsed via addition. Otherwise, duplicated edges or removed and only the first instance of a duplicated edge is used. Defaults to "add".}

\item{is_bipartite}{Indicate whether input data is bipartite (rather than unipartite/one-mode). Defaults to TRUE.}
}
\value{
A sparse matrix of class dgCMatrix.
}
\description{
Converts edge lists (class data.frame) to sparse matrices (class "dgCMatrix"). For unipartite edge lists that contain any nodes with outdegrees or indegrees equal to 0, it is recommended that users append self-ties to the edge list to ensure that the IDs of the rows and columns are ordered intuitively to the user.
}
\examples{
#make edge.list
   df <- data.frame(
     id1 = sample(x = 1:20, size = 100, replace = TRUE),
     id2 = sample(x = 1:10, size = 100, replace = TRUE),
     weight = sample(x = 1:10, size = 100, replace = TRUE)
   )
#convert to sparsematrix
   sparsematrix_from_edgelist(data = df)
}
\keyword{dgCMatrix}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/02_sparsematrix_rm_weights.R
\name{sparsematrix_rm_weights}
\alias{sparsematrix_rm_weights}
\title{Remove sparse matrix edge weights}
\usage{
sparsematrix_rm_weights(adj_mat)
}
\arguments{
\item{adj_mat}{Sparse matrix of class dgCMatrix}
}
\value{
A sparse matrix of class dgCMatrix.
}
\description{
Removes edge weights from sparse matrices.
}
\examples{
#make matrix
   my_matrix <- sparseMatrix(
       i = c(1, 1, 2, 3, 4, 4, 5, 6, 7, 7), 
       j = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), 
       x = c(1, 1, 3, 1, 2, 1, 1, 1, 2, 1)
   )
#remove weights
   sparsematrix_rm_weights(my_matrix)
}
\keyword{dgCMatrix}
\keyword{matrix}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/07_HITS.R
\name{br_hits}
\alias{br_hits}
\title{HITS Ranks}
\usage{
br_hits(
  data,
  sender_name = NULL,
  receiver_name = NULL,
  weight_name = NULL,
  rm_weights = FALSE,
  duplicates = c("add", "remove"),
  return_mode = c("rows", "columns", "both"),
  return_data_frame = TRUE,
  alpha = 0.85,
  beta = 0.85,
  max_iter = 200,
  tol = 1e-04,
  verbose = FALSE
)
}
\arguments{
\item{data}{Data to use for estimating HITS. Must contain bipartite graph data, either formatted as an edge list (class data.frame, data.table, or tibble (tbl_df)) or as an adjacency matrix (class matrix or dgCMatrix).}

\item{sender_name}{Name of sender column. Parameter ignored if data is an adjacency matrix. Defaults to first column of edge list.}

\item{receiver_name}{Name of sender column. Parameter ignored if data is an adjacency matrix. Defaults to the second column of edge list.}

\item{weight_name}{Name of edge weights. Parameter ignored if data is an adjacency matrix. Defaults to edge weights = 1.}

\item{rm_weights}{Removes edge weights from graph object before estimating HITS. Parameter ignored if data is an edge list. Defaults to FALSE.}

\item{duplicates}{How to treat duplicate edges if any in data. Parameter ignored if data is an adjacency matrix. If option "add" is selected, duplicated edges and corresponding edge weights are collapsed via addition. Otherwise, duplicated edges are removed and only the first instance of a duplicated edge is used. Defaults to "add".}

\item{return_mode}{Mode for which to return HITS ranks. Defaults to "rows" (the first column of an edge list).}

\item{return_data_frame}{Return results as a data frame with node names in the first column and ranks in the second column. If set to FALSE, the function just returns a named vector of ranks. Defaults to TRUE.}

\item{alpha}{Dampening factor for first mode of data. Defaults to 0.85.}

\item{beta}{Dampening factor for second mode of data. Defaults to 0.85.}

\item{max_iter}{Maximum number of iterations to run before model fails to converge. Defaults to 200.}

\item{tol}{Maximum tolerance of model convergence. Defaults to 1.0e-4.}

\item{verbose}{Show the progress of this function. Defaults to FALSE.}
}
\value{
A dataframe containing each node name and node rank. If return_data_frame changed to FALSE or input data is classed as an adjacency matrix, returns a vector of node ranks. Does not return node ranks for isolates.
}
\description{
Estimate HITS ranks of nodes from an edge list or adjacency matrix. Returns a vector of ranks or (optionally) a list containing a vector for each mode.
}
\details{
If input data is an edge list, this function returns ranks ordered by the unique values in the supplied edge list. Data inputted as an edge list are always assumed to contain named vertex IDs rather than to reflect an index of vertex positions in a network matrix. Users who wish for their edge lists to reflect vertex indices are recommended to input their data as a matrix or as a sparse matrix. \cr \cr
Network isolates are assigned a value of \eqn{(1 - alpha) / (n\_columns)} or \eqn{(1 - beta) / (n\_rows)} depending on their mode in the network. These values will always be smaller than the minimum value assigned to non-isolated nodes in the given mode. However, estimates on network isolates are non-meaningful. Users are advised to treat isolates with caution. \cr \cr
Although originally designed for estimating ranks in unipartite graphs, HITS (Hyperlink-Induced Topic Search) is also one of the earliest bipartite ranking algorithms. Created by Jon Kleinberg (2009) \doi{10.1145/324133.324140} as an alternative to PageRank, HITS takes better account of the topology of bipartite networks by iteratively ranking nodes according to their role as an "Authority" and as a "Hub". Nodes with authority have high indegree from high ranking hubs; high ranking hubs have high outdegree to nodes with high authority. This function provides a slightly expanded version of HITS that only interfaces with bipartite networks and that allows for weighted edges. In general, HITS ranks tend to be more sensitive to user query than PageRanks, but HITS is substantially less efficient in ranking large graphs. HITS is likely less preferable than the other bipartite ranking algorithms in most applications. There are a number of contexts where HITS performs poorly, such as in graphs with extreme outliers.
}
\examples{
#create edge list between patients and providers
    df <- data.table(
      patient_id = sample(x = 1:10000, size = 10000, replace = TRUE),
      provider_id = sample(x = 1:5000, size = 10000, replace = TRUE)
    )

#estimate HITS ranks
    HITS <- br_hits(data = df)
}
\references{
Jon M. Kleinberg. "Authoritative sources in a hyperlinked environment". \emph{J. ACM}, 46(5):604-632, September 1999.
}
\keyword{Bipartite}
\keyword{HITS}
\keyword{centrality}
\keyword{rank}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/01_project_to_one_mode.R
\name{project_to_one_mode}
\alias{project_to_one_mode}
\title{Create a one-mode projection of a two mode graph}
\usage{
project_to_one_mode(adj_mat, mode = c("rows", "columns"))
}
\arguments{
\item{adj_mat}{Sparse matrix of class dgCMatrix}

\item{mode}{Mode to return. Defaults to projecting by rows.}
}
\value{
A double or complex matrix, with appropriate dimnames taken from x and y.
}
\description{
Create a one-mode projection of a two mode graph. Converts a rectangular matrix to a square one by taking the cross product of the input matrix. The edge weights in the resulting matrix are equal to the number of transitive ties of each node in the input matrix.
}
\examples{
#make matrix
   my_matrix <- sparseMatrix(i = c(1, 1, 2, 3, 4, 4, 5, 6, 7, 7), 
       j = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), x = 1
   )
#project to one mode
   project_to_one_mode(adj_mat = my_matrix, mode = "rows")
}
\keyword{dgCMatrix}
\keyword{matrix}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/09_BGRM.R
\name{br_bgrm}
\alias{br_bgrm}
\title{BGRM Ranks}
\usage{
br_bgrm(
  data,
  sender_name = NULL,
  receiver_name = NULL,
  weight_name = NULL,
  rm_weights = FALSE,
  duplicates = c("add", "remove"),
  return_mode = c("rows", "columns", "both"),
  return_data_frame = TRUE,
  alpha = 0.85,
  beta = 0.85,
  max_iter = 200,
  tol = 1e-04,
  verbose = FALSE
)
}
\arguments{
\item{data}{Data to use for estimating BGRM. Must contain bipartite graph data, either formatted as an edge list (class data.frame, data.table, or tibble (tbl_df)) or as an adjacency matrix (class matrix or dgCMatrix).}

\item{sender_name}{Name of sender column. Parameter ignored if data is an adjacency matrix. Defaults to first column of edge list.}

\item{receiver_name}{Name of sender column. Parameter ignored if data is an adjacency matrix. Defaults to the second column of edge list.}

\item{weight_name}{Name of edge weights. Parameter ignored if data is an adjacency matrix. Defaults to edge weights = 1.}

\item{rm_weights}{Removes edge weights from graph object before estimating BGRM. Parameter ignored if data is an edge list. Defaults to FALSE.}

\item{duplicates}{How to treat duplicate edges if any in data. Parameter ignored if data is an adjacency matrix. If option "add" is selected, duplicated edges and corresponding edge weights are collapsed via addition. Otherwise, duplicated edges are removed and only the first instance of a duplicated edge is used. Defaults to "add".}

\item{return_mode}{Mode for which to return BGRM ranks. Defaults to "rows" (the first column of an edge list).}

\item{return_data_frame}{Return results as a data frame with node names in the first column and ranks in the second column. If set to FALSE, the function just returns a named vector of ranks. Defaults to TRUE.}

\item{alpha}{Dampening factor for first mode of data. Defaults to 0.85.}

\item{beta}{Dampening factor for second mode of data. Defaults to 0.85.}

\item{max_iter}{Maximum number of iterations to run before model fails to converge. Defaults to 200.}

\item{tol}{Maximum tolerance of model convergence. Defaults to 1.0e-4.}

\item{verbose}{Show the progress of this function. Defaults to FALSE.}
}
\value{
A dataframe containing each node name and node rank. If return_data_frame changed to FALSE or input data is classed as an adjacency matrix, returns a vector of node ranks. Does not return node ranks for isolates.
}
\description{
Estimate BGRM ranks of nodes from an edge list or adjacency matrix. Returns a vector of ranks or (optionally) a list containing a vector for each mode.
}
\details{
If input data is an edge list, this function returns ranks ordered by the unique values in the supplied edge list. Data inputted as an edge list are always assumed to contain named vertex IDs rather than to reflect an index of vertex positions in a network matrix. Users who wish for their edge lists to reflect vertex indices are recommended to input their data as a matrix or as a sparse matrix. \cr \cr
Network isolates are assigned a value of \eqn{(1 - alpha) / (n\_columns)} or \eqn{(1 - beta) / (n\_rows)} depending on their mode in the network. These values will always be smaller than the minimum value assigned to non-isolated nodes in the given mode. However, estimates on network isolates are non-meaningful. Users are advised to treat isolates with caution. \cr \cr
Created by Rui et. al (2007) \doi{10.1145/1291233.1291378}, BGRM (Bipartite Graph Reinforcement Model) was developed explicitly for use in bipartite graphs. Like every bipartite ranking algorithm in this package, BGRM simultaneously estimates ranks across each mode of the input data. BGRM primarily differs from CoHITS and HITS by symmetrically normalizing the transition matrix, both by the out-degree of the source node and the indegree of the target node.
}
\examples{
#create edge list between patients and providers
    df <- data.table(
      patient_id = sample(x = 1:10000, size = 10000, replace = TRUE),
      provider_id = sample(x = 1:5000, size = 10000, replace = TRUE)
    )

#estimate BGRM ranks
    BGRM <- br_bgrm(data = df)
}
\references{
Xiaoguang Rui, Mingjing Li, Zhiwei Li, Wei-Ying Ma, and Nenghai Yu. "Bipartite graph reinforcement model for web image annotation". In \emph{Proceedings of the 15th ACM International Conference on Multimedia}, MM '07, pages 585-594, New York, NY, USA, 2007. ACM.
}
\keyword{BGRM}
\keyword{Bipartite}
\keyword{centrality}
\keyword{rank}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/03_sparsematrix_from_matrix.R
\name{sparsematrix_from_matrix}
\alias{sparsematrix_from_matrix}
\title{Convert matrix to sparse matrix}
\usage{
sparsematrix_from_matrix(adj_mat)
}
\arguments{
\item{adj_mat}{Adjacency matrix.}
}
\value{
A sparse matrix of class dgCMatrix.
}
\description{
Converts adjacency matrices (class "matrix") to a sparse matrices (class "dgCMatrix").
}
\examples{
#make matrix
   my_matrix <- rep(0, 100)
   my_matrix[c(1, 11, 22, 33, 44, 54, 65, 76, 87, 97)] <- 1
   my_matrix <- matrix(data = my_matrix, nrow = 10, ncol = 10)
#convert to sparsematrix
   sparsematrix_from_matrix(adj_mat = my_matrix)
}
\keyword{dgCMatrix}
\keyword{matrix}
