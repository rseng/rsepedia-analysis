### 0.3.3
- ROV (Range of Value) Method implemented.

### 0.3.2
- SD method implemented for determining weights.

### 0.3.1
- fix Copeland.
- add Moora Ratio method with new tests.

### 0.3.0
- entropy() returns a result even though there are NaNs for some criterion.
- rwrapper.R added so the library is callable from within R
- copeland() method added for combining multiple ordering results.
  

### 0.2.9
- Default optimizer is now GLPK (Cbc removed)


### 0.2.8
- Direction of optimization added for nds()
- New tests added

### 0.2.7
- Bug in Moore fixed.
- Bug in Marcos fixed.


### 0.2.6
- Bug in Electre fixed.
- New tests added.
- Tests were divided into several files
- Dependencies upgraded


### 0.2.5
- On/Off switch for tests. 
- New tests.

### 0.2.4
- Base.show(io:IO, MCDMResult) implementations for pretty printing all of the results
- I() implemented, LinearAlgebra package removed from dependencies.
- mean(), geomean(), std(), and cor() are implemented and StatsBase & Statistics packages are removed
- using keyword replaced by import and only the needed functions are loaded at startup.
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03430/status.svg)](https://doi.org/10.21105/joss.03430)

# JMcDM
A package for Multiple-criteria decision-making techniques in Julia.

## The problem

Suppose a decision process has n alternatives and m criteria  which are either to be maximized or minimized. Each single criterion has a weight _0 ≤ wᵢ ≤ 1_ where sum of _wᵢ_ is 1. _fᵢ_ is either maximum or minimum. _gⱼ(.)_ is the evolution function and it is choosen as _gⱼ(x) = x_ in many methods. A multiple criteria decision problem can be represented using the decision table 

<img src="https://github.com/jbytecode/JMcDM/blob/gh-pages/images/generalformula.png" width = "50%"/>

<!--
   | **Criteria**  |   C_1    |   C_2    |  ...  |   C_m    |
   | :-----------: | :------: | :------: | :---: | :------: |
   |  **Weights**  |   w_1    |   w_2    |  ...  |   w_m    |
   | **Functions** |   f_1    |   f_2    |  ...  |   f_m    |
   |      A_1      | g_1(A_1) | g_2(A_1) |  ...  | g_m(S_A) |
   |      A_2      | g_1(A_2) | g_2(A_2) |  ...  | g_m(A_2) |
   |       ⋮       |    ⋮     |    ⋮     |  ...  |    ⋮     |
   |      A_n      | g_1(A_n) | g_2(A_n) |  ...  | g_m(A_n) |
-->

without loss of generality. When _A₁, A₂, ..., Aₙ_ are alternatives and _C₁, C₂, ..., Cₙ_ are different situations of a single criterion then the decision problem is said to be single criterion decision problem. If _Cⱼ_ are strategies of two game players then _gⱼ(Aᵢ)_ is the gain of the row player when she selects the strategy _i_ and the column player selects the strategy _Cⱼ_.


The package mainly focuses on solving these kinds of decision problems.

## For whom?

Multiple-criteria decision-making is an inter-discipline subject and there is a vast amount of research in the literature in this area. However, the existing software packages in this area are generally focused on a small subset of tools. JMcDM is a developer and researcher-friendly Julia package that combines the developed methods, utility functions for implementing new ones, and serves an environment for comparing results of multiple analyses.  

## Installation

Please type 

```julia
julia> ]
(@v1.5) pkg> add JMcDM
```

or

```julia
julia> using Pkg
julia> Pkg.add("JMcDM")
```

in Julia REPL. 


## Package Dependencies

Since the Julia package manager installs all of the dependencies automatically, a standard user doesn't need to
install them manually. The package dependencies are listed below:

- DataFrames
- GLPK
- JuMP


## Documentation

Please check out the reference manual [here](https://jbytecode.github.io/JMcDM/docs/build/).


## Implemented methods

### MCDM Tools

- TOPSIS (Technique for Order Preference by Similarity to Ideal Solutions)
- ELECTRE (Elimination and Choice Translating Reality)
- DEMATEL (The Decision Making Trial and Evaluation Laboratory)
- MOORA Reference (Multi-Objective Optimization By Ratio Analysis)
- MOORA Ratio
- VIKOR (VlseKriterijumska Optimizcija I Kaompromisno Resenje in Serbian)
- AHP (Analytic Hierarchy Process)
- DEA (Data Envelopment Analysis)
- GRA (Grey Relational Analysis)
- Non-dominated Sorting 
- SAW (Simple Additive Weighting) (aka WSM)
- ARAS (Additive Ratio Assessment)
- WPM (Weighted Product Model)
- WASPAS (Weighted Aggregated Sum Product ASsessment)
- EDAS (Evaluation based on Distance from Average Solution)
- MARCOS (Measurement Alternatives and Ranking according to COmpromise Solution)
- MABAC (Multi-Attributive Border Approximation area Comparison)
- MAIRCA (Multi Attributive Ideal-Real Comparative Analysis)
- COPRAS (COmplex PRoportional ASsessment)
- PROMETHEE (Preference Ranking Organization METHod for Enrichment of Evaluations)
- CoCoSo (Combined Compromise Solution)
- CRITIC (CRiteria Importance Through Intercriteria Correlation)
- Entropy
- CODAS (COmbinative Distance-based ASsessment)
- Copeland (For combining multiple ordering results)
- SD Method for determining weights of criteria
- ROV (Range of Value) Method

### SCDM Tools

- minimax
- maximin
- minimin
- maximax
- Savage
- Hurwicz
- MLE
- Laplace
- Expected Regret

### Game

- Game solver for zero sum games


## Unimplemented methods
- UTA
- MAUT
- STEM
- PAPRIKA
- ANP (Analytical Network Process)
- Goal Programming
- MACBETH
- COMET

- will be updated soon. 

## Example

```julia
julia> using JMcDM
julia> df = DataFrame(
:age        => [6.0, 4, 12],
:size       => [140.0, 90, 140],
:price      => [150000.0, 100000, 75000],
:distance   => [950.0, 1500, 550],
:population => [1500.0, 2000, 1100]);
```


```julia
julia> df
3×5 DataFrame
 Row │ age      size     price     distance  population 
     │ Float64  Float64  Float64   Float64   Float64    
─────┼──────────────────────────────────────────────────
   1 │     6.0    140.0  150000.0     950.0      1500.0
   2 │     4.0     90.0  100000.0    1500.0      2000.0
   3 │    12.0    140.0   75000.0     550.0      1100.0
```


```julia
julia> w  = [0.35, 0.15, 0.25, 0.20, 0.05];
julia> fns = makeminmax([minimum, maximum, minimum, minimum, maximum]);
julia> result = topsis(df, w, fns);
julia> result.scores
3-element Array{Float64,1}:
0.5854753145549456
0.6517997936899308
0.41850223305822903

julia> result.bestIndex
2
```

alternatively

```julia
julia> result = mcdm(df, w, fns, TopsisMethod())
```

or 

```julia
julia> setting = MCDMSetting(df, w, fns)
julia> result = topsis(setting)
```

or

```julia
julia> setting = MCDMSetting(df, w, fns)
julia> result = mcdm(setting, TopsisMethod())
```

### Jupyter Notebook

Here is a Jupyter Notebook for basic usage: 

https://github.com/jbytecode/JMcDM/blob/main/notebook/basic-usage.ipynb


## Community guidelines

### How to cite 

Please use the BibTeX entry:

```bibtex
@article{Satman2021,
  doi = {10.21105/joss.03430},
  url = {https://doi.org/10.21105/joss.03430},
  year = {2021},
  publisher = {The Open Journal},
  volume = {6},
  number = {65},
  pages = {3430},
  author = {Mehmet Hakan Satman and Bahadır Fatih Yıldırım and Ersagun Kuruca},
  title = {JMcDM: A Julia package for multiple-criteria decision-making tools},
  journal = {Journal of Open Source Software}
}
```

or citation string

Satman et al., (2021). JMcDM: A Julia package for multiple-criteria decision-making tools. Journal of Open Source Software, 6(65), 3430, https://doi.org/10.21105/joss.03430

to cite this software.

### Contribute to software
Do you want to contribute?

- Please create an issue first. In this issue, please specify the idea.
- Let the community discuss the new contribution in our Slack channel or the created issue.

If the community decision is yes, please

- Fork the repository
- Add the new code to this forked repository
- Make sure the tests are passed 
- Send a pull request with a good description of functionality.

### Where to start?
The TOPSIS method, defined in [topsis.jl](https://github.com/jbytecode/JMcDM/blob/main/src/topsis.jl), is a basis for many methods and it can be followed before implementing a new one. 

### The design pattern
- ```topsis()``` takes the decision matrix, weights, and vector of directions of optimization as arguments. This function is defined in ```topsis.jl```.

```julia
   function topsis(decisionMat::DataFrame, weights::Array{Float64,1}, fns::Array{Function,1})::TopsisResult
```

- ```topsis()``` method has a return type of ```TopsisResult```. This ```struct``` is defined in ```types.jl```

```julia
  struct TopsisResult <: MCDMResult
    decisionMatrix::DataFrame
    weights::Array{Float64,1}
    normalizedDecisionMatrix::DataFrame
    normalizedWeightedDecisionMatrix::DataFrame 
    bestIndex::Int64 
    scores::Array{Float64,1}
end
```

- Optionally, a ```show``` function can be derived for pretty-printing the result. These functions are defined in ```print.jl```

```julia
function Base.show(io::IO, result::TopsisResult)
    println(io, "Scores:")
    println(io, result.scores)
    println(io, "Best indices:")
    println(io, result.bestIndex)
end
```

Please read the issue [Welcome to newcomers!](https://github.com/jbytecode/JMcDM/issues/3) for other implementation details.

### Report Issues

If you find a bug or error, first report the problem in a new issue. If the problem is already addressed
in an existing issue please follow the existing one.

### Seek Support
Our Slack channel is [JMcDM Slack Channel](https://julialang.slack.com/archives/C01MJ0VF1U3). Please feel free to ask about any problem using our Slack channel or issues. [Julia Discourse](https://discourse.julialang.org/t/jmcdm-a-julia-package-for-multiple-criteria-decision-making-tools/54942) is the JMcDM entry in Julia Discourse site and any thoughts, problems, and issues can also be discussed there.


Welcome!
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
reported by contacting the project team at mhsatman@gmail.com. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see
https://www.contributor-covenant.org/faq
Do you want to contribute?

- Please create an Issue first. In this issue, please specify the idea.
- Let the community discuss the new contribution in our Slack channel or in the created issue.

If the community decision is yes, please

- Fork the repository
- Send a pull request.

Please read the issue [Welcome to newcomers!](https://github.com/jbytecode/JMcDM/issues/3) for implementation details.

Our Slack channel is [JMcDM Slack Channel](https://julialang.slack.com/archives/C01MJ0VF1U3).

Welcome!


---
title: 'JMcDM: A Julia package for multiple-criteria decision-making tools'
tags:
  - Julia
  - decision making
  - multiple criteria
  - outranking
authors:
  - name: Mehmet Hakan Satman
    orcid: 0000-0002-9402-1982
    affiliation: 1
  - name: Bahadır Fatih Yıldırım
    orcid: 0000-0002-0475-741X
    affiliation: 2
  - name: Ersagun Kuruca
    orcid: 0000-0002-2552-7701
    affiliation: 3
affiliations:
 - name: Department of Econometrics, Istanbul University, Istanbul, Turkey
   index: 1
 - name: Department of Transportation and Logistics, Istanbul University, Istanbul, Turkey
   index: 2
 - name: Independent researcher
   index: 3


date: 8 May 2021
bibliography: paper.bib
---

# Summary
```JMcDM``` is a ```Julia``` package that implements some leading multiple-criteria decision-making tools for both researchers and developers. By having a REPL tool, ```Julia``` is well suited for researchers to perform their analysis using different methods and comparing their results. ```JMcDM``` also provides the necessary infrastructure, utility functions, and a standardized API for implementing recently published methods.  The package brings MCDM (Multiple-Criteria Decision-Making) tools to a relatively new language such as ```Julia``` with its significant performance promises. Besides Julia being a new language, the methods developed in the package are designed to be familiar to users who previously used the ```R``` and ```Python``` languages. This paper presents the basics of the design, example usage, and code snippets.

# Introduction
The one-dimensional array $a$ is in ascending order if and only if $a_i \le a_{i+1}$ where $i = 1, 2, \dots, n-1$, and $n$ is the length of the array. In other terms, the process of ordering numbers requires the logical $\le$ operator to be perfectly defined. Since the operator $\le$ is not defined for any set of points in higher dimensions, $\mathbb{R}^p$ for $p \ge 2$, there is not a unique ordering of points. In the multi-dimensional case, the binary domination operator $\succ$ applied on points $a$ and $b$, $a \succ b$, is true if each item in $a$ is not worse than the corresponding item in $b$ and at least one item is better than the corresponding item in $b$ [@Deb_2002]. On the other hand, the more relaxed operator $\succeq$ returns true if each item in $a$ is as good as the corresponding item in $b$ [@greco2016multiple]. Several outranking methods in MCDM (Multiple-Criteria Decision Making) define a unique ranking mechanism to select the best alternative among others.

Suppose a decision process has $n$ alternatives and $m$ criteria that are either to be maximized or minimized. Each single criterion has a weight $0 \le w_i \le 1$ where $\sum_i^m w_i = 1$ and is represented by a function $f_i$ which is either maximum or minimum. $g_j(.)$ is an evolution function and it is taken as $g_j(x) = x$ in many methods. A multiple criteria decision problem can be represented using the decision table shown in Table \ref{decision_table} without loss of generality. When $A_1$, $A_2$, $\dots$, $A_n$ are alternatives and $C_1$, $C_2$, $\dots$, $C_m$ are different situations of a single criterion then the decision problem is said to be a single criterion decision problem. If $A_i$ and $C_j$ are strategies of two game players then $g_j(A_i)$ is the gain of the row player when she selects the strategy $i$ and the column player selects the strategy $C_j$. 

Table: Decision table \label{decision_table}

| **Criteria**  |   $C_1$    |   $C_2$    | $\dots$  |   $C_m$    |
| :-----------: | :--------: | :--------: | :------: | :--------: |
|  **Weights**  |   $w_1$    |   $w_2$    | $\dots$  |   $w_m$    |
| **Functions** |   $f_1$    |   $f_2$    | $\dots$  |   $f_m$    |
|     $A_1$     | $g_1(A_1)$ | $g_2(A_1)$ | $\dots$  | $g_m(S_A)$ |
|     $A_2$     | $g_1(A_2)$ | $g_2(A_2)$ | $\dots$  | $g_m(A_2)$ |
|       ⋮       |     ⋮      |     ⋮      | $\ddots$  |     ⋮      |
|     $A_n$     | $g_1(A_n)$ | $g_2(A_n)$ | $\dots$  | $g_m(A_n)$ |


# State of the field

Multiple-criteria decision-making (MCDM) tools provide several algorithms for ordering or  selecting alternatives and/or determining the weights when there is uncertainty. Although some algorithms are suitable for hand calculations, computer software is often required. While some previous applications only focused on a single method, some applications appear to include multiple methods. ```PyTOPS``` is a Python tool for TOPSIS [@PyTOPS]. ```Super Decisions``` is a software package that is mainly focused on AHP (Analytic Hierarchy Process) and ANP (Analytic Network Process) [@superdecision]. ```Visual Promethee``` implements the Promethee method on Windows platforms [@visualpromethee]. ```M-BACBETH``` is another commercial software product that implements MACBETH with an easy to use GUI [@macbeth]. ```Sanna``` is a standard ```MS Excel``` add-in application that supports several basic methods for multi-criteria evaluation of alternatives (WSA, TOPSIS, ELECTRE I and III, PROMETHEE I and II, MAPPAC and ORESTE) [@sanna]. The ```DEAFrontier``` software requires an ```Excel``` add-in that can solve up to 50 DMUs with unlimited number of inputs and outputs (subject to the capacity of the standard ```MS Excel Solver```) [@deafrontier]. 



# Statement of need 

While the applications mentioned above are lacking in features such as the number of methods included, being programmable, being free, and the results being comparable by the researcher, ```JMcDM``` clearly differs as it has all of these features.
```JMcDM``` is designed to provide a developer-friendly library for solving multiple-criteria decision problems in ```Julia``` [@julia]. Since ```Julia``` is a dynamic language, it is also useful for researchers who are familiar with REPL (Read-Eval-Print-Loop) environments. The package includes multi-criteria decision methods as well as a game solver for zero-sum games, and methods for single criterion methods. 

The package implements methods for 
AHP [@ahp],
ARAS [@aras],
COCOSO [@cocoso],
CODAS [@codas],
COPRAS [@copras], 
CRITIC [@critic],
DEMATEL [@dematel], 
EDAS [@edas], 
ELECTRE [@electre], 
Entropy [@entropy],
GRA [@gra], 
MABAC [@mabac], 
MAIRCA [@mairca], 
MARCOS [@marcos], 
MOORA [@moora], 
NDS [@Deb_2002], 
PROMETHEE [@promethee], 
SAW [@saw; @wsm_wpm],  
TOPSIS [@topsis],
VIKOR [@vikor_1; @vikor_2],  
WASPAS [@waspas], 
and
WPM [@wsm_wpm]
for multiple-criteria tools. This list of selected methods includes both classical (TOPSIS, ELECTRE, PROMETHEE, etc.) and modern (COCOSO, MABAC, MARCOS, etc.) tools of the relevant literature. 

The package also performs Data Envelopment Analysis (DEA) [@dea] and includes a method for solving zero-sum games. Although these methods may seem different from the methods mentioned above, they are basically members of the same method family and solve similar problems. DEA differs from the above methods in that it is not an outranking method but compares efficiencies of decision units. Solving zero-sum games is also a multi-criteria decision-making problem, but this time, unlike outranking methods, both the rows and columns of the decision matrix show alternative strategies. 

The full set of other tools and utility functions are listed and documented in the source code as well as in the online documentation.

# Installation and basic usage

`JMcDM` can be downloaded and installed using the Julia package manager by typing

```julia
julia> using Pkg
julia> Pkg.add("JMcDM")
```

and can be loaded before using any functions by typing

```julia
julia> using JMcDM
```

in ```Julia``` REPL.

Suppose a decision problem is given in the table below.

| **Criteria**  |  Age   |  Size  |  Price   | Distance | Population |
| :-----------: | :----: | :----: | :------: | :------: | :--------: |
|  **Weights**  | $0.35$ | $0.15$ |  $0.25$  |  $0.20$  |   $0.05$   |
| **Functions** |  min   |  max   |   min    |   min    |    max     |
|     $A_1$     |  $6$   | $140$  | $150000$ |  $950$   |   $1500$   |
|     $A_2$     |  $4$   |  $90$  | $100000$ |  $1500$  |   $2000$   |
|     $A_3$     |  $12$  | $140$  | $75000$  |  $550$   |   $1100$   |

In this sample problem, a decision maker is subject to select an apartment by considering the age of the building, size (in $m^2$s), price (in \$), distance to city centre (in $m$s), and nearby population.
The data can be entered as a two-dimensional array (matrix) or as a DataFrame object:

```julia
julia> using JMcDM
julia> df = DataFrame(
:age        => [6.0, 4, 12],
:size       => [140.0, 90, 140],
:price      => [150000.0, 100000, 75000],
:distance   => [950.0, 1500, 550],
:population => [1500.0, 2000, 1100]);
```
The weight vector ```w```, vector of directions ```fns```, and ```topsis()``` function call can be performed using the ```Julia``` REPL.

```julia
julia> w  = [0.35, 0.15, 0.25, 0.20, 0.05];
julia> fns = makeminmax([minimum, maximum, minimum, minimum, maximum]);
julia> result = topsis(df, w, fns);
julia> result.scores
3-element Array{Float64,1}:
0.5854753145549456
0.6517997936899308
0.41850223305822903

julia> result.bestIndex
2
```

In the output above, it is shown that the alternative $A_2$ has a score of $0.65179$ and it is selected as the best. The same analysis can be performed using ```saw()``` for the method of Simple Additive Weighting

```julia
julia> result = saw(df, w, fns);
julia> result.bestIndex
2
```

as well as using ```wpm``` for the method of Weighted Product Method 

```julia
julia> result = wpm(df, w, fns);
julia> result.bestIndex
2
```

For any method, ```?methodname``` shows the documentation as in the same way in other ```Julia``` packages.

# References
# Contents 

[The GitHub Repository](https://github.com/jbytecode/JMcDM)
```@contents
Pages = ["index.md", "mcdms.md", "game.md", "dataenvelop.md", "utility.md"]
Depth = 3
```
# Data Envelopment Analysis

## dataenvelop
```@docs
JMcDM.dataenvelop
```


# Multiple Criteria Decision Making Tools

## TOPSIS
```@docs
JMcDM.topsis
```


## ELECTRE
```@docs
JMcDM.electre
```


## DEMATEL
```@docs
JMcDM.dematel
```


## MOORA
```@docs
JMcDM.moora
```


## VIKOR
```@docs
JMcDM.vikor
```


## AHP
```@docs
JMcDM.ahp
```


## Grey Relational Analysis
```@docs
JMcDM.grey
```

## Non-dominated Sorting
```@docs
JMcDM.nds
```

## SAW
```@docs
JMcDM.saw
```

## ARAS
```@docs
JMcDM.aras
```

## WPM
```@docs
JMcDM.wpm
```

## WASPAS
```@docs
JMcDM.waspas
```


## EDAS
```@docs
JMcDM.edas
```

## MARCOS
```@docs
JMcDM.marcos
```

## MABAC
```@docs
JMcDM.mabac
```

## MAIRCA
```@docs
JMcDM.mairca
```


## COPRAS
```@docs
JMcDM.copras
```

## PROMETHEE
```@docs
JMcDM.promethee
```


## CoCoSo
```@docs
JMcDM.cocoso
```

## Critic
```@docs
JMcDM.critic
```


## CODAS
```@docs
JMcDM.codas
```





# Single-criterion decision making tools

## laplace
```@docs
JMcDM.laplace
```

## maximin
```@docs
JMcDM.maximin
```

## maximax
```@docs
JMcDM.maximax
```


## minimax
```@docs
JMcDM.minimax
```

## minimin
```@docs
JMcDM.minimin
```

## savage
```@docs
JMcDM.savage
```

## hurwicz
```@docs
JMcDM.hurwicz
```

## mle
```@docs
JMcDM.mle
```


## expectedregret
```@docs
JMcDM.expectedregret
```



# Utility functions

## MCDMSetting
```@docs
JMcDM.MCDMSetting
```

## mcdm
```@docs
JMcDM.mcdm
```

## summary
```@docs
JMcDM.summary
```


## makeDecisionMatrix
```@docs
JMcDM.makeDecisionMatrix
```# Zero-sum Game Solver

## game
```@docs
JMcDM.game
```

