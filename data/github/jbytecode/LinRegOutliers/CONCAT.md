# v0.8.9

- LAD (Least Absolute Deviations) is now exact and uses a linear programming based model
- Dependencies for JuMP and GLPK are added 
- Dependency for Optim removed

[![Build Status](https://travis-ci.org/jbytecode/LinRegOutliers.svg?branch=master)](https://travis-ci.org/jbytecode/LinRegOutliers) [![DOI](https://joss.theoj.org/papers/10.21105/joss.02892/status.svg)](https://doi.org/10.21105/joss.02892)

# LinRegOutliers

A Julia package for outlier detection in linear regression.

## Implemented Methods
- Ordinary Least Squares, Weighted Least Squares, Basic diagnostics
- Hadi & Simonoff (1993)
- Kianifard & Swallow (1989)
- Sebert & Montgomery & Rollier (1998)
- Least Median of Squares 
- Least Trimmed Squares 
- Minimum Volume Ellipsoid (MVE)
- MVE & LTS Plot 
- Billor & Chatterjee & Hadi (2006)
- Pena & Yohai (1995)
- Satman (2013)
- Satman (2015)
- Setan & Halim & Mohd (2000)
- Least Absolute Deviations (LAD)
- Least Trimmed Absolute Deviations (LTA)
- Hadi (1992)
- Marchette & Solka (2003) Data Images
- Satman's GA based LTS estimation (2012)
- Fischler & Bolles (1981) RANSAC Algorithm
- Minimum Covariance Determinant Estimator
- Imon (2005) Algorithm
- Barratt & Angeris & Boyd (2020) CCF algorithm
- Atkinson (1994) Forward Search Algorithm
- BACON Algorithm (Billor & Hadi & Velleman (2000))
- Hadi (1994) Algorithm
- Chatterjee & Mächler (1997)
- Summary


## Unimplemented Methods
- Depth based estimators (Regression depth, deepest regression, etc.)
- Theil & Sen estimator for multiple regression


## Installation

```LinRegOutliers``` can be installed using the ```Julia``` REPL.  

```julia
julia> ]
(@v1.5) pkg> add LinRegOutliers
```

or

```julia
julia> using Pgk
julia> Pkg.add("LinRegOutliers")
```

then

```julia
julia> using LinRegOutliers
```

to make all the stuff be ready!


## Examples
We provide some examples [here](https://github.com/jbytecode/LinRegOutliers/blob/master/examples.md).
 
## Documentation
Please check out the reference manual [here](https://jbytecode.github.io/LinRegOutliers/docs/build/).

## News
- We implemented algorithm(X, y) style calls for all of the algorithms where X is the design matrix and y is the response vector. 
- We implemented ~25 outlier detection algorithms which covers a high percentage of the literature.


## Contributions
You are probably the right contributor

- If you have statistics background
- If you like Julia

However, the second condition is more important because an outlier detection algorithm is just an algorithm. Reading the implemented methods is enough to implement new ones. Please follow the issues. [Here is the a bunch of first shot introductions for new comers](https://github.com/jbytecode/LinRegOutliers/issues/3). Welcome and thank you in advance!


## Citation
Please refer our original paper if you use the package in your research using

```
Satman et al., (2021). LinRegOutliers: A Julia package for detecting outliers in linear regression. Journal of Open Source Software, 6(57), 2892, https://doi.org/10.21105/joss.02892
```

or the bibtex entry

```
@article{Satman2021,
  doi = {10.21105/joss.02892},
  url = {https://doi.org/10.21105/joss.02892},
  year = {2021},
  publisher = {The Open Journal},
  volume = {6},
  number = {57},
  pages = {2892},
  author = {Mehmet Hakan Satman and Shreesh Adiga and Guillermo Angeris and Emre Akadal},
  title = {LinRegOutliers: A Julia package for detecting outliers in linear regression},
  journal = {Journal of Open Source Software}
}
```


## Contact & Communication
- Please use issues for a new feature request or bug reports.
- We are in #linregoutliers channel on [Julia Slack](http://julialang.slack.com/) for any discussion requires online chatting. 
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
## welcome, contributer!

Please read the [pull requests](https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/about-pull-requests) 
section in GitHub before preparing any feauture for this library.
# Some examples


## A brief introduction

Suppose the linear regression model is 

*y = Xβ + ε*

where *y* is the vector of dependent variable, *X* is the *n ˟ p* matrix of design, *ε* is i.i.d error term with zero mean, *n* is the number of observations, and *p* is the number of regression parameters. 

When a single regressor exists in the model, it can be basically written as 

*y = β₀ + β₁ x +  ε*

where *β₀* and *β₁* are unknown intercept and slope parameters. In `R` and `Julia` we can represent this model in a similar form. Specifically, in Julia, the simple model can be expressed using the `@formula` macro as

```julia
@formula(y ~ x)
```

where ```~``` operator seperates the dependent and independent variables. When the model includes more than one regressors, the model can similarly be expressed as

```julia
@formula(y ~ x1 + x2 + x3)
```

```LinRegOutliers``` follows this convention for expressing linear models. 

_________________

## Sebert & Montgomery & Rollier (1998) Algorithm
Sebert & Montgometry & Rollier (smr98) algorithm starts with an ordinary least squares estimation for a given model and data. Residuals and fitted responses are calculated using the estimated model. A hierarchical clustering analysis is applied using standardized residuals and standardized fitted responses. The tree structure of clusters are cut using a threshold, e.g Majona criterion, as suggested by the authors. It is expected that the subtrees with relatively small number of observations are declared to be clusters of outliers.

Hawkings & Bradu & Kass dataset has 4 variables and 75 observations. The observations 1-14 are known to be outliers. In the example below, we create an regression setting using the formula ```y ~ x1 + x2 + x3``` and ```hbk``` dataset. ```smr98``` is directly applied on this setting.  

```julia
julia> using LinRegOutliers
julia> # Regression setting for Hawkins & Bradu & Kass data
julia> reg = createRegressionSetting(@formula(y ~ x1 + x2 + x3), hbk)
julia> smr98(reg)
Dict{String,Array{Int64,1}} with 1 entry:
  "outliers" => [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]
```

The Julia method ```smr98()``` returns a ```Dict``` object with produced output. In this case, the single output is indices of detected outlying observations. 
The algorithm successfully detects the outliers.

________________________
## Peña and Yohai (1995) 

Peña and Yohai (```py1995```) algorithm starts by constructing an influence matrix using results of an ordinary least squares estimate for a given model and data. In the second stage, the eigen structure of the influence matrix is examined for detecting subset of potential outliers of data. 

Here is an example of ```py95``` method applied on the ```hbk``` data. The method returns a ```Dict```  object with keys ```outliers``` and ```suspected.sets```. An usual researcher may directly focus on the ```outliers``` indices. The method reports the observations 1-14 are outliers.

```julia
julia> py95(reg)
Dict{Any,Any} with 2 entries:
  "outliers"       => [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]
  "suspected.sets" => Set([[14, 13], [43, 54, 24, 38, 22], Int64[], [58, 66, 32, 28, 65, 36], [62], [73, 2…
```


________________________________

## Least Trimmed Squares Regression
Least Trimmed Squares (LTS) is a robust regression estimator with high break-down point. LTS searches for the parameter estimates that minimize sum of the *h* smallest squared residuals where *h* is a constant larger than *n/2*.    

Phone data is a regression data with a single regressor variable. The independent and dependent variables are ```year``` and ```calls``` and have 24 observations. Observations 14-21 are said to be outliers. 

Since LTS is a robust method, the parameter estimates are of interest. However, we provide indices of outliers in results for diagnostic purposes only. 

The method ```lts``` also reports standard deviation of estimate, scaled residuals and LTS objective function as well. 

```julia
julia> reg = createRegressionSetting(@formula(calls ~ year), phones);

julia> lts(reg)
Dict{Any,Any} with 6 entries:
  "betas"            => [-56.5219, 1.16488]
  "S"                => 1.10918
  "hsubset"          => [11, 10, 5, 6, 23, 12, 13, 9, 24, 7, 3, 4, 8]
  "outliers"         => [14, 15, 16, 17, 18, 19, 20, 21]
  "scaled.residuals" => [2.41447, 1.63472, 0.584504, 0.61617, 0.197052, -0.222066, -0.551027, -0.970146, -0.397538, -0.185558  …  91.0312, 94.4889, 109.667, 123.943, 143.629, …
  "objective"        => 3.43133
```  



```julia
using Plots
x = phones[:,"year"]
y = phones[:,"calls"]
f(x) = -56.5219 +  1.16488x 
scatter(x, y, label=false, title="Phone Data")
px = [x[1], x[end]]
py = map(f, px)
plot!(px, py, label=false, color=:red, width=2) 
```


<img src="https://github.com/jbytecode/jbytecode/blob/master/images/ltsandphonedata.png" alt="dataimages" width="500"/>

Figure 1 - Phone Data and estimated LTS line



_________________

## Data Images

The method ```dataimage``` implements the Data Image algorithm and serves a visual tool as an outlier detection algorithm for multivariate data only. The algorithm generates a color matrix with each single cell represents a proper distance between observations. Since 

```dataimage(data, distance = :euclidean)```

defines color using the Euclidean distance, whereas

```dataimage(data, distance = :mahalabobis)```

uses Mahalanobis distances for determining color values. The default distance metric is Euclidean distance. 

In the example below, the distances between observations are calculated and drawn using corresponding colors. Since the method is for multivariate data, only the desing matrix is used. In other terms, the response vector is omitted. 

```julia
julia> # Matrix of independent variables of Hawkins & Bradu & Kass data
julia> data = hcat(hbk.x1, hbk.x2, hbk.x3);
julia> dataimage(data)
``` 


<img src="https://github.com/jbytecode/jbytecode/blob/master/images/dataimages.png" alt="dataimages" width="500"/>

Figure 2 - Data Image of Design Matrix of ```hbk``` Data


_________________________
## Atkinson's Stalactite Plot

Atkinson's Stalactite Plot serves a visual method for detecting outliers in linear regression. Despite it shares the same calling convention with the other methods, the method ```atkinsonstalactiteplot``` generates a text based plot. The method performs a robust regression estimator many times and residuals higher than some threshold are labelled using ```+``` and ```*```. After many iterations, the observations with many labels are considered as suspected or outlying observations.  

```julia
julia> using LinRegOutliers
julia> reg = createRegressionSetting(@formula(calls ~ year), phones);
julia> atkinsonstalactiteplot(reg)
m           1         2
   123456789012345678901234
 2              ********   
 3              ********   
 4 +            ********   
 5 +            ********   
 6 +            ********   
 7 +            ********   
 8 +            ********   
 9 +            ********+  
10 +            ********+  
11 +            ********   
12              ********   
13 +            ********   
14              ********   
15              ********   
16              ********   
17              ********   
18               *******+  
19               ****** +++
20               ****** +++
21               ****** +++
22               ****** +++
23               +++*** +++
24                  ++* +++
   123456789012345678901234
            1         2

```

The output above can be considered as an evidence that the observations 14-21 are suspected. Observations 1, 22, 23, 24 are also labelled as ```+``` in some iterations. However, the frequency of labels of these observations are relatively small. 


____________________
## Other algorithms
```LinRegOutliers``` implements more than 20 outlier detection methods in linear regression and covers a big proportion of the classical literature in this subject. The documentation of the package includes the referenced citations. Any researcher can follow the details of algorithms using these information. 

_____________________
## Other calling conventions
The calling convention 

```julia
julia> setting = createRegressionSetting(@formula(...), data)
julia> method(setting) 
```

is the preferred way of calling implemented methods in ```LinRegOutliers```, we multiple dispatch the methods using the syntax

```julia
julia> method(X, y) 
```

where *X* is the design matrix and *y* is the response vector. This calling convention may be more suitable for those who iteratively calls the methods possibly in a simulation study or other kinds of researching stuff. 

For example, we can perform ```hs93``` on the Phones data using 

```julia
julia> hs93(reg)
Dict{Any,Any} with 3 entries:
  "outliers" => [14, 15, 16, 17, 18, 19, 20, 21]
  "t"        => -3.59263
  "d"        => [2.04474, 1.14495, -0.0633255, 0.0632934, -0.354349, -0.766818, -1.06862, -1.47638, -0.710…
```

as well as 

```julia
julia> X = hcat(ones(24), phones[:, "year"]);

julia> y = phones[:, "calls"];

julia> hs93(X, y)
Dict{Any,Any} with 3 entries:
  "outliers" => [14, 15, 16, 17, 18, 19, 20, 21]
  "t"        => -3.59263
  "d"        => [2.04474, 1.14495, -0.0633255, 0.0632934, -0.354349, -0.766818, -1.06862, -1.47638, -0.710…
```

__________________________
## Multiple Methods in a single shot!

We also provide ```detectOutliers``` method for data scientist for performing many methods and presenting the summarized results. 

The method can be called using default arguments only by feeding a regression setting object:

```julia
julia> detectOutliers(aSettingObject)
```

The method generates a console output:

<img src="https://github.com/jbytecode/jbytecode/blob/master/images/detectoutliers.png" alt="dataimages" width="500"/>
---
title: 'LinRegOutliers: A Julia package for detecting outliers in linear regression'
tags:
  - Julia
  - linear regression
  - outlier detection
  - robust statistics
authors:
  - name: Mehmet Hakan Satman
    orcid: 0000-0002-9402-1982
    affiliation: 1
  - name: Shreesh Adiga
    orcid: 0000-0002-1818-6961
    affiliation: 2
  - name: Guillermo Angeris
    orcid: 0000-0002-4950-3990
    affiliation: 3
  - name: Emre Akadal
    orcid: 0000-0001-6817-0127 
    affiliation: 4
affiliations:
 - name: Department of Econometrics, Istanbul University, Istanbul, Turkey
   index: 1
 - name: Department of Electronics and Communication Engineering, RV College of Engineering, Bengaluru, India
   index: 2
 - name: Department of Electrical Engineering, Stanford University, Stanford, California, USA
   index: 3
 - name: Department of Informatics, Istanbul University, Istanbul, Turkey
   index: 4

date: 26 November 2020
bibliography: paper.bib
---

# Summary

`LinRegOutliers` is a Julia package that implements a number of outlier detection algorithms for linear regression. The package also implements robust covariance matrix estimation and graphing functions which can be used to visualize the regression residuals and distances between observations, with many possible metrics (*e.g.*, the Euclidean or Mahalanobis distances with either given or estimated covariance matrices). Our package implements many algorithms and diagnostics for model fitting with outliers under a single interface, which allows users to quickly try many different methods with reasonable default settings, while also providing a good starting framework for researchers who may want to extend the package with novel methods.


# State of the field
In linear regression, we are given a number of data points (say, $n$) where each data point is represented by a vector $x_i$, with $p$ entries, and a dependent variable that corresponds to each of these data points, represented by the scalar $y_i$, for $i=1, 2, \dots, n$. We then seek to find a linear model which best describes the data (up to some error term, $\epsilon_i$):

$$
y_i = \beta_1 (x_{i})_1+ \dots + \beta_{p} (x_i)_p +  \epsilon_i,
$$

where $\beta_1, \dots, \beta_p$ are the $p$ unknown parameters. We will assume that the $\epsilon_i$ are independent and identically-distributed (i.i.d.) error terms with zero mean. Note that, if $(x_i)_1 = 1$ for all $i=1, \dots, n$, this is equivalent to having an intercept term given by $\beta_1$.

We can write this more conveniently by letting $X$ be the *design matrix* of size $n\times p$, whose $i$th row is given by the vectors $x_i$ (where $(x_i)_1=1$ if the model has an intercept), while $y$ is an $n$-vector of observations, whose entries are $y_i$, and similarly for $\epsilon$:
$$
y = X\beta + \epsilon.
$$
The usual approach to finding an estimate for $\beta$, which we call $\hat \beta$, is the Ordinary Least Squares (OLS) estimator given by $\hat{\beta} = (X^TX)^{-1}X^Ty$, which is efficient and has good statistical properties when the error terms are all of roughly the same magnitude (*i.e.*, there are no outliers). On the other hand, the OLS estimator is very sensitive to outliers: even if a single  observation lies far from the regression hyperplane, OLS will often fail to find a good estimate for the parameters, $\beta$.

To solve this problem, a number of methods have been developed in the literature. These methods can be roughly placed in one or more of the five following categories: diagnostics, direct methods, robust methods, multivariate methods, and visual methods. *Diagnostics* are methods which attempt to find points that significantly affect the fit of a model (often, such points can be labeled as outliers). Diagnostics can then be used to initialize *direct methods*, which fit a (usually non-robust) model to a subset of points suspected to be clear of outliers; remaining points which are not outliers with respect to this fit are continually added to this subset until all points not in the subset are deemed outliers. *Robust methods*, on the other hand, find a best-fit model by approximately minimizing a loss function that is resistant to outliers. Some of the proposed methods are also *multivariate methods*, which can accommodate obtaining robust location and scale measures of multivariate data. *Visual methods* generally work on the principle of visualizing the statistics obtained from these mentioned methods. As an example, the method `mveltsplot` constructs a 2-D plot using robust distances and scaled residuals obtained from `mve` and `lts` which are multivariate data and robust regression methods, respectively. Many direct and robust methods for regression select an initial basic or clean subset of observations using the results of diagnostics and methods for multivariate data. This is why methods that are not directly related to regression are included in the package. 

# Statement of need 

In practice, many of the proposed methods have reasonable performance and yield similar results for most datasets, but sometimes differ widely in specific circumstances by means of masking and swamping ratios. Additionally, some of the methods are relatively complicated and, if canonical implementations are available, they are often out of date or only found in specific languages of the author's choice, making it difficult for researchers to compare the performance of these algorithms on their datasets.

We have reimplemented many of the algorithms available in the literature in Julia [@julia], an open-source, high performance programming language designed primarily for scientific computing. Our package, `LinRegOutliers`, is a comprehensive and simple-to-use Julia package that includes many of the algorithms in the literature for detecting outliers in linear regression. The implemented `Julia` methods for diagnostics, direct methods, robust methods, multivariate methods, and visual diagnostics are shown in **Table 1**, **Table 2**, **Table 3**, **Table 4**, and **Table 5**, respectively. 
 

| Algorithm             | Reference      | Method                     |
| :-------------------- | :------------- | :------------------------- |
| Hadi Measure          | [@hadimeasure] | `hadimeasure`              |
| Covariance Ratio      | [@diagnostics] | `covratio`                 |
| DFBETA                | [@diagnostics] | `dfbeta`                   |
| DFFIT                 | [@diagnostics] | `dffit`                    |
| Mahalanobis Distances | [@mahalanobis] | `mahalanobisSquaredMatrix` |
| Cook Distances        | [@cooks]       | `cooks`                    |

Table: Regression Diagnostics

| Algorithm   | Reference     | Method       |
| :---------- | :------------ | :----------- |
| Ransac      | [@ransac]     | `ransac`     |
| KS-89       | [@ks89]       | `ks89`       |
| HS-93       | [@hs93]       | `hs93`       |
| Atkinson-94 | [@atkinson94] | `atkinson94` |
| PY-95       | [@py95]       | `py95`       |
| SMR-98      | [@smr98]      | `smr98`      |
| ASM-2000    | [@asm2000]    | `asm2000`    |
| BACON       | [@bacon]      | `bacon`      |
| Imon-2005   | [@imon2005]   | `imon2005`   |
| bch         | [@bch]        | `bch`        |

Table: Direct Methods


| Algorithm                         | Reference     | Method       |
| :-------------------------------- | :------------ | :----------- |
| Least Absolute Deviations         | [@lad]        | `lad`        |
| Least Absolute Trimmed Deviations | [@lta]        | `lta`        |
| Least Median of Squares           | [@lms]        | `lms`        |
| Least Trimmed Squares             | [@lts]        | `lts`        |
| CM-97                             | [@cm97]       | `cm97`       |
| ga-lts                            | [@galts]      | `galts`      |
| Satman-2013                       | [@satman2013] | `satman2013` |
| Satman-2015                       | [@satman2015] | `satman2015` |
| CCF                               | [@ccf]        | `ccf`        |

Table: Robust Methods


| Algorithm                      | Reference   | Method     |
| :----------------------------- | :---------- | :--------- |
| Hadi-1992                      | [@hadi1992] | `hadi1992` |
| Hadi-1994                      | [@hadi1994] | `hadi1994` |
| Minimum Volume Ellipsoid       | [@mve]      | `mve`      |
| Minimum Covariance Determinant | [@mcd]      | `mcd`      |

Table: Multivariate Methods


| Algorithm       | Reference     | Method                   |
| :-------------- | :------------ | :----------------------- |
| BCH Plot        | [@bch]        | `bchplot`                |
| MVE-LTS Plot    | [@mve]        | `mveltsplot`             |
| Data Images     | [@dataimage]  | `dataimage`              |
| Stalactite Plot | [@atkinson94] | `atkinsonstalactiteplot` |

Table: Visual Methods

# Installation and basic usage

`LinRegOutliers` can be downloaded and installed using the Julia package manager by typing

```julia
julia> using Pkg
julia> Pkg.add("LinRegOutliers")
```

in the Julia console. The regression methods follow a uniform call convention. For instance, a user can type

```julia
julia> setting = createRegressionSetting(@formula(calls ~ year), phones);
julia> smr98(setting)
Dict{String,Array{Int64,1}} with 1 entry:
  "outliers" => [15, 16, 17, 18, 19, 20, 21, 22, 23, 24]
```

or

```julia
julia> X = hcat(ones(24), phones[:, "year"]);
julia> y = phones[:, "calls"];
julia> smr98(X, y)
Dict{String,Array{Int64,1}} with 1 entry:
  "outliers" => [15, 16, 17, 18, 19, 20, 21, 22, 23, 24]
```

to apply *smr98* [@smr98] on the Telephone dataset [@lms], where $X$ is the design matrix with ones in its first column. In this case, observations 15 to 24 are reported as outliers by the method. Some methods may also return additional information specific to the method which is passed back in a ```Dict``` object. For example, the *ccf* function returns a ```Dict``` object containing *betas*, *outliers*, *lambdas*, and *residuals*:

```julia
julia> ccf(X, y)
Dict{Any,Any} with 4 entries:
  "betas"     => [-63.4816, 1.30406]
  "outliers"  => [15, 16, 17, 18, 19, 20]
  "lambdas"   => [1.0, 1.0, 1.0, 1.0, 1.0, ...
  "residuals" => [-2.67878, -1.67473, -0.37067, -0.266613, …
```

Indices of outliers can be accessed using standard ```Dict``` operations like

```julia
julia> result = ccf(X, y)
julia> result["outliers"]
6-element Array{Int64,1}:
 15
 16
 17
 18
 19
 20
```
 
# Acknowledgements

Guillermo Angeris is supported by the National Science Foundation Graduate Research Fellowship under Grant No. DGE-1656518. 




# References
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: ''
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

# Contents 
```@contents
Pages = ["index.md", "datasets.md", "types.md", "diagnostics.md", "algorithms.md"]
Depth = 3
```





# Algorithms

## Hadi & Simonoff (1993)
```@docs
LinRegOutliers.hs93
```

## Kianifard & Swallow (1989)
```@docs
LinRegOutliers.ks89
```

## Sebert & Montgomery & Rollier (1998)
```@docs
LinRegOutliers.smr98
```

## Least Median of Squares
```@docs
LinRegOutliers.lms
```

## Least Trimmed Squares
```@docs
LinRegOutliers.lts
```

## Minimum Volume Ellipsoid (MVE)
```@docs
LinRegOutliers.mve
```

## MVE & LTS Plot
```@docs
LinRegOutliers.mveltsplot
```

## Billor & Chatterjee & Hadi (2006)
```@docs
LinRegOutliers.bch
```

## Pena & Yohai (1995)
```@docs
LinRegOutliers.py95
```

## Satman (2013)
```@docs
LinRegOutliers.satman2013
```

## Satman (2015)
```@docs
LinRegOutliers.satman2015
```

## Setan & Halim & Mohd (2000)
```@docs
LinRegOutliers.asm2000
```

## Least Absolute Deviations (LAD)
```@docs
LinRegOutliers.lad
```

## Least Trimmed Absolute Deviations (LTA)
```@docs
LinRegOutliers.lta
```

## Hadi (1992)
```@docs
LinRegOutliers.hadi1992
```

## Marchette & Solka (2003) Data Images
```@docs
LinRegOutliers.dataimage
```

## Satman's GA based LTS estimation (2012)
```@docs
LinRegOutliers.galts
```

## Fischler & Bolles (1981) RANSAC Algorithm
```@docs
LinRegOutliers.ransac
```

## Minimum Covariance Determinant Estimator (MCD)
```@docs
LinRegOutliers.mcd
```

## Imon (2005) Algorithm
```@docs
LinRegOutliers.imon2005
```

## Barratt & Angeris & Boyd (2020) CCF algorithm
```@docs
LinRegOutliers.ccf
```

## Atkinson (1994) Forward Search Algorithm
```@docs
LinRegOutliers.atkinson94
```

## BACON Algorithm (Billor & Hadi & Velleman (2000))
```@docs
LinRegOutliers.bacon
```

## Hadi (1994) Algorithm
```@docs
LinRegOutliers.hadi1994
```

## Chatterjee & Mächler (1997)
```@docs
LinRegOutliers.cm97
```















# Types

## RegressionSetting
```@docs
LinRegOutliers.RegressionSetting
```

## OLS
```@docs
LinRegOutliers.OLS
```
# Datasets

## Phone data
```@docs
LinRegOutliers.phones
```


## Hawkings & Bradu & Kass data
```@docs
LinRegOutliers.hbk
```

## Animals data
```@docs
LinRegOutliers.animals
```

## Weight Loss data
```@docs
LinRegOutliers.weightloss
```


## Stack Loss data
```@docs
LinRegOutliers.stackloss
```

## Hadi & Simonoff (1993) random data
```@docs
LinRegOutliers.hs93randomdata
```

## Modified Wood Gravity data
```@docs
LinRegOutliers.woodgravity
```

## Scottish Hill Races data
```@docs
LinRegOutliers.hills
```

## Soft Drink Delivery data
```@docs
LinRegOutliers.softdrinkdelivery
```

# Diagnostics

## ols 
```@docs
LinRegOutliers.ols
```

## wls
```@docs
LinRegOutliers.wls
```


## dffit
```@docs
LinRegOutliers.dffit
```

## hatmatrix
```@docs
LinRegOutliers.hatmatrix
```

## studentizedResiduals
```@docs
LinRegOutliers.studentizedResiduals
```

## adjustedResiduals
```@docs
LinRegOutliers.adjustedResiduals
```

## jacknifedS
```@docs
LinRegOutliers.jacknifedS
```

## cooks
```@docs
LinRegOutliers.cooks
```

## mahalanobisSquaredMatrix
```@docs
LinRegOutliers.mahalanobisSquaredMatrix
```

## dfbeta
```@docs
LinRegOutliers.dfbeta
```

## covratio
```@docs
LinRegOutliers.covratio
```

## hadimeasure
```@docs
LinRegOutliers.hadimeasure
```

