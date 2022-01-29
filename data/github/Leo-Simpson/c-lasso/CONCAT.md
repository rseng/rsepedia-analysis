[![arXiv](https://img.shields.io/badge/arXiv-2011.00898-b31b1b.svg)](https://arxiv.org/abs/2011.00898)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02844/status.svg)](https://doi.org/10.21105/joss.02844)

<img src="https://i.imgur.com/2nGwlux.png" alt="c-lasso" height="145" align="right"/>

# c-lasso: a Python package for constrained sparse regression and classification 


c-lasso is a Python package that enables sparse and robust linear regression and classification with linear equality
constraints on the model parameters. For detailed info, one can check the [documentation](https://c-lasso.readthedocs.io/en/latest/).

The forward model is assumed to be: 

<img src="https://latex.codecogs.com/gif.latex?y=X\beta&plus;\sigma\epsilon\qquad\text{s.t.}\qquad&space;C\beta=0" title="y=X\beta+\sigma\epsilon\qquad\text{s.t.}\qquad C\beta=0" />

Here, y and X are given outcome and predictor data. The vector y can be continuous (for regression) or binary (for classification). C is a general constraint matrix. The vector &beta; comprises the unknown coefficients and &sigma; an 
unknown scale.

The package handles several different estimators for inferring &beta; (and &sigma;), including 
the constrained Lasso, the constrained scaled Lasso, sparse Huber M-estimation with linear equality constraints, and regularized Support Vector Machines.
Several different algorithmic strategies, including path and proximal splitting algorithms, are implemented to solve 
the underlying convex optimization problems.

We also include two model selection strategies for determining the sparsity of the model parameters: k-fold cross-validation and stability selection.   

This package is intended to fill the gap between popular python tools such as [scikit-learn](https://scikit-learn.org/stable/) which CANNOT solve sparse constrained problems and general-purpose optimization solvers that do not scale well or are inaccurate (see [benchmarks](./benchmark/README.md)) for the considered problems. In its current stage, however, c-lasso is not yet compatible with the scikit-learn API but rather a stand-alone tool.

Below we show several use cases of the package, including an application of sparse *log-contrast*
regression tasks for *compositional* microbiome data.

The code builds on results from several papers which can be found in the [References](#references). We also refer to the accompanying [JOSS paper submission](https://github.com/Leo-Simpson/c-lasso/blob/master/paper/paper.md), also available on [arXiv](https://arxiv.org/pdf/2011.00898.pdf).

## Table of Contents

* [Installation](#installation)
* [Regression and classification problems](#regression-and-classification-problems)
* [Getting started](#getting-started)
* [Log-contrast regression for microbiome data](#log-contrast-regression-for-microbiome-data)
* [Optimization schemes](#optimization-schemes)


* [References](#references)


##  Installation

c-lasso is available on pip. You can install the package
in the shell using

```shell
pip install c-lasso
```
To use the c-lasso package in Python, type 

```python

from classo import classo_problem 
# one can add auxiliary functions as well such as random_data or csv_to_np
```

The `c-lasso` package depends on the following Python packages:

- `numpy`; 
- `matplotlib`; 
- `scipy`; 
- `pandas`; 
- `pytest` (for tests)

##  Regression and classification problems

The c-lasso package can solve six different types of estimation problems: 
four regression-type and two classification-type formulations.

#### [R1] Standard constrained Lasso regression:             

<img src="https://latex.codecogs.com/gif.latex?\arg\min_{\beta\in&space;R^d}&space;||&space;X\beta-y&space;||^2&space;&plus;&space;\lambda&space;||\beta||_1&space;\qquad\mbox{s.t.}\qquad&space;C\beta=0" />

This is the standard Lasso problem with linear equality constraints on the &beta; vector. 
The objective function combines Least-Squares for model fitting with l1 penalty for sparsity.   

#### [R2] Constrained sparse Huber regression:                   

<img src="https://latex.codecogs.com/gif.latex?\arg\min_{\beta\in&space;R^d}&space;h_{\rho}(X\beta-y&space;)&space;&plus;&space;\lambda&space;||\beta||_1&space;\qquad\mbox{s.t.}\qquad&space;C\beta=0" />

This regression problem uses the [Huber loss](https://en.wikipedia.org/wiki/Huber_loss) as objective function 
for robust model fitting with l1 and linear equality constraints on the &beta; vector. The parameter &rho;=1.345.

#### [R3] Constrained scaled Lasso regression: 

<img src="https://latex.codecogs.com/gif.latex?\arg&space;\min_{\beta&space;\in&space;\mathbb{R}^d,&space;\sigma&space;>&space;0}&space;\frac{||&space;X\beta&space;-&space;y||^2}{\sigma}&space;&plus;&space;\frac{n}{2}&space;\sigma&plus;&space;\lambda&space;||\beta||_1&space;\qquad&space;\mbox{s.t.}&space;\qquad&space;C\beta&space;=&space;0" title="\arg \min_{\beta \in \mathbb{R}^d, \sigma > 0} \frac{|| X\beta - y||^2}{\sigma} + \frac{n}{2} \sigma+ \lambda ||\beta||_1 \qquad \mbox{s.t.} \qquad C\beta = 0" />

This formulation is similar to [R1] but allows for joint estimation of the (constrained) &beta; vector and 
the standard deviation &sigma; in a concomitant fashion (see [References](#references) [4,5] for further info).
This is the default problem formulation in c-lasso.

#### [R4] Constrained sparse Huber regression with concomitant scale estimation:        

<img src="https://latex.codecogs.com/gif.latex?\arg&space;\min_{\beta&space;\in&space;\mathbb{R}^d,&space;\sigma&space;>&space;0}&space;\left(&space;h_{\rho}&space;\left(&space;\frac{&space;X\beta&space;-&space;y}{\sigma}&space;\right)&plus;&space;n&space;\right)&space;\sigma&plus;&space;\lambda&space;||\beta||_1&space;\qquad&space;\mbox{s.t.}&space;\qquad&space;C\beta&space;=&space;0" title="\arg \min_{\beta \in \mathbb{R}^d, \sigma > 0} \left( h_{\rho} \left( \frac{ X\beta - y}{\sigma} \right)+ n \right) \sigma+ \lambda ||\beta||_1 \qquad \mbox{s.t.} \qquad C\beta = 0" />

This formulation combines [R2] and [R3] to allow robust joint estimation of the (constrained) &beta; vector and 
the scale &sigma; in a concomitant fashion (see [References](#references) [4,5] for further info).

#### [C1] Constrained sparse classification with Square Hinge loss: 

<img src="https://latex.codecogs.com/gif.latex?\arg&space;\min_{\beta&space;\in&space;\mathbb{R}^d}&space;\sum_{i=1}^n&space;l(y_i&space;x_i^\top&space;\beta)&space;&plus;&space;\lambda&space;\left\lVert&space;\beta\right\rVert_1&space;\qquad&space;s.t.&space;\qquad&space;C\beta&space;=&space;0" title="\arg \min_{\beta \in \mathbb{R}^d} \sum_{i=1}^n l(y_i x_i \beta) + \lambda \left\lVert \beta\right\rVert_1 \qquad s.t. \qquad C\beta = 0" />

where the x<sub>i</sub> are the rows of X and l is defined as:

<img src="https://latex.codecogs.com/gif.latex?l(r)&space;=&space;\begin{cases}&space;(1-r)^2&space;&&space;if&space;\quad&space;r&space;\leq&space;1&space;\\&space;0&space;&if&space;\quad&space;r&space;\geq&space;1&space;\end{cases}" title="l(r) = \begin{cases} (1-r)^2 & if \quad r \leq 1 \\ 0 &if \quad r \geq 1 \end{cases}" />

This formulation is similar to [R1] but adapted for classification tasks using the Square Hinge loss
with (constrained) sparse &beta; vector estimation.

#### [C2] Constrained sparse classification with Huberized Square Hinge loss:        

<img src="https://latex.codecogs.com/gif.latex?\arg&space;\min_{\beta&space;\in&space;\mathbb{R}^d}&space;\sum_{i=1}^n&space;l_{\rho}(y_i&space;x_i^\top\beta)&space;&plus;&space;\lambda&space;\left\lVert&space;\beta\right\rVert_1&space;\qquad&space;s.t.&space;\qquad&space;C\beta&space;=&space;0" title="\arg \min_{\beta \in \mathbb{R}^d} \sum_{i=1}^n l_{\rho}(y_i x_i\beta) + \lambda \left\lVert \beta\right\rVert_1 \qquad s.t. \qquad C\beta = 0" />

where the x<sub>i</sub> are the rows of X and l<sub>ρ</sub> is defined as:

<img src="https://latex.codecogs.com/gif.latex?l_{\rho}(r)&space;=&space;\begin{cases}&space;(1-r)^2&space;&if&space;\quad&space;\rho&space;\leq&space;r&space;\leq&space;1&space;\\&space;(1-\rho)(1&plus;\rho-2r)&space;&&space;if&space;\quad&space;r&space;\leq&space;\rho&space;\\&space;0&space;&if&space;\quad&space;r&space;\geq&space;1&space;\end{cases}" title="l_{\rho}(r) = \begin{cases} (1-r)^2 &if \quad \rho \leq r \leq 1 \\ (1-\rho)(1+\rho-2r) & if \quad r \leq \rho \\ 0 &if \quad r \geq 1 \end{cases}" />


This formulation is similar to [C1] but uses the Huberized Square Hinge loss for robust classification 
with (constrained) sparse &beta; vector estimation.


## Getting started

#### Basic example

We begin with a basic example that shows how to run c-lasso on synthetic data. This example and the next one can be found on the notebook 'Synthetic data Notebook.ipynb'

The c-lasso package includes
the routine ```random_data``` that allows you to generate problem instances using normally distributed data.

```python
m, d, d_nonzero, k, sigma = 100, 200, 5, 1, 0.5
(X, C, y), sol = random_data(m, d, d_nonzero, k, sigma, zerosum=True, seed=1)
```
This code snippet generates a problem instance with sparse &beta; in dimension
d=100 (sparsity d_nonzero=5). The design matrix X comprises n=100 samples generated from an i.i.d standard normal
distribution. The dimension of the constraint matrix C is d x k matrix. The noise level is &sigma;=0.5. 
The input ```zerosum=True``` implies that C is the all-ones vector and C&beta;=0. The n-dimensional outcome vector y
and the regression vector &beta; is then generated to satisfy the given constraints. 

Next we can define a default c-lasso problem instance with the generated data:
```python
problem = classo_problem(X, y, C) 
```
You can look at the generated problem instance by typing:

```python
print(problem)
```

This gives you a summary of the form:

```
FORMULATION: R3
 
MODEL SELECTION COMPUTED:  
     Stability selection
 
STABILITY SELECTION PARAMETERS: 
     numerical_method : not specified
     method : first
     B = 50
     q = 10
     percent_nS = 0.5
     threshold = 0.7
     lamin = 0.01
     Nlam = 50
```
As we have not specified any problem, algorithm, or model selection settings, this problem instance
represents the *default* settings for a c-lasso instance: 
- The problem is of regression type and uses formulation [R3], i.e. with concomitant scale estimation. 
- The *default* optimization scheme is the path algorithm (see [Optimization schemes](#optimization-schemes) for further info). 
- For model selection, stability selection at a theoretically derived &lambda; value is used (see [Reference](#references) [4] for details). Stability selection comprises a relatively large number of parameters. For a description of the settings, we refer to the more advanced examples below and the API.

You can solve the corresponding c-lasso problem instance using

```python
problem.solve()
```

After completion, the results of the optimization and model selection routines 
can be visualized using

```python
print(problem.solution)
```

The command shows the running time(s) for the c-lasso problem instance, and the selected variables for sability selection

```
STABILITY SELECTION : 
   Selected variables :  7    63    148    164    168    
   Running time :  1.546s

```

Here, we only used stability selection as *default* model selection strategy. 
The command also allows you to inspect the computed stability profile for all variables 
at the theoretical &lambda; 

![1.StabSel](https://github.com/Leo-Simpson/c-lasso/blob/master/figures/basic/StabSel.png)


The refitted &beta; values on the selected support are also displayed in the next plot

![beta](https://github.com/Leo-Simpson/c-lasso/blob/master/figures/basic/beta.png)


#### Advanced example             

In the next example, we show how one can specify different aspects of the problem 
formulation and model selection strategy.

```python
m,  d,  d_nonzero,  k, sigma = 100, 200, 5, 0, 0.5
(X, C, y), sol = random_data(m, d, d_nonzero, k, sigma, zerosum = True, seed = 4)
problem                                     = classo_problem(X, y, C)
problem.formulation.huber                   = True
problem.formulation.concomitant             = False
problem.model_selection.CV                  = True
problem.model_selection.LAMfixed            = True
problem.model_selection.PATH                = True
problem.model_selection.StabSelparameters.method = 'max'
problem.model_selection.CVparameters.seed = 1
problem.model_selection.LAMfixedparameters.rescaled_lam = True
problem.model_selection.LAMfixedparameters.lam = .1

problem.solve()
print(problem)

print(problem.solution)

```

Results : 
```
     FORMULATION: R2
     
     MODEL SELECTION COMPUTED:  
          Lambda fixed
          Path
          Cross Validation
          Stability selection
     
     LAMBDA FIXED PARAMETERS: 
          numerical_method = Path-Alg
          rescaled lam : True
          threshold = 0.09
          lam = 0.1
          theoretical_lam = 0.224
     
     PATH PARAMETERS: 
          numerical_method : Path-Alg
          lamin = 0.001
          Nlam = 80
     
     
     CROSS VALIDATION PARAMETERS: 
          numerical_method : Path-Alg
          one-SE method : True
          Nsubset = 5
          lamin = 0.001
          Nlam = 80
     
     
     STABILITY SELECTION PARAMETERS: 
          numerical_method : Path-Alg
          method : max
          B = 50
          q = 10
          percent_nS = 0.5
          threshold = 0.7
          lamin = 0.01
          Nlam = 50

     LAMBDA FIXED : 
     Selected variables :  17    59    123    
     Running time :  0.104s

     PATH COMPUTATION : 
     Running time :  0.638s

     CROSS VALIDATION : 
     Selected variables :  16    17    57    59    64    73    74    76    93    115    123    134    137    181    
     Running time :  2.1s

     STABILITY SELECTION : 
     Selected variables :  17    59    76    123    137    
     Running time :  6.062s

```


![2.StabSel](https://github.com/Leo-Simpson/c-lasso/blob/master/figures/advanced/StabSel.png)

![2.StabSel-beta](https://github.com/Leo-Simpson/c-lasso/blob/master/figures/advanced/StabSel-beta.png)

![2.CV-beta](https://github.com/Leo-Simpson/c-lasso/blob/master/figures/advanced/CVbeta.png)

![2.CV-graph](https://github.com/Leo-Simpson/c-lasso/blob/master/figures/advanced/CV.png)

![2.LAM-beta](https://github.com/Leo-Simpson/c-lasso/blob/master/figures/advanced/beta.png)

![2.Path](https://github.com/Leo-Simpson/c-lasso/blob/master/figures/advanced/Beta-path.png)


## Log-contrast regression for microbiome data

In the [the accompanying notebook](./examples/example-notebook.ipynb) we study several microbiome data sets. We showcase two examples below.

#### BMI prediction using the COMBO dataset 

We first consider the [COMBO data set](./examples/COMBO_data) and show how to predict Body Mass Index (BMI) from microbial genus abundances and two non-compositional covariates  using "filtered_data".

```python
from classo import csv_to_np, classo_problem, clr

# Load microbiome and covariate data X
X0  = csv_to_np('COMBO_data/complete_data/GeneraCounts.csv', begin = 0).astype(float)
X_C = csv_to_np('COMBO_data/CaloriData.csv', begin = 0).astype(float)
X_F = csv_to_np('COMBO_data/FatData.csv', begin = 0).astype(float)

# Load BMI measurements y
y   = csv_to_np('COMBO_data/BMI.csv', begin = 0).astype(float)[:, 0]
labels = csv_to_np('COMBO_data/complete_data/GeneraPhylo.csv').astype(str)[:, -1]


# Normalize/transform data
y   = y - np.mean(y) #BMI data (n = 96)
X_C = X_C - np.mean(X_C, axis = 0)  #Covariate data (Calorie)
X_F = X_F - np.mean(X_F, axis = 0)  #Covariate data (Fat)
X0 = clr(X0, 1 / 2).T

# Set up design matrix and zero-sum constraints for 45 genera
X     = np.concatenate((X0, X_C, X_F, np.ones((len(X0), 1))), axis = 1) # Joint microbiome and covariate data and offset
label = np.concatenate([labels, np.array(['Calorie', 'Fat', 'Bias'])])
C = np.ones((1, len(X[0])))
C[0, -1], C[0, -2], C[0, -3] = 0., 0., 0.


# Set up c-lassso problem
problem = classo_problem(X, y, C, label = label)


# Use stability selection with theoretical lambda [Combettes & Müller, 2020b]
problem.model_selection.StabSelparameters.method      = 'lam'
problem.model_selection.StabSelparameters.threshold_label = 0.5

# Use formulation R3
problem.formulation.concomitant = True

problem.solve()
print(problem)
print(problem.solution)

# Use formulation R4
problem.formulation.huber = True
problem.formulation.concomitant = True

problem.solve()
print(problem)
print(problem.solution)

```

![3.Stability profile R3](https://github.com/Leo-Simpson/c-lasso/blob/master/figures/exampleFilteredCOMBO/R3-StabSel.png)

![3.Beta solution R3](https://github.com/Leo-Simpson/c-lasso/blob/master/figures/exampleFilteredCOMBO/R3-StabSel-beta.png)

![3.Stability profile R4](https://github.com/Leo-Simpson/c-lasso/blob/master/figures/exampleFilteredCOMBO/R4-StabSel.png)

![3.Beta solution R4](https://github.com/Leo-Simpson/c-lasso/blob/master/figures/exampleFilteredCOMBO/R4-StabSel-beta.png)


<!---
<img src="https://i.imgur.com/8tFmM8T.png" alt="Central Park Soil Microbiome" height="250" align="right"/>
#### pH prediction using the Central Park soil dataset 
The next microbiome example considers the [Central Park Soil dataset](./examples/pH_data) from [Ramirez et al.](https://royalsocietypublishing.org/doi/full/10.1098/rspb.2014.1988). The sample locations are shown in the Figure on the right.)
-->

#### pH prediction using the 88 soils dataset

The next microbiome example considers the [88 soils dataset](./examples/pH_data) from [Lauber et al., 2009](https://pubmed.ncbi.nlm.nih.gov/19502440/).

The task is to predict pH concentration in the soil from microbial abundance data. A similar analysis is available
in [Tree-Aggregated Predictive Modeling of Microbiome Data](https://www.biorxiv.org/content/10.1101/2020.09.01.277632v1) 
with Central Park soil data from [Ramirez et al.](https://royalsocietypublishing.org/doi/full/10.1098/rspb.2014.1988).

Code to run this application is available in [the accompanying notebook](./examples/example-notebook.ipynb) under `pH data`. Below is a summary of a c-lasso problem instance (using the R3 formulation).
 
```
FORMULATION: R3
 
MODEL SELECTION COMPUTED:  
     Lambda fixed
     Path
     Stability selection
 
LAMBDA FIXED PARAMETERS: 
     numerical_method = Path-Alg
     rescaled lam : True
     threshold = 0.004
     lam : theoretical
     theoretical_lam = 0.2182
 
PATH PARAMETERS: 
     numerical_method : Path-Alg
     lamin = 0.001
     Nlam = 80
 
 
STABILITY SELECTION PARAMETERS: 
     numerical_method : Path-Alg
     method : lam
     B = 50
     q = 10
     percent_nS = 0.5
     threshold = 0.7
     lam = theoretical
     theoretical_lam = 0.3085
```

The c-lasso estimation results are summarized below:

```
LAMBDA FIXED : 
   Sigma  =  0.198
   Selected variables :  14    18    19    39    43    57    62    85    93    94    104    107    
   Running time :  0.008s

 PATH COMPUTATION : 
   Running time :  0.12s

 STABILITY SELECTION : 
   Selected variables :  2    12    15    
   Running time :  0.287s
```

![Ex4.1](https://github.com/Leo-Simpson/c-lasso/blob/master/figures/examplePH/R3-Beta-path.png)

![Ex4.2](https://github.com/Leo-Simpson/c-lasso/blob/master/figures/examplePH/R3-Sigma-path.png)

![Ex4.3](https://github.com/Leo-Simpson/c-lasso/blob/master/figures/examplePH/R3-StabSel.png)

![Ex4.4](https://github.com/Leo-Simpson/c-lasso/blob/master/figures/examplePH/R3-StabSel-beta.png)

![Ex4.5](https://github.com/Leo-Simpson/c-lasso/blob/master/figures/examplePH/R3-beta.png)


## Optimization schemes

The available problem formulations [R1-C2] require different algorithmic strategies for 
efficiently solving the underlying optimization problem. We have implemented four 
algorithms (with provable convergence guarantees) that vary in generality and are not 
necessarily applicable to all problems. For each problem type, c-lasso has a default algorithm 
setting that proved to be the fastest in our numerical experiments.

### Path algorithms (Path-Alg) 
This is the default algorithm for non-concomitant problems [R1,R3,C1,C2]. 
The algorithm uses the fact that the solution path along &lambda; is piecewise-
affine (as shown, e.g., in [1]). When Least-Squares is used as objective function,
we derive a novel efficient procedure that allows us to also derive the 
solution for the concomitant problem [R2] along the path with little extra computational overhead.

### Projected primal-dual splitting method (P-PDS):
This algorithm is derived from [2] and belongs to the class of 
proximal splitting algorithms. It extends the classical Forward-Backward (FB) 
(aka proximal gradient descent) algorithm to handle an additional linear equality constraint
via projection. In the absence of a linear constraint, the method reduces to FB.
This method can solve problem [R1]. For the Huber problem [R3], 
P-PDS can solve the mean-shift formulation of the problem (see [6]).

### Projection-free primal-dual splitting method (PF-PDS):
This algorithm is a special case of an algorithm proposed in [3] (Eq.4.5) and also belongs to the class of 
proximal splitting algorithms. The algorithm does not require projection operators 
which may be beneficial when C has a more complex structure. In the absence of a linear constraint, 
the method reduces to the Forward-Backward-Forward scheme. This method can solve problem [R1]. 
For the Huber problem [R3], PF-PDS can solve the mean-shift formulation of the problem (see [6]).

### Douglas-Rachford-type splitting method (DR)
This algorithm is the most general algorithm and can solve all regression problems 
[R1-R4]. It is based on Doulgas Rachford splitting in a higher-dimensional product space.
It makes use of the proximity operators of the perspective of the LS objective (see [4,5])
The Huber problem with concomitant scale [R4] is reformulated as scaled Lasso problem 
with the mean shift (see [6]) and thus solved in (n + d) dimensions. 



## References 

* [1] B. R. Gaines, J. Kim, and H. Zhou, [Algorithms for Fitting the Constrained Lasso](https://www.tandfonline.com/doi/abs/10.1080/10618600.2018.1473777?journalCode=ucgs20), J. Comput. Graph. Stat., vol. 27, no. 4, pp. 861–871, 2018.

* [2] L. Briceno-Arias and S.L. Rivera, [A Projected Primal–Dual Method for Solving Constrained Monotone Inclusions](https://link.springer.com/article/10.1007/s10957-018-1430-2?shared-article-renderer), J. Optim. Theory Appl., vol. 180, Issue 3, March 2019.

* [3] P. L. Combettes and J.C. Pesquet, [Primal-Dual Splitting Algorithm for Solving Inclusions with Mixtures of Composite, Lipschitzian, and Parallel-Sum Type Monotone Operators](https://arxiv.org/pdf/1107.0081.pdf), Set-Valued and Variational Analysis, vol. 20, pp. 307-330, 2012.

* [4] P. L. Combettes and C. L. Müller, [Perspective M-estimation via proximal decomposition](https://arxiv.org/abs/1805.06098), Electronic Journal of Statistics, 2020, [Journal version](https://projecteuclid.org/euclid.ejs/1578452535) 

* [5] P. L. Combettes and C. L. Müller, [Regression models for compositional data: General log-contrast formulations, proximal optimization, and microbiome data applications](https://arxiv.org/abs/1903.01050), Statistics in Bioscience, 2020.

* [6] A. Mishra and C. L. Müller, [Robust regression with compositional covariates](https://arxiv.org/abs/1909.04990), arXiv, 2019.

* [7] S. Rosset and J. Zhu, [Piecewise linear regularized solution paths](https://projecteuclid.org/euclid.aos/1185303996), Ann. Stat., vol. 35, no. 3, pp. 1012–1030, 2007.

* [8] J. Bien, X. Yan, L. Simpson, and C. L. Müller,   [Tree-Aggregated Predictive Modeling of Microbiome Data](https://www.biorxiv.org/content/10.1101/2020.09.01.277632v1), biorxiv, 2020.


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
reported by contacting the project team at leo.bill.simpson@gmail.com. All
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
---
title: 'c-lasso - a Python package for constrained sparse and robust regression and classification'
tags:
  - Python
  - regression
  - classification
  - constrained regression
  - Lasso
  - Huber function
  - Square Hinge SVM
  - convex optimization
  - perspective function
authors:
  - name: Léo Simpson
    affiliation: 1
  - name: Patrick L. Combettes
    affiliation: 2
  - name: Christian L. Müller
    orcid: 0000-0002-3821-7083
    affiliation: "3,4,5"

affiliations:
  - name: Technische Universität München
    index: 1
  - name: Department of Mathematics, North Carolina State University, Raleigh
    index: 2
  - name: Center for Computational Mathematics, Flatiron Institute, New York
    index: 3
  - name: Institute of Computational Biology, Helmholtz Zentrum München
    index: 4
  - name: Department of Statistics, Ludwig-Maximilians-Universität München
    index: 5
date: 09 October 2020
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:

---

# Summary

We introduce `c-lasso`, a Python package that enables sparse and robust linear regression and classification with linear equality constraints. 
The underlying statistical forward model is assumed to be of the following form:

$$
y = X \beta + \sigma \epsilon \qquad \textrm{subject to} \qquad C\beta=0
$$

Here, $X \in \mathbb{R}^{n\times d}$ is a given design matrix and the vector $y \in \mathbb{R}^{n}$ is a continuous or binary response vector. The matrix $C$ is a general
constraint matrix. The vector $\beta \in \mathbb{R}^{d}$ contains the unknown coefficients and $\sigma$ an unknown scale. Prominent use cases are (sparse) log-contrast
regression with compositional data $X$,  requiring the constraint  $1_d^T \beta = 0$ [@Aitchison:1984] and the Generalized Lasso which is a *special case* of the described problem (see, e.g, [@James:2020], Example 3). The `c-lasso` package provides estimators for inferring unknown coefficients and scale (i.e., perspective M-estimators [@Combettes:2020a]) of the form

$$
    \min_{\beta \in \mathbb{R}^d, \sigma \in \mathbb{R}_{0}} f\left(X\beta - y,{\sigma} \right) + \lambda \left\lVert \beta\right\rVert_1 \qquad \textrm{subject to} \qquad  C\beta = 0
$$

for several convex loss functions $f(\cdot,\cdot)$. This includes the constrained Lasso, the constrained scaled Lasso, sparse Huber M-estimators with linear equality constraints, and constrained (Huberized) Square Hinge Support Vector Machines (SVMs) for classification.

# Statement of need 

Currently, there is no Python package available that can solve these ubiquitous statistical estimation problems in a fast and efficient manner. 
`c-lasso` provides algorithmic strategies, including path and proximal splitting algorithms, to solve the underlying convex optimization problems with provable convergence guarantees. The `c-lasso` package is intended to fill the gap between popular Python tools such as [`scikit-learn`](https://scikit-learn.org/stable/) which <em>cannot</em> solve these constrained problems and general-purpose optimization solvers such as [`cvxpy`](https://www.cvxpy.org) that do not scale well for these problems and/or are inaccurate. `c-lasso` can solve the estimation problems at a single regularization level, across an entire regularization path, and includes three model selection strategies for determining the regularization parameter: a theoretically-derived fixed regularization, k-fold cross-validation, and stability selection. We show several use cases of the package, including an application of sparse log-contrast regression tasks for compositional microbiome data, and highlight the seamless integration into `R` via [`reticulate`](https://rstudio.github.io/reticulate/).

# Functionalities

## Installation and problem instantiation {#gettingstarted}

`c-lasso` is available on pip and can be installed in the shell using

```shell
pip install c-lasso
```

`c-lasso` is a stand-alone package and not yet compatible with the `scikit-learn` API.
The central object in the `c-lasso` package is the instantiation of a `c-lasso` problem. 

```python
# Import the main class of the package
from classo import classo_problem

# Define a c-lasso problem instance with default setting, 
# given data X, y, and constraints C.
problem  = classo_problem(X, y, C)
```

We next describe what type of problem instances are available and how to solve them.

## Statistical problem formulations {#formulations}

Depending on the type of and the prior assumptions on the data, the noise $\epsilon$, and the model parameters, `c-lasso` allows 
for different estimation problem formulations. More specifically, the package can solve the following 
four regression-type and two classification-type formulations:


### *R1* Standard constrained Lasso regression: {#R1}           

$$
    \min_{\beta \in \mathbb{R}^d} \left\lVert X\beta - y \right\rVert^2 + \lambda \left\lVert \beta\right\rVert_1 \qquad \textrm{subject to} \qquad  C\beta = 0
$$

This is the standard Lasso problem with linear equality constraints on the $\beta$ vector. 
The objective function combines Least-Squares (LS) for model fitting with the $L_1$-norm penalty for sparsity.   

```python
# Formulation R1
problem.formulation.huber = False
problem.formulation.concomitant = False
problem.formulation.classification = False
```

### *R2* Constrained sparse Huber regression: {#R2}                  

$$
    \min_{\beta \in \mathbb{R}^d} h_{\rho} (X\beta - y) + \lambda \left\lVert \beta\right\rVert_1 \qquad \textrm{subject to} \qquad  C\beta = 0
$$

This regression problem uses the [Huber loss](https://en.wikipedia.org/wiki/Huber_loss) $h_{\rho}$ as objective function 
for robust model fitting with an $L_1$ penalty and linear equality constraints on the $\beta$ vector. The default parameter $\rho$ is set to $1.345$ [@Huber:1981].

```python
# Formulation R2
problem.formulation.huber = True
problem.formulation.concomitant = False
problem.formulation.classification = False
```

### *R3* Constrained scaled Lasso regression: {#R3}

$$
    \min_{\beta \in \mathbb{R}^d, \sigma \in \mathbb{R}_{0}} \frac{\left\lVert X\beta - y \right\rVert^2}{\sigma} + \frac{n}{2} \sigma + \lambda \left\lVert \beta\right\rVert_1 \qquad \textrm{subject to} \qquad  C\beta = 0
$$

This formulation is the default problem formulation in `c-lasso`. It is similar to [*R1*](#R1) but allows for joint estimation of the (constrained) $\beta$ vector and the standard deviation $\sigma$ in a concomitant fashion [@Combettes:2020a;@Combettes:2020b].

```python
# Formulation R3
problem.formulation.huber = False
problem.formulation.concomitant = True
problem.formulation.classification = False
```

### *R4* Constrained sparse Huber regression with concomitant scale estimation: {#R4}       

$$
    \min_{\beta \in \mathbb{R}^d, \sigma \in  \mathbb{R}_{0}} \left( h_{\rho} \left( \frac{X\beta - y}{\sigma} \right) + n \right) \sigma + \lambda \left\lVert \beta\right\rVert_1 \qquad \textrm{subject to} \qquad  C\beta = 0
$$

This formulation combines [*R2*](#R2) and [*R3*](#R3) allowing robust joint estimation of the (constrained) $\beta$ vector and 
the scale $\sigma$ in a concomitant fashion [@Combettes:2020a;@Combettes:2020b].

```python
# Formulation R4
problem.formulation.huber = True
problem.formulation.concomitant = True
problem.formulation.classification = False
```

### *C1* Constrained sparse classification with Square Hinge loss: {#C1}

$$
    \min_{\beta \in \mathbb{R}^d} \sum_{i=1}^n l(y_i x_i^\top\beta) + \lambda \left\lVert \beta\right\rVert_1 \qquad \textrm{subject to} \qquad  C\beta = 0
$$

where $x_i$ denotes the $i^{th}$ row of $X$, $y_i \in \{-1,1\}$, and $l(\cdot)$ is defined for $r \in \mathbb{R}$ as:

$$
l(r) = \begin{cases} (1-r)^2 & if \quad r \leq 1 \\ 0 &if \quad r \geq 1 \end{cases}
$$

This formulation is similar to [*R1*](#R1) but adapted for classification tasks using the Square Hinge loss with (constrained) sparse $\beta$ vector estimation [@Lee:2013].

```python
# Formulation C1
problem.formulation.huber = False
problem.formulation.concomitant = False
problem.formulation.classification = True
```

### *C2* Constrained sparse classification with Huberized Square Hinge loss: {#C2}       

$$
    \min_{\beta \in \mathbb{R}^d}  \sum_{i=1}^n  l_{\rho}(y_i x_i^\top\beta) + \lambda \left\lVert \beta\right\rVert_1 \qquad \textrm{subject to} \qquad  C\beta = 0 \,.
$$

This formulation is similar to [*C1*](#C1) but uses the Huberized Square Hinge loss $l_{\rho}$ for robust classification with (constrained) sparse $\beta$ vector estimation [@Rosset:2007]:

$$
l_{\rho}(r) = \begin{cases} (1-r)^2 &if \quad \rho \leq r \leq 1 \\ (1-\rho)(1+\rho-2r) & if \quad r \leq \rho \\ 0 &if \quad r \geq 1 \end{cases}
$$

This formulation can be selected in `c-lasso` as follows:

```python
# Formulation C2
problem.formulation.huber = True
problem.formulation.concomitant = False
problem.formulation.classification = True
```

## Optimization schemes {#method}

The problem formulations *R1*-*C2* require different algorithmic strategies for efficiently solving the underlying optimization problems. The `c-lasso` package implements four published algorithms with provable convergence guarantees. The package also includes novel algorithmic extensions to solve Huber-type problems using the mean-shift formulation [@Mishra:2019]. The following algorithmic schemes are implemented: 

- Path algorithms (*Path-Alg*): 
This algorithm follows the proposal in [@Gaines:2018;@Jeon:2020] and uses the fact that the solution path along $\lambda$ is piecewise-affine [@Rosset:2007]. We also provide a novel efficient procedure that allows to derive the solution for the concomitant problem *R3* along the path with little computational overhead.

- Douglas-Rachford-type splitting method (*DR*): 
This algorithm can solve all regression problems *R1-R4*. It is based on Doulgas-Rachford splitting in a higher-dimensional product space and 
makes use of the proximity operators of the perspective of the LS objective [@Combettes:2020a;@Combettes:2020b]. The Huber problem with concomitant scale *R4* is reformulated as scaled Lasso problem with mean shift vector [@Mishra:2019] and thus solved in (n + d) dimensions.

- Projected primal-dual splitting method (*P-PDS*): 
This algorithm is derived from [@Briceno:2020] and belongs to the class of proximal splitting algorithms, extending the classical Forward-Backward (FB) 
(aka proximal gradient descent) algorithm to handle an additional linear equality constraint via projection. In the absence of a linear constraint, the method reduces to FB.

- Projection-free primal-dual splitting method (*PF-PDS*):
This algorithm is a special case of an algorithm proposed in [@Combettes:2012] (Eq. 4.5) and also belongs to the class of 
proximal splitting algorithms. The algorithm does not require projection operators which may be beneficial when C has a more complex structure. 
In the absence of a linear constraint, the method reduces to the Forward-Backward-Forward scheme.

The following table summarizes the available algorithms and their recommended use for each problem: 

|             |*Path-Alg*| *DR* | *P-PDS* | *PF-PDS* |
|-|:-:|:-:|:-:|:-:|
| [*R1*](#R1) | use for large $\lambda$ and path computation | use for small $\lambda$ | possible | use for complex constraints |
| [*R2*](#R2) | use for large $\lambda$ and path computation | use for small $\lambda$ | possible | use for complex constraints |
| [*R3*](#R3) | use for large $\lambda$ and path computation | use for small $\lambda$ |    -      |        -           |
| [*R4*](#R4) |           -                                                         | only option                |    -      |  - |-
| [*C1*](#C1) | only option                                                      |                 -                |   -       |  - |
| [*C2*](#C2) | only option                                                     |                -                 |     -     |  - |



The following Python snippet shows how to select a specific algorithm: 
```python
problem.numerical_method = "Path-Alg" 
# Alternative options: "DR", "P-PDS", and "PF-PDS" 
```

## Computation modes and model selection {#model}

The `c-lasso` package provides several computation modes and model selection schemes for tuning the regularization parameter.

- *Fixed Lambda*: This setting lets the user choose a fixed parameter $\lambda$ or a proportion $l \in [0,1]$ such that $\lambda = l\times \lambda_{\max}$. 
The default value is a scale-dependent tuning parameter that has been derived in [@Shi:2016] and applied in [@Combettes:2020b].

- *Path Computation*: This setting allows the computation of a solution path for $\lambda$ parameters in an interval $[\lambda_{\min}, \lambda_{\max}]$. The solution path is computed via the *Path-Alg* scheme or via warm-starts for other optimization schemes. 

[comment]: <> (This can be done much faster than by computing separately the solution for each $\lambda$ of the grid, by using the Path-alg algorithm. One can also use warm starts: starting with $\beta_0 = 0$ for $\lambda_0 = \lambda_{\max}$, and then iteratvely compute $\beta_{k+1}$ using one of the optimization schemes with $\lambda = \lambda_{k+1} := \lambda_{k} - \epsilon$ and with a warm start set to $\beta_{k}$. )

- *Cross Validation*: This setting allows the selection of the regularization parameter $\lambda$ via k-fold cross validation for $\lambda \in [\lambda_{\min}, \lambda_{\max}]$. Both the Minimum Mean Squared Error (or Deviance) (MSE)  and the "One-Standard-Error rule" (1SE) are available [@Hastie:2009].

- *Stability Selection*: This setting allows the selection of the $\lambda$ via stability selection [@Meinshausen:2010;@Lin:2014;@Combettes:2020b]. Three modes are available: selection at a fixed $\lambda$ [@Combettes:2020b], selection of the q *first* variables entering the path (default setting), and of the q *largest coefficients* (in absolute value) across the path [@Meinshausen:2010].

The Python syntax to use a specific computation mode and model selection is exemplified below:

```python
# Example how to perform ath computation and cross-validation:
problem.model_selection.LAMfixed = False
problem.model_selection.PATH = True
problem.model_selection.CV = True
problem.model_selection.StabSel = False

# Example how to add stability selection to the problem instance
problem.model_selection.StabSel = True
```
Each model selection procedure has additional meta-parameters that are described in the [Documentation](https://c-lasso.readthedocs.io/en/latest/).

# Numerical benchmarks

To evaluate optimization accuracy and running time of the different algorithms available in `c-lasso`, we provide [micro-benchmark](https://github.com/Leo-Simpson/c-lasso/tree/master/benchmark) experiments which also include [cvxpy](https://www.cvxpy.org), an open source convex optimization software, for baseline comparison. All experiments have been computed using Python 3.9.1 on a `MacBook Air` with a `1.8 GHz Intel Core i5` processor and `8 Gb 1600 MHz DDR3` memory, operating on macOS High Sierra. 

Figure 1 summarizes the results for the *Path-Alg*, *DR*, and *P-PDS* algorithms solving the regression formulation [R1](#R1) for different samples sizes $n$ and problem dimensions $p$ on synthetic data (using `c-lasso`'s data generator). We observe that `c-lasso`'s algorithms are faster and more accurate than the `cvx` baseline. For instance, for $d=500$ features and $n=500$ samples, the *Path-Alg* algorithm is about $70$ times faster than `cvx`.

![Average running times (left panel) of *Path-Alg* (blue), *P-PDS* (yellow), *DR* (green), and cvx (red) at fixed $\lambda = 0.1$ and corresponding average objective function value differences (with respect to the function value obtained by the *Path-Alg* solution as baseline) (right panel). Mean (and standard deviation) running time is calculated over 20 data replications for each sample size/dimension scenario $(n,d)$. On a single data set, the reported running time of an algorithm is the average time of five algorithm runs (to guard against system background process fluctuations).](figures/figure_benchmark.png)

The complete reproducible micro-benchmark is avaialable [here](https://github.com/Leo-Simpson/c-lasso/tree/master/benchmark).

# Computational examples  

## Toy example using synthetic data

We illustrate the workflow of the `c-lasso` package on synthetic data using the built-in routine ```random_data``` which enables the generation of test 
problem instances with normally distributed data $X$, sparse coefficient vectors $\beta$, and constraints $C \in \mathbb{R}^{k\times d}$.

Here, we use a problem instance with $n=100$, $d=100$, a $\beta$ with five non-zero components, $\sigma=0.5$, and a zero-sum contraint.


```python






from classo import classo_problem, random_data

n, d, d_nonzero, k, sigma = 100, 100, 5, 1, 0.5
(X, C, y), sol = random_data(
  n, d, d_nonzero, k, sigma,
  zerosum = True, seed = 123
)
print("Relevant variables  : {}".format(numpy.nonzero(sol)[0]))

problem = classo_problem(X, y, C)

problem.formulation.huber = True
problem.formulation.concomitant = False
problem.formulation.rho = 1.5

problem.model_selection.LAMfixed = True
problem.model_selection.PATH = True
problem.model_selection.LAMfixedparameters.rescaled_lam = True
problem.model_selection.LAMfixedparameters.lam = 0.1

problem.solve()

print(problem.solution)
```

We use [formulation](#formulations) [*R2*](#R2) with $\rho=1.5$, [computation mode and model selections](#model) *Fixed Lambda* with $\lambda = 0.1\lambda_{\max}$, *Path computation*, and *Stability Selection* (as per default). 

The corresponding output reads: 

```shell
Relevant variables  : [43 47 74 79 84]

 LAMBDA FIXED : 
   Selected variables :  43    47    74    79    84    
   Running time :  0.294s

 PATH COMPUTATION : 
   Running time :  0.566s

 STABILITY SELECTION : 
   Selected variables :  43    47    74    79    84    
   Running time :  5.3s
```

`c-lasso` allows standard visualization of the computed solutions, e.g., coefficient plots at fixed $\lambda$, the solution path, the stability selection 
profile at the selected $\lambda$, and the stability selection profile across the entire path. 

![Visualizations after calling problem.solution](figures/synthetic.png)

For this tuned example, the solutions at the fixed lambda and with stability selection recover the oracle solution. The solution vectors
are stored in ```problem.solution``` and can be directly acccessed for each mode/model selection. 

```python
# Access to the estimated coefficient vector at a fixed lambda 
problem.solution.LAMfixed.beta
```
Note that the run time for this $d=100$-dimensional example for a single path computation is about 0.5 seconds on a standard laptop.

## Log-contrast regression on gut microbiome data

We next illustrate the application of `c-lasso` on the [`COMBO` microbiome dataset](https://github.com/Leo-Simpson/c-lasso/tree/master/examples/COMBO_data) [@Lin:2014;@Shi:2016;@Combettes:2020b]. Here, the task is to predict the Body Mass Index (BMI) of $n=96$ participants from $d=45$ relative abundances of bacterial genera, and absolute calorie and fat intake measurements. The code snippet for this example is available in the [`README.md`](https://github.com/Leo-Simpson/c-lasso/README.md) and the [example notebook](https://github.com/Leo-Simpson/c-lasso/blob/master/examples/example-notebook.ipynb).

![Stability selection profiles of problems R3/R4 on the COMBO data](figures/StabSelFilteredCOMBO.png)

Stability selection profiles using [formulation](#formulations) [*R3*](#R3) (left) and [*R4*](#R4)(right) on the COMBO dataset, reproducing Figure 5a in [@Combettes:2020b].

## Calling `c-lasso` in R 

The `c-lasso` package also integrates with `R` via the `R` package [`reticulate`](https://rstudio.github.io/reticulate/). We refer to 
`reticulate`'s manual for technical details about connecting `python` environments and `R`. A successful use case of `c-lasso` is available in the `R` package [`trac`](https://github.com/jacobbien/trac) [@Bien:2020], enabling tree-structured aggregation of predictors when features are rare.

The code snippet below shows how `c-lasso` is called in `R` to perform regression at a fixed $\lambda$ $\lambda = 0.1\lambda_{\max}$. 
In `R`, X and C need to be of ```matrix``` type, and y of ```array``` type.

```r
problem <- classo$classo_problem(X = X, C = C, y = y) 
problem$model_selection$LAMfixed <- TRUE
problem$model_selection$StabSel <- FALSE
problem$model_selection$LAMfixedparameters$rescaled_lam <- TRUE
problem$model_selection$LAMfixedparameters$lam <- 0.1
problem$solve()

# Extract coefficent vector with tidy-verse
beta <- as.matrix(map_dfc(problem$solution$LAMfixed$beta, as.numeric))
```

# Acknowledgements

The work of LS was conducted at and financially supported by the Center for Computational Mathematics (CCM), Flatiron Institute, New York, and the Institute of Computational Biology, Helmholtz Zentrum München. We thank Dr. Leslie Greengard (CCM and Courant Institute, NYU) for facilitating the initial contact between LS and CLM. The work of PLC was supported by the National Science Foundation under grant DMS-1818946. 

# References



---
title: 'c-lasso - a package for constrained sparse and robust regression and classification in Python'
tags:
  - Python
  - Linear regression
  - Optimization
authors:
  - name: Léo Simpson
    affiliation: 1
  - name: Patrick L. Combettes
    affiliation: 2
  - name: Christian L. Müller
    affiliation: 3,4,5
affiliations:
 - name: Technische Universität München 
   index: 1
 - name: Department of Mathematics, North Carolina State University
   index: 2
 - name: Center for Computational Mathematics, Flatiron Institute
   index: 3
 - name: Institute of Computational Biology, Helmholtz Zentrum München
   index: 4
 - name: Department of Statistics, Ludwig-Maximilians-Universität München
   index: 5
   
date: 01 October 2020
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:

---

# Summary

We introduce c-lasso, a Python package that enables sparse and robust linear
regression and classification with linear equality constraints. 


The forward model is assumed to be:

$$
y = X \beta + \sigma \epsilon \qquad \textrm{s.t.} \qquad C\beta=0
$$

Here, $X \in R^{n\times d}$ is a given design matrix and the vector $y \in R^{n}$ is a continuous or binary response vector. The matrix $C$ is a general constraint matrix. The vector $\beta \in R^{d}$ contains the unknown coefficients and $\sigma$ an unknown scale.


# Statement of need 

The package handles several estimators for inferring coefficients and scale, including the constrained Lasso, the constrained scaled Lasso, and sparse Huber M-estimators with linear equality constraints, none of which can be currently solved in Python in an efficient manner. Several algorithmic strategies, including path and proximal splitting algorithms, are implemented to solve the underlying convex optimization problems at fixed regularization and across an entire regularization path. We include three model selection strategies for determining model parameter regularization levels: a theoretically derived fixed regularization, k-fold cross-validation, and stability selection. The c-lasso package is intended to fill the gap between popular python tools such as `scikit-learn` which <em>cannot</em> solve sparse constrained problems and general-purpose optimization solvers such as `cvx` that do not scale well for the considered problems and/or are inaccurate. We show several use cases of the package, including an application of sparse log-contrast regression tasks for compositional microbiome data. We also highlight the seamless integration of the solver into `R` via the `reticulate` package. 


# Current functionalities

## Formulations 

Depending on the prior on the solution $\beta, \sigma$ and on the noise $\epsilon$, the previous forward model can lead to different types of estimation problems. 
The package can solve four regression-type and two classification-type formulations:


### *R1* Standard constrained Lasso regression:             

$$
    \arg \min_{\beta \in \mathbb{R}^d} \left\lVert X\beta - y \right\rVert^2 + \lambda \left\lVert \beta\right\rVert_1 \qquad s.t. \qquad  C\beta = 0
$$

This is the standard Lasso problem with linear equality constraints on the $\beta$ vector. 
The objective function combines Least-Squares for model fitting with $L_1$ penalty for sparsity.   
 

### *R2* Contrained sparse Huber regression:                   

$$
    \arg \min_{\beta \in \mathbb{R}^d} h_{\rho} (X\beta - y) + \lambda \left\lVert \beta\right\rVert_1 \qquad s.t. \qquad  C\beta = 0
$$

This regression problem uses the [Huber loss](https://en.wikipedia.org/wiki/Huber_loss) as objective function 
for robust model fitting with an $L_1$ penalty and linear equality constraints on the $\beta$ vector. The parameter $\rho$ is set to $1.345$ by default [@Huber:1981].

### *R3* Contrained scaled Lasso regression: 

$$
    \arg \min_{\beta \in \mathbb{R}^d} \frac{\left\lVert X\beta - y \right\rVert^2}{\sigma} + \frac{n}{2} \sigma + \lambda \left\lVert \beta\right\rVert_1 \qquad s.t. \qquad  C\beta = 0
$$

This formulation is similar to *R1* but allows for joint estimation of the (constrained) $\beta$ vector and 
the standard deviation $\sigma$ in a concomitant fashion [@Combettes:2020.1; @Combettes:2020.2].
This is the default problem formulation in c-lasso.

### *R4* Contrained sparse Huber regression with concomitant scale estimation:        

$$
    \arg \min_{\beta \in \mathbb{R}^d} \left( h_{\rho} \left( \frac{X\beta - y}{\sigma} \right) + n \right) \sigma + \lambda \left\lVert \beta\right\rVert_1 \qquad s.t. \qquad  C\beta = 0
$$

This formulation combines *R2* and *R3* allowing robust joint estimation of the (constrained) $\beta$ vector and 
the scale $\sigma$ in a concomitant fashion [@Combettes:2020.1; @Combettes:2020.2].

### *C1* Contrained sparse classification with Square Hinge loss: 

$$
    \arg \min_{\beta \in \mathbb{R}^d} L(y^T X\beta - y) + \lambda \left\lVert \beta\right\rVert_1 \qquad s.t. \qquad  C\beta = 0
$$

where $L \left((r_1,...,r_n)^T \right) := \sum_{i=1}^n l(r_i)$ and $l$ is defined as :

$$
l(r) = \begin{cases} (1-r)^2 & if \quad r \leq 1 \\ 0 &if \quad r \geq 1 \end{cases}
$$

This formulation is similar to *R1* but adapted for classification tasks using the Square Hinge loss with (constrained) sparse $\beta$ vector estimation [@Lee:2013].

### *C2* Contrained sparse classification with Huberized Square Hinge loss:        

$$
    \arg \min_{\beta \in \mathbb{R}^d} L_{\rho}(y^T X\beta - y) + \lambda \left\lVert \beta\right\rVert_1 \qquad s.t. \qquad  C\beta = 0
$$

where $L_{\rho} \left((r_1,...,r_n)^T \right) := \sum_{i=1}^n l_{\rho}(r_i)$ and $l_{\rho}$ is defined as :

$$
l_{\rho}(r) = \begin{cases} (1-r)^2 &if \quad \rho \leq r \leq 1 \\ (1-\rho)(1+\rho-2r) & if \quad r \leq \rho \\ 0 &if \quad r \geq 1 \end{cases}
$$

This formulation is similar to *C1* but uses the Huberized Square Hinge loss for robust classification with (constrained) sparse $\beta$ vector estimation [@Rosset:2007].




## Optimization schemes

The available problem formulations *R1-C2* require different algorithmic strategies for 
efficiently solving the underlying optimization problem. We have implemented four 
algorithms (with provable convergence guarantees) that vary in generality and are not 
necessarily applicable to all problems. For each problem type, c-lasso has a default algorithm 
setting that proved to be the fastest in our numerical experiments.

- **Path algorithms (*Path-Alg*)** : 
This is the default algorithm for non-concomitant problems *R1,R2,C1,C2*. 
The algorithm uses the fact that the solution path along &lambda; is piecewise-affine (as shown, e.g., in [@Gaines:2018]). When Least-Squares is used as objective function, we provide a novel efficient procedure that also allows to derive the solution for the concomitant problem *R3* along the path with little computational overhead.

- **Projected primal-dual splitting method (*P-PDS*)** : 
This algorithm is derived from [@Briceno:2020] and belongs to the class of 
proximal splitting algorithms. It extends the classical Forward-Backward (FB) 
(aka proximal gradient descent) algorithm to handle an additional linear equality constraint
via projection. In the absence of a linear constraint, the method reduces to FB.
This method can solve problem *R1*. For the Huber problem *R2*, 
P-PDS can solve the mean-shift formulation of the problem [@Mishra:2019].

- **Projection-free primal-dual splitting method (*PF-PDS*)** :
This algorithm is a special case of an algorithm proposed in [@Combettes:2011] (Eq.4.5) and also belongs to the class of 
proximal splitting algorithms. The algorithm does not require projection operators 
which may be beneficial when C has a more complex structure. In the absence of a linear constraint, 
the method reduces to the Forward-Backward-Forward scheme. This method can solve problem *R1*. 
For the Huber problem *R2*, PF-PDS can solve the mean-shift formulation of the problem [@Mishra:2019].

- **Douglas-Rachford-type splitting method (*DR*)** : 
This algorithm is the most general algorithm and can solve all regression problems 
*R1-R4*. It is based on Doulgas-Rachford splitting in a higher-dimensional product space and 
makes use of the proximity operators of the perspective of the LS objective [@Combettes:2020.1; @Combettes:2020.2].
The Huber problem with concomitant scale *R4* is reformulated as scaled Lasso problem 
with the mean shift [@Mishra:2019] and thus solved in (n + d) dimensions.

## Model selections and computation modes

Different models are implemented together with the optimization schemes, to overcome the difficulty of choosing the penalization free parameter $\lambda$. 

- *Fixed Lambda*: This setting lets the user choose a parameter $\lambda$, or a proportion $l \in [0,1]$ such that $\lambda = l\times \lambda_{\max}$. 
The default value is a scale-dependent tuning parameter that has been derived in [@Shi:2016] and applied in [Combettes:2020.2].

- *Path Computation*: This setting allows the computation of a solution path for $\lambda$ parameters in an interval $[\lambda_{\min}, \lambda_{\max}]$. The solution path is computed via the *Path-Alg* scheme or via warm-starts for other optimization schemes. 

[comment]: <> (This can be done much faster than by computing separately the solution for each $\lambda$ of the grid, by using the Path-alg algorithm. One can also use warm starts : starting with $\beta_0 = 0$ for $\lambda_0 = \lambda_{\max}$, and then iteratvely compute $\beta_{k+1}$ using one of the optimization schemes with $\lambda = \lambda_{k+1} := \lambda_{k} - \epsilon$ and with a warm start set to $\beta_{k}$. )

- *Cross Validation*: This setting allows the selection of the regularization parameter $\lambda$ via k-fold cross validation for $\lambda \in [\lambda_{\min}, \lambda_{\max}]$. Both the Minimum Mean Squared Error (or Deviance) (MSE)  and the "One-Standard-Error rule" (1SE) are available [@Hastie:2009].

- *Stability Selection*: This setting allows the selection of the $\lambda$ via stability selection [@Lin:2014; @Meinshausen:2010; Combettes:2020.2].

# Basic workflow

Here is a basic example how to install and run c-lasso on synthetic data.

#### Installation 
c-lasso is available on pip. You can install the package
in the shell using

```shell
pip install c_lasso
```
To use the c-lasso package in Python, type 

```python
from classo import *
```

The c-lasso package depends on the following Python packages,
included in the package: 

`numpy` ; 
`matplotlib` ; 
`scipy` ; 
`pandas` ; 
`h5py` . 


#### Generate random data
The c-lasso package includes
the routine ```random_data``` that allows to generate problem instances with normally distributed data $X$, sparse $\beta$, and constraints $C$.

```python
n,d,d_nonzero,k,sigma =100,100,5,1,0.5
(X,C,y),sol = random_data(n,d,d_nonzero,k,sigma,zerosum=True, seed = 123 )
```
This code snippet generates randomly the vectors $\beta \in R^d$ , $X \in R^{n\times d}$ , $C \in R^{k\times d}$ (here it is the all-one vector instead because of the input ```zerosum```), and $y \in R^n$ normally distributed with respect to the model $C\beta=0$, $y-X\beta \sim N(0,I_n\sigma^2)$ and $\beta$ has only d_nonzero non-null componant.

#### Use of c-lasso on the data
Here is an example of problem instance one can create with those set of data. 

```python
# Define a c-lasso problem instance with default setting
problem  = classo_problem(X,y,C)


# Example how to change the default formulation R3 to formulation R2 with a new rho parameter
problem.formulation.huber  = True
problem.formulation.concomitant = False
problem.formulation.rho = 1.5


# Example how to add the computation of beta at a fixed lambda (as a proportion of lambdamax) 
problem.model_selection.LAMfixed = True
# and set it to to 0.1*lambdamax
problem.model_selection.LAMfixedparameters.rescaled_lam = True
problem.model_selection.LAMfixedparameters.lam = 0.1

# Example how to a computation of the lambda-path
problem.model_selection.PATH = True


# Solve the specified problem instance
problem.solve()
```

Here, we modified the [formulation](##formulations) of the problem in order to use [*R2*](###*R2*-contrained-sparse-Huber-regression), with $\rho=1.5$. 

Then, we chose the [model selections](##model-selections) we want to compute: *Fixed Lambda* with $\lambda = 0.1\lambda_{\max}$ and *Path computation*. 
*Stability Selection* is also computed by default. 

Finally, these problem specifications are solved using the method `solve`. 

#### Visualizing the problem setup and solutions 
c-lasso enables the visualization of the problem specifications and (after solve) the problem solutions: 

```python
>>> # Visualizing the main parameter specifications in the problem instance
>>> problem

FORMULATION: R2
 
MODEL SELECTION COMPUTED:  
     Lambda fixed
     Path
     Stability selection
 
LAMBDA FIXED PARAMETERS: 
     numerical_method = DR
     rescaled lam : True
     threshold = 0.177
     lam = 0.1
     theoretical_lam = 0.1994
 
PATH PARAMETERS: 
     numerical_method : Path-Alg
     Npath = 40
     lamin = 0.013
     lamax = 1.0
 
STABILITY SELECTION PARAMETERS: 
     numerical_method : Path-Alg
     method : first
     B = 50
     q = 10
     percent_nS = 0.5
     threshold = 0.7
     lamin = 0.01
     Nlam = 50

>>> # Visualizing the solutions
>>> problem.solution

 LAMBDA FIXED : 
   Selected variables :  43    47    74    79    84    
   Running time :  0.094 s

 PATH COMPUTATION : 
   Running time :  0.221 s

 STABILITY SELECTION : 
   Selected variables :  43    47    74    79    84    
   Running time :  2.468 s

```

![Graphics plotted after calling problem.solution ](figures/_figure-concat.png)


In the present synthetic setup, one can also plot the ground truth $\beta$ solution of the problem:

```python
>>> print( list(numpy.nonzero(sol)) )
[43, 47, 74, 79, 84]
```

It is indeed the variables that have been selected with the solution threshold for a fixed lambda, and with stability selection. Let us also note that the running time is still very low in our example. 

Those remarks are comforting, but not surprising because in this example the noise is little and the number of variable is still small. 



# Acknowledgements

This work was supported by the Flatiron Institute and the Helmholtz Zentrum Munchen. 

# References



---
title: 'c-lasso - a package for constrained sparse and robust regression and classification in Python'
tags:
  - Python
  - Linear regression
  - Optimization
authors:
  - name: Léo Simpson
    affiliation: 1
  - name: Patrick L. Combettes
    affiliation: 2
  - name: Christian L. Müller
    affiliation: 3
affiliations:
 - name: TUM  
   index: 1
 - name: NC State
   index: 2
 - name: LMU
   index: 3
date: 13 August 2020
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:

---

# Summary
c-lasso is a Python package that enables sparse and robust linear
regression and classification with linear equality constraints. 


The forward model is assumed to be:

$$
y = X \beta + \sigma \epsilon \qquad \textrm{s.t.} \qquad C\beta=0
$$

Here, $X$ and $y$ are given outcome and predictor data. The vector $y$ can be continuous (for regression) or binary (for classification). $C$ is a general constraint matrix. The vector $\beta$ comprises the unknown coefficients $\epsilon$ an unknown noise and $\sigma$ an unknown scale.

Depending on the prior we assume on those unknown variables, this forward model can lead to different types of estimation. Our package can solve six of those : four regression-type and two classification-type formulations. Those are all variants of the standard formulation "*R1*" : 

$$
    \arg \min_{\beta \in \mathbb{R}^d} \left\lVert X\beta - y \right\rVert^2 + \lambda \left\lVert \beta\right\rVert_1 \qquad s.t. \qquad  C\beta = 0
$$


### Formulations
- ***R1* Standard constrained Lasso regression**: 
This is the standard Lasso problem with linear equality constraints on the $\beta$ vector. 
The objective function combines Least-Squares for model fitting with $L_1$ penalty for sparsity.   

- ***R2* Contrained sparse Huber regression**: 
This regression problem uses the [Huber loss](https://en.wikipedia.org/wiki/Huber_loss) as objective function 
for robust model fitting with $L_1$ and linear equality constraints on the $\beta$ vector. The parameter $\rho$ is set to $1.345$ by default [@Aigner:1976]

- ***R3* Contrained scaled Lasso regression**: 
This formulation is similar to *R1* but allows for joint estimation of the (constrained) $\beta$ vector and 
the standard deviation $\sigma$ in a concomitant fashion [@Combettes:2020.1; @Combettes:2020.2].
This is the default problem formulation in c-lasso.

- ***R4* Contrained sparse Huber regression with concomitant scale estimation**: 
This formulation combines *R2* and *R3* to allow robust joint estimation of the (constrained) $\beta$ vector and 
the scale $\sigma$ in a concomitant fashion [@Combettes:2020.1; @Combettes:2020.2].

- ***C1* Contrained sparse classification with Square Hinge loss**: 
This formulation is similar to *R1* but adapted for classification tasks using the Square Hinge loss with (constrained) sparse $\beta$ vector estimation [@Lee:2013].

- ***C2* Contrained sparse classification with Huberized Square Hinge loss**:        
This formulation is similar to *C1* but uses the Huberized Square Hinge loss for robust classification with (constrained) sparse $\beta$ vector estimation [@Rosset:2007].


### Model selections
Different models are implemented together with the optimization schemes, to overcome the difficulty of choosing the penalization free parameter $\lambda$. 

- *Fixed Lambda* : This approach is simply letting the user choose the parameter $\lambda$, or to choose $l \in [0,1]$ such that $\lambda = l\times \lambda_{\max}$. 
The default value is a scale-dependent tuning parameter that has been proposed in [Combettes:2020.2] and derived in [@Shi:2016].

- *Path Computation* :The package also leaves the possibility to us to compute the solution for a range of $\lambda$ parameters in an interval $[\lambda_{\min}, \lambda_{\max}]$. It can be done using *Path-Alg* or warm-start with any other optimization scheme. 

[comment]: <> (This can be done much faster than by computing separately the solution for each $\lambda$ of the grid, by using the Path-alg algorithm. One can also use warm starts : starting with $\beta_0 = 0$ for $\lambda_0 = \lambda_{\max}$, and then iteratvely compute $\beta_{k+1}$ using one of the optimization schemes with $\lambda = \lambda_{k+1} := \lambda_{k} - \epsilon$ and with a warm start set to $\beta_{k}$. )

- *Cross Validation* : Then one can use a model selection, to choose the appropriate penalisation. This can be done by using k-fold cross validation to find the best $\lambda \in [\lambda_{\min}, \lambda_{\max}]$ with or without "one-standard-error rule" [@Hastie:2009].

- *Stability Selection* : Another variable selection model than can be used is stability selection [@Lin:2014; @Meinshausen:2010; Combettes:2020.2].



# Statement of need 
The package handles several estimators for inferring location and scale, including the constrained Lasso, the constrained scaled Lasso, and sparse Huber M-estimation with linear equality constraints Several algorithmic strategies, including path and proximal splitting algorithms, are implemented to solve the underlying convex optimization problems. We also include two model selection strategies for determining the sparsity of the model parameters: k-fold cross-validation and stability selection. This package is intended to fill the gap between popular python tools such as `scikit-learn` which <em>cannot</em> solve sparse constrained problems and general-purpose optimization solvers such as `cvx` that do not scale well for the considered problems or are inaccurate. We show several use cases of the package, including an application of sparse log-contrast regression tasks for compositional microbiome data. We also highlight the seamless integration of the solver into `R` via the `reticulate` package. 




# Basic workflow
Here is a basic example that shows how to run c-lasso on synthetic data.

c-lasso is available on pip, one can install it using ```pip install c_lasso```. Then on python, to import the package, one should use ```import classo```

Let us now begin the tutorial. 
Firstly, let us generate a dataset using
the routine ```random_data``` included in the c-lasso package, that allows you to generate instances using normally distributed data.

```python
>>> n,d,d_nonzero,k,sigma =100,100,5,1,0.5
>>> (X,C,y),sol = random_data(n,d,d_nonzero,k,sigma,zerosum=True, seed = 123 )
>>> list(numpy.nonzero(sol))
[43, 47, 74, 79, 84]
```
This code snippet generates randomly the vectors $\beta \in R^d$ , $X \in R^{n\times d}$ , $C \in R^{k\times d}$ (here it is the all-one vector instead because of the input ```zerosum```), and $y \in R^n$ normally distributed with respect to the model $C\beta=0$, $y-X\beta \sim N(0,I_n\sigma^2)$ and $\beta$ has only d_nonzero non-null componant (which are plot above).


Then, let us define a ```classo_problem``` instance with the generated dataset in order to formulate the optimization problem we want to solve. 

```python
# to define a c-lasso problem instance with default setting : 
>>> problem  = classo_problem(X,y,C)  
# to change the formulation of the problem :
>>> problem.formulation.huber  = True
>>> problem.formulation.concomitant = False
>>> problem.formulation.rho = 1.5  
# to add the computation for a fixed lambda :
>>> problem.model_selection.LAMfixed = True 
# to set lambda to 0.1*lambdamax : 
>>> problem.model_selection.LAMfixedparameters.rescaled_lam = True
>>> problem.model_selection.LAMfixedparameters.lam = 0.1 
# to add the computation of the lambda-path : 
>>> problem.model_selection.PATH = True 
# to solve our optimization problem : 
>>> problem.solve() 
```


Here, we have modified the [formulation](###formulations) of the problem in order to use *R2*, with $\rho=1.5$. 
We have chosen the following [model selections](###model-selections) : *Fixed Lambda* with $\lambda = 0.1\lambda_{\max}$ ; *Path computation* and *Stability Selection* which is computed by default. 
Then, those problems are solved using the recommanded optimization scheme on each model according to the formulation and the size of the parameter $\lambda$

Finally, one can visualize the solutions and see the running time, and the name of the selected variables by calling the instance ```problem.solution```. Note that by calling directly the instance ```problem``` one could also visualize the main parameters of the optimization problems one is solving. In our case, the running time is in the order of 0.1sec for the fixed lambda and path computation, but vary from 2sec to 4sec for the stability selection computation.  


![Graphics plotted after calling ```problem.solution``` ](figures/_figure-concat.png)




# ReferencesFolder containing previous paper versions
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior:
1. Go to '...'
2. Click on '....'
3. Scroll down to '....'
4. See error

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Desktop (please complete the following information):**
 - OS: [e.g. iOS]
 - Browser [e.g. chrome, safari]
 - Version [e.g. 22]

**Smartphone (please complete the following information):**
 - Device: [e.g. iPhone6]
 - OS: [e.g. iOS8.1]
 - Browser [e.g. stock browser, safari]
 - Version [e.g. 22]

**Additional context**
Add any other context about the problem here.
<img src="https://i.imgur.com/2nGwlux.png" alt="c-lasso" height="150" align="right"/>

# Numerical benchmarks for c-lasso 

We provide numerical benchmarks for the [c-lasso package](https://c-lasso.readthedocs.io/en/latest/) in comparison to [cvxpy](https://www.cvxpy.org). 
We report run times, achieved minimum function values (with the function value found by the path algorithm as baseline), and constraint satisfaction quality of the zero-sum constraint.  

## Table of Contents

* [Benchmark setup](#installation)
* [Results](#getting-started)
* [Optimization schemes](#optimization-schemes)

* [References](#references)


##  Benchmark set-up

###  Tested Regression and classification problems

Here, we consider the following problem formulations (see [here](../README.md) for detailed infos):

#### [R1] Standard constrained Lasso regression:             

<img src="https://latex.codecogs.com/gif.latex?\arg\min_{\beta\in&space;R^d}&space;||&space;X\beta-y&space;||^2&space;&plus;&space;\lambda&space;||\beta||_1&space;\qquad\mbox{s.t.}\qquad&space;C\beta=0" />

#### [R2] Contrained sparse Huber regression:                   

<img src="https://latex.codecogs.com/gif.latex?\arg\min_{\beta\in&space;R^d}&space;h_{\rho}(X\beta-y&space;)&space;&plus;&space;\lambda&space;||\beta||_1&space;\qquad\mbox{s.t.}\qquad&space;C\beta=0" />

#### [C1] Contrained sparse classification with Square Hinge loss (L1-Regularized Square-Hinge SVM): 

<img src="https://latex.codecogs.com/gif.latex?\arg&space;\min_{\beta&space;\in&space;\mathbb{R}^d}&space;\sum_{i=1}^n&space;l(y_i&space;x_i^\top&space;\beta)&space;&plus;&space;\lambda&space;\left\lVert&space;\beta\right\rVert_1&space;\qquad&space;s.t.&space;\qquad&space;C\beta&space;=&space;0" title="\arg \min_{\beta \in \mathbb{R}^d} \sum_{i=1}^n l(y_i x_i \beta) + \lambda \left\lVert \beta\right\rVert_1 \qquad s.t. \qquad C\beta = 0" />

where the x<sub>i</sub> are the rows of X and l is defined as:

<img src="https://latex.codecogs.com/gif.latex?l(r)&space;=&space;\begin{cases}&space;(1-r)^2&space;&&space;if&space;\quad&space;r&space;\leq&space;1&space;\\&space;0&space;&if&space;\quad&space;r&space;\geq&space;1&space;\end{cases}" title="l(r) = \begin{cases} (1-r)^2 & if \quad r \leq 1 \\ 0 &if \quad r \geq 1 \end{cases}" />

###  Synthetic data generation and test problem set-ups 

We use the `random_data` function in c-lasso to generate X and y. We use the standard `zeroSum` constraint. We vary the number of samples n and dimensionionality d of the problems. The regularization parameter is fixed to &lambda;=0.1. This setting does not favor the path algorithm. The reported performance is thus rather a lower bound on the actual speed-up. Since, for most model selection schemes, the computation of the entire solution path is required, the path algorithm formulation would be even more preferable then.  

## Results

The running times of the micro-benchmark has been computed using Python 3.9.1 on a laptop `MacBook Air`, operating on macOS high Sierra with the processor `1.8 GHz Intel Core i5`, with memory of `8 Gb 1600 MHz DDR3`.

#### R1

Run times for R1. 

![Run times on R1](./output/bm-R1-times.png)

Achieved minimum function values on R1. We observe considerable inaccuracies in cvx solutions.

![Achieved function values on C1](./output/bm-R1-losses.png)

Satistifaction of the zero-sum constraint.

![Satistifaction of the zero-sum constraint](./output/bm-R1-constraint.png)


#### R2

Run times for R2. 

![Run times on R2](./output/bm-R2-times.png)

Achieved minimum function values on R2. 

![Achieved function values on R2](./output/bm-R2-losses.png)

Satistifaction of the zero-sum constraint

![Satistifaction of the zero-sum constraint](./output/bm-R2-constraint.png)

#### C1

Run times for C1. 

![Run times on C1](./output/bm-C1-times.png)

Achieved minimum function values on C1. 

![Achieved function values on C1](./output/bm-C1-losses.png)

Satistifaction of the zero-sum constraint

![Satistifaction of the zero-sum constraint](./output/bm-C1-constraint.png)



## Optimization schemes

We consider the following schemes in the benchmark.

### Path algorithms (Path-Alg, pa) 
This is the default algorithm for non-concomitant problems [R1,R3,C1,C2]. 
The algorithm uses the fact that the solution path along &lambda; is piecewise-
affine (as shown, e.g., in [1]). When Least-Squares is used as objective function,
we derive a novel efficient procedure that allows us to also derive the 
solution for the concomitant problem [R2] along the path with little extra computational overhead.

### Projected primal-dual splitting method (pds):
This algorithm is derived from [2] and belongs to the class of 
proximal splitting algorithms. It extends the classical Forward-Backward (FB) 
(aka proximal gradient descent) algorithm to handle an additional linear equality constraint
via projection. In the absence of a linear constraint, the method reduces to FB.
This method can solve problem [R1]. 

### Douglas-Rachford-type splitting method (dr)
This algorithm is the most general algorithm and can solve all regression problems 
[R1-R4]. It is based on Doulgas Rachford splitting in a higher-dimensional product space.
It makes use of the proximity operators of the perspective of the LS objective (see [4,5])


### Operator splitting conic solver (SCS) in CVX (cvx)
For external comparison, we use cvx and its underlying conic solver (scs). For more info, see [6].



## References 

* [1] B. R. Gaines, J. Kim, and H. Zhou, [Algorithms for Fitting the Constrained Lasso](https://www.tandfonline.com/doi/abs/10.1080/10618600.2018.1473777?journalCode=ucgs20), J. Comput. Graph. Stat., vol. 27, no. 4, pp. 861–871, 2018.

* [2] L. Briceno-Arias and S.L. Rivera, [A Projected Primal–Dual Method for Solving Constrained Monotone Inclusions](https://link.springer.com/article/10.1007/s10957-018-1430-2?shared-article-renderer), J. Optim. Theory Appl., vol. 180, Issue 3, March 2019.

* [3] S. Rosset and J. Zhu, [Piecewise linear regularized solution paths](https://projecteuclid.org/euclid.aos/1185303996), Ann. Stat., vol. 35, no. 3, pp. 1012–1030, 2007.

* [4] P. L. Combettes and C. L. Müller, [Perspective M-estimation via proximal decomposition](https://arxiv.org/abs/1805.06098), Electronic Journal of Statistics, 2020, [Journal version](https://projecteuclid.org/euclid.ejs/1578452535) 

* [5] P. L. Combettes and C. L. Müller, [Regression models for compositional data: General log-contrast formulations, proximal optimization, and microbiome data applications](https://arxiv.org/abs/1903.01050), Statistics in Bioscience, 2020.

* [6] B. O’Donoghue, E. Chu, N. Parikh, and S. Boyd. [Conic optimization via operator splitting and homogeneous self-dual embedding.](https://link.springer.com/article/10.1007/s10957-016-0892-3), Journal of Optimization Theory and Applications 169, no. 3 (2016): 1042-1068.


