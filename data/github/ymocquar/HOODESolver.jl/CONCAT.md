# HOODESolver.jl

*A Julia package for solving numerically highly-oscillatory ODE problems*

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ymocquar.github.io/HOODESolver.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ymocquar.github.io/HOODESolver.jl/dev)
[![Build Status](https://github.com/ymocquar/HOODESolver.jl/workflows/CI/badge.svg)](https://github.com/ymocquar/HOODESolver.jl/actions)
[![Coverage](https://codecov.io/gh/ymocquar/HOODESolver.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ymocquar/HOODESolver.jl)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![status](https://joss.theoj.org/papers/816cd9b9f4815a78a08ede5e46fd2978/status.svg)](https://joss.theoj.org/papers/816cd9b9f4815a78a08ede5e46fd2978)

## Installation

HOODESolver.jl is a registered package and can be installed using the Julia package manager. From the Julia REPL, enter Pkg mode (by pressing `]`)

```julia
julia>]
(@v1.5) pkg> add HOODESolver
```

## Usage

The following is an example with the system of Hénon-Heiles. Please see the [documentation](https://ymocquar.github.io/HOODESolver.jl/stable/) for further usage, tutorials, and api reference.

```julia
using HOODESolver
using Plots

epsilon= 0.0001

A = [ 0 0 1 0 ; 
      0 0 0 0 ; 
     -1 0 0 0 ; 
      0 0 0 0 ]

f1 = LinearHOODEOperator( epsilon, A)

f2 = (u,p,t) ->  [ 0, u[4], 2*u[1]*u[2], -u[2] - u[1]^2 + u[2]^2 ] 

tspan = (0.0, 3.0)

u0 = [0.55, 0.12, 0.03, 0.89]
prob = SplitODEProblem(f1, f2, u0, tspan)
```

solve the defined problem

```julia
sol = solve(prob, HOODEAB()) 
plot(sol) 
```
![](docs/src/img/example.png)

For support with using HOODESolver.jl, please open an [issue](https://github.com/ymocquar/HOODESolver.jl/issues/new/) describing the problem and steps to reproduce it.

## How to Contribute

Here's an outline of the workflow you should use if you want to make contributions to this package.

1. Fork this repository
2. Make a new branch on your fork, named after whatever changes you'll be making
3. Apply your code changes to the branch on your fork
4. When you're done, submit a PR to `HOODESolver.jl` to merge your fork into master branch.


This package is licensed under the MIT Expat license. See [LICENSE](LICENSE) for more information.
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: 'bug'
assignees: '@pnavaro'
---

A clear and concise description of what the bug is.

Give us your environment

```julia
using HOODESolver
println(versioninfo())
using Pkg; 
println(Pkg.status("HOODESolver"))
```

A clear and concise description of what you expected to happen.
A minimal example of the input text should be supplied,
together with the output.

Thanks
# HOODESolver.jl

The objective of this Julia package is to valorize the recent developments carried out within [INRIA team MINGuS](https://team.inria.fr/mingus/) on Uniformly Accurate numerical methods (UA) for highly oscillating problems. We propose to solve the following equation 

$$\frac{d u(t)}{dt} = \frac{1}{\varepsilon} A u(t) + f(t, u(t)), \qquad u(t=t_{start})=u_{in}, \qquad \varepsilon\in ]0, 1], \qquad (1)$$

with 
-  $u : t\in [t_{start}, t_{end}] \mapsto u(t)\in \mathbb{R}^n, \quad t_{start}, t_{end}\in \mathbb{R}$, 
-  $u_{in}\in \mathbb{R}^n$, 
-  $A\in {\mathcal{M}}_{n,n}(\mathbb{R})$ is such that $\tau \mapsto \exp(\tau A)$ is $2 \pi$-periodic,  
-  $f : (t, u) \in  \mathbb{R}\times \mathbb{R}^n \mapsto \mathbb{R}^n$.

The purpose here is to write an explanatory documentation of the *Julia* package containing the two-scale method (see [Chartier2020](@cite), [Chartier2015](@cite) and [Crouseilles2013](@cite). This package is inspired by the Differential Equations package [SciML](https://diffeq.sciml.ai/dev/index.html).


## References

```@bibliography
```

## Index

```@index
```
# Future work

- precision on parameter.
- exponential runge-kutta method

## The "classical" Runge-Kutta method (order 4) adapted to the exponential (not yet implemented)

Notations : 
- We denote by $G$ the function which transforms $\hat{U}$ to $\hat{f}$, so $G(\hat{U}) = \hat{f}$.
- We denote by $S_{t_0}^{t_1}(t_2,\ell)$ the intégral $S_{t_0}^{t_1}(t_2,\ell) = \int_{t_0}^{t_1} e^{- i \ell (t_2 - s)/\varepsilon} ds = ( i \varepsilon / \ell) ( e^{- i \ell (t_2 - t_1)/\varepsilon}-e^{- i \ell (t_2 - t_0)/\varepsilon})$


Here are the calculations

-  $u_{1,\ell} = \hat{U}_{n, \ell}$
-  $u_{2,\ell} = e^{- i \ell h_n /(2 \varepsilon)}\hat{U}_{n, \ell} + S_0^{h_n /2} ( h_n /2,\ell ) G_{\ell}(u_1)$
-  $u_{3,\ell} = e^{- i \ell h_n /(2 \varepsilon)}\hat{U}_{n, \ell} + S_0^{h_n /2} ( h_n /2,\ell )  G_{\ell}(u_2)$

-  $u_{4,\ell} = e^{- i \ell h_n /(2\varepsilon)}u_{2,\ell} + S_0^{h_n/2} ( h_n/2,\ell )[ 2 G_{\ell}(u_3)-G_{\ell}(u_1)]$ (see (28) of [Cox2002](@cite), with $c=-i \ell h_n /\varepsilon$)   

From (29) of [Cox2002](@cite), with $c=-i \ell h_n /\varepsilon$, we have

$$\hat{U}_{n+1, \ell} = e^{- i \ell h_n /\varepsilon}\hat{U}_{n, \ell} +  G_{\ell}(u_1) [-4+i \ell h_n /\varepsilon + e^{-i \ell h_n /\varepsilon}(4+3i \ell h_n /\varepsilon+(i \ell h_n /\varepsilon)^2]\\+ (2 G_{\ell}(u_2) + G_{\ell}(u_3) )[-2-i \ell h_n /\varepsilon+e^{-i \ell h_n /\varepsilon}(2-i \ell h_n /\varepsilon)]\\ + G_{\ell}(u_4)[-4+3i \ell h_n /\varepsilon -(i \ell h_n /\varepsilon)^2 + e^{-i \ell h_n /\varepsilon}(4+i \ell h_n /\varepsilon)]/(h_n^2 (i \ell h_n /\varepsilon)^3)$$

# Numerical method

## Two-scale formulation 

First, rewrite equation (1) using the variable change $w(t)=\exp(-(t-t_{start})A/\varepsilon) u(t)$ to obtain

$$\frac{d w(t)}{dt} = F\Big(\frac{t-t_{start}}{\varepsilon}, w(t) \Big), \;\;\; w(t_{start})=u_{in}, \;\; \varepsilon\in ]0, 1],$$

where the function $F$ is expressed from the data of the original problem (1)

$$F\Big( \frac{s}{\varepsilon}, w \Big) = \exp(-sA/\varepsilon) \; f( \exp(sA/\varepsilon), \; w).$$

We then introduce the function $U(t, \tau), \tau\in [0, 2 \pi]$ such that $U(t, \tau=(t-t_{start})/\varepsilon) = w(t)$. The two-scale function is then the solution of the following equation

$$\frac{\partial U}{\partial t} + \frac{1}{\varepsilon} \frac{\partial U}{\partial \tau} =  F( \tau, U), \;\;\; U(t=t_{start}, \tau)=\Phi(\tau), \;\; \varepsilon\in ]0, 1], \;\;\;\;\;\;\;\;\;\; (2)$$

where $\Phi$ is a function checking $\Phi(\tau=0)=u_{0}$ chosen so that the $U$ solution of (2) is smooth (see [Chartier2015](@cite) and [Chartier2020](@cite).

## Discretization  

The numerical method is based on a discretization of equation (2). In the direction $\tau$, a spectral method is used, while for the time $t$, an exponential Adams-Bashforth method allows to build a high order method (see [Chartier2020](@cite)). The initialization is based on a "butterfly" technique (going back and forth around the initial time).

### Initialization 

Let r be the order of the method $AB_r$.\
Let $\Delta t$ the time step, for $i \in \{r, -(r-1), \ldots, r-1, r\}$,
we note $u_i = u(t_{start}+i \Delta t)$.\
Let $r'$ be the orders of the intermediate AB methods we will use.\
If $u_{k}$ is known with a precision of ${\mathcal O}(\Delta t^{r'+1})$, and for $r' \geq 2, u_{k-1}, \ldots, u_{k-r'+1}$ are known with a precision of ${\mathcal O}(\Delta t^{r'})$ then we can calculate $u_{k+1}$ with a precision of ${\mathcal O}(\Delta t^{r'+1})$ with the method $AB_{r'}$.\
Similarly, if $u_{k}$ is known with a precision of ${\mathcal O}(\Delta t^{r'+1})$, and for $r' \geq 2, u_{k+1}, \ldots, u_{k+r'-1}$ are known with a precision of ${\mathcal O}(\Delta t^{r'})$ then we can calculate $u_{k-1}$ with a precision of ${\mathcal O}(\Delta t^{r'+1})$ with the method $AB_{r'}$.

### Algorithm

- With the method $AB_1$, from $u_0$ we calculate $u_{-1}$ with a precision of ${\mathcal O}(\Delta t^2)$
- With the method $AB_2$, starting from $u_{0}$ and $u_{-1}$, we calculate $u_{1}$ with a precision of ${\mathcal O}(\Delta t^3)$
- For $r' = $3 to $r' = r$.
    - For $k=1$ to $k=r'-1$
        - With the method $AB_{r'-1}$, from $u_{1-k}, u_{2-k}, \ldots,u_{r'-1-k}$, we calculate $u_{-k}$ with a precision of ${\mathcal O}(\Delta t^{r'})$
    - For $k=1$ to $k=r'-1$
         - With the method $AB_{r'}$, from $u_{k-1}, u_{k-2}, \ldots,u_{k-r'}$, we calculate $u_{k}$ with a precision of ${\mathcal O}(\Delta t^{r'+1})$

At the end of this algorithm, the values $u_0, u_1, \ldots u_{r-1}$ are known with a precision of ${\mathcal O}(\Delta t^{r+1})$, we can launch the algorithm $AB_r$. Note that this technique requires that the function $f$ has to be smooth enough (namely $f\in {\mathcal C}^r([t_{start} - r \Delta t, t_{end}])$).  


### The Adams-Bashforth Method
For an Adams-Bashforth method of order $r$ in time and spectral $\tau$, we first introduce a mesh in the $\tau$ direction.\
$\tau_{\ell} = \ell \Delta \tau, \ell = 0, \ldots, N_{\tau}-1$. Where $N_{\tau}$ is the number of points of discretization. If we apply the Fourier transform to the two-scale equation, we obtain 

$$\frac{\partial \hat{U}_\ell}{\partial t} + \frac{i\ell}{\varepsilon}\hat{U}_\ell = \hat{F}_\ell(t), \;\; \ell=-N_\tau/2, \dots, N_\tau/2-1,$$

with

$$U(t, \tau_k) = \sum_{\ell=-N_{\tau}/2}^{N_{\tau}/2-1} \hat{U}_{\ell}(t) e^{i\ell k 2\pi/N_{\tau}} \;\;\; \text{ and } F(\tau_k, U(t, \tau_k)) = \sum_{\ell=-N_{\tau}/2}^{N_{\tau}/2-1} \hat{F}_{\ell}(t) e^{i\ell k 2\pi/N_\tau}.$$

and the inverse (discrete) Fourier transform formulae

$$\hat{U}_{\ell}(t) = \frac{1}{N_{\tau}}\sum_{k=0}^{N_{\tau}-1} U(t, \tau_k) e^{-i\ell k 2\pi/N_{\tau}}\;\;\; \text{ and } \hat{F}_{\ell}(t) = \frac{1}{N_{\tau}}\sum_{k=0}^{N_{\tau}-1} F(\tau_k, U(t, \tau_k))e^{-i\ell k 2\pi/N_{\tau}}.$$

If we wish to calculate $\hat{F}_{\ell}$ from $\hat{U}_{\ell}$ we have the following formula 

$$F(\tau_k,U(t, \tau_k)) = e^{-\tau_k A}f(e^{\tau_k A}U(t, \tau_k)) \;\;$$

from which the Fourier transform is calculated in $\tau$ from the discrete Fourier transform formulas above.

Now, given a time step $\Delta t>0$ and a discretization in time $t_n=n\Delta t$, we can write the following Duhamel formula ($n\geq 0$)

$$\hat{U}_{\ell}(t_{n+1}) 
= e^{-i\ell\Delta t/\varepsilon}\hat{U}_{\ell}(t_{n})  + \int_0^{\Delta t} e^{-i\ell(\Delta t -s)/\varepsilon} \hat{F}_\ell(t_n+s)ds.$$

Thus, to obtain a $(r+1)$ scheme, we can 
approaches the function $\hat{F}_\ell(t_n+s)$ by the Lagrange polynomial of order $r$ interpolator at points $t_{n-j}, j=0, \dots, r$. This polynomial is written 

$$\hat{F}_\ell(t_n+s) \approx \sum_{k=0}^r \Big(\Pi_{j=0, j\neq k}^r \frac{s+j \Delta t}{(j-k)\Delta t} \Big) \hat{F}_\ell(t_n-t_j), \;\; n\geq 0.$$

Thus, from this approximation, we integrate exactly, which requires the following formulas

$$p^{[r]}_{\ell, j} = \int_0^{\Delta t}e^{-i\ell(\Delta t -s)/\varepsilon}\Big( \Pi_{j=0, j\neq k}^r \frac{s+j \Delta t}{(j-k)\Delta t}\Big) ds,$$

for each $j$ and $\ell$ such that $0\leq j\leq r, \; \ell=-N_\tau/2, \dots, N_\tau/2-1$. These coefficients $p^{[r]}_{\ell, j}$ can be pre-calculated and stored once and for all. Thus, the schema is finally written

$$\hat{U}_{\ell}^{n+1}= e^{-i\ell\Delta t/\varepsilon}\hat{U}_{\ell}^n + \sum_{j=0}^r p^{[r]}_{\ell, j} \hat{F}_\ell^{n-j},$$

with $\hat{U}_{\ell}^n \approx \hat{U}_{\ell}(t_n)$ and $\hat{F}_\ell^{n-j}\approx \hat{F}_\ell(t_{n-j})$. 

We can verify that the truncation error in this schema is  ${\mathcal O}(\Delta t^{r+1})$, once the initial values  $\hat{U}_\ell^1, \dots, \hat{U}_\ell^r$ have been calculated.  

## Non-homogeneous case $f(t, u)$

Here we consider the case where $f$ depends on the variable $t$.

$$\frac{d u(t)}{dt} = \frac{1}{\varepsilon} A u(t) + f(t, u(t)), \;\;\; u(t=t_{start})=u_{in}, \;\; \varepsilon\in ]0, 1] \;\;\;\; (3)$$

The non-homogeneous case (3) falls under (1), by entering the variable 
$\theta : t \in [t_{start}, t_{end}] \mapsto \theta(t) = t\in \mathbb{R}$ 
which allows us to reformulate the non-homogeneous case (3) into a homogeneous problem of the form (1).
Indeed, we rephrase (3) as follows

$$\begin{aligned}
\frac{d u(t) }{dt}  & = \frac{1}{\varepsilon} Au(t) + f(\theta(t), u(t)), \\
\frac{d \theta(t) }{dt}  & = 1
\end{aligned}\;\;\;\;(4)$$

with the initial condition  $u(t_{start})=u_{in}, \theta(t_{start})=t_{start}$. Thus, the problem (4) is rewritten into an equation\
satisfied by $y: t\in [t_{start}, t_{end}] \mapsto y(t) =(u(t), \theta(t))\in \mathbb{R}^{n+1}$

$$\frac{d y}{dt} = \frac{1}{\varepsilon} \tilde{A} y + g(y), \;\; y(t_{start})=(u_{in}, t_{start}),$$

with $\tilde{A}\in{\mathcal M}_{n+1, n+1}(\mathbb{R})$

$$\tilde{A}=
\left(
\begin{array}{cccc}
   &    &     & 0 \\
   & A &     & 0 \\
   &    &    & 0 \\
0 & 0 & 0 &  0 
\end{array}
\right) \;\;\;\; \text{ and } \;\;\;\;
g(y)=g(u, \theta) = \left(
\begin{array}{cccccc}
f(\theta, u) \\
1
\end{array}
\right) \in \mathbb{R}^{n+1}.$$
# Charged particle Example

A system of charged particles under the effect of an external electro-magnetic field is considered to be 
$(E(t, x), B(t, x))\in \mathbb{R}^6$.\
Particles are dynamically described by their position 
$x(t)\in\mathbb{R}^3$ and their speed $v(t)\in\mathbb{R}^3$. We'll index by $i$ the $i$-th component of a vector.
Newton's equations applied to a particle can be written as

$$\begin{aligned}
\frac{d x(t) }{dt}&= v(t) \\
\frac{d v(t) }{dt}&= \frac{e}{m} \left[E(t, x(t)) + v(t)\times B(t, x(t))\right]. 
\end{aligned}$$

We will assume that the magnetic field is written $B(t, x)=(0, 0, 1)^T$ and under a certain scaling, we consider the following equation

$$\begin{aligned}
\frac{d x_1(t) }{dt} &= \frac{1}{\varepsilon}v_1(t) \\
\frac{d x_2(t) }{dt} &= \frac{1}{\varepsilon} v_2(t) \\
\frac{d x_3(t) }{dt} &= v_3(t) \\
\frac{d v_1(t) }{dt} &= E_1(t, x(t)) + \frac{1}{\varepsilon}v_2(t)\\
\frac{d v_2(t) }{dt} &= E_2(t, x(t)) - \frac{1}{\varepsilon}v_1(t)\\
\frac{d v_3(t) }{dt} &= E_3(t, x(t)) 
\end{aligned}$$

which is rewritten as follows

$$\frac{d u(t) }{dt}= \frac{1}{\varepsilon}A u(t) + F(t, u(t)),$$

where the unknown vector $u(t)=(x(t), v(t))\in\mathbb{R}^6$, $A$ is a square matrix of size $6\times 6$
and $F$ is a function with a value in $\mathbb{R}^6$. $A$ and $F$ are given by

$$A=
\left(
\begin{array}{cccccc}
0 & 0 & 0 & 1 & 0 & 0 \\
0 & 0 & 0 & 0 & 1 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 1 & 0 \\
0 & 0 & 0 & -1 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 
\end{array}
\right) \;\;\;\; \text{ and } \;\;\;\;
F(t, u(t)) = \left(
\begin{array}{cccccc}
0 \\
0\\
u_6(t)\\
E_1(t, u_1(t), u_2(t), u_3(t))\\
E_2(t,  u_1(t), u_2(t), u_3(t)\\
E_3(t,  u_1(t), u_2(t), u_3(t)
\end{array}
\right).$$

We can consider the following $E=(E_1, E_2, E_3)$ function

$$E(t, x) =
\left(
\begin{array}{ccc}
\cos(x_1/2)\sin(x_2)\sin(x_3)/2\\
\sin(x_1/2)\cos(x_2)\sin(x_3)\\
\sin(x_1/2)\sin(x_2)\cos(x_3)
\end{array}
\right)$$

```@setup 12
using HOODESolver
using Plots
```

```@example 12
epsilon = 0.05

A = [0 0 0  1 0 0; 
     0 0 0  0 1 0;
     0 0 0  0 0 0; 
     0 0 0  0 1 0; 
     0 0 0 -1 0 0; 
     0 0 0  0 0 0]

f1 = LinearHOODEOperator( epsilon, A)

function f2(u, p, t)
    s1, c1 = sincos(u[1]/2)
    s2, c2 = sincos(u[2])
    s3, c3 = sincos(u[3])
    return [0, 0, u[6], c1*s2*s3/2, s1*c2*s3, s1*s2*c3]
end

tspan = (0.0, 1.0)
u0 = [1.0, 1.5, -0.5, 0, -1.2, 0.8]
prob = SplitODEProblem(f1, f2, u0, tspan)
sol = solve(prob, HOODEAB() )
plot(sol)
```


```@example 12
plot(sol,vars=(1,2,3))
```
# Functions

```@autodocs
Modules = [HOODESolver]
Order   = [:function]
```
# Types

```@autodocs
Modules = [HOODESolver]
Order   = [:type]
```
# SplitODEProblem

The `SplitODEProblem` type from package [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/types/split_ode_types/) offers a interface for problems similar to the ones we are trying to solve.

Consider the Henon-Heiles system we can use one of the solvers dedicated to split problems, see [Split ODE Solvers](https://diffeq.sciml.ai/stable/solvers/split_ode_solve/#split_ode_solve):

```@example 4
using Plots
```

```@example 4
using DifferentialEquations

epsilon = 0.002
A = [ 0 0 1 0 ;
      0 0 0 0 ;
     -1 0 0 0 ;
      0 0 0 0 ]

f1 = DiffEqArrayOperator( A ./ epsilon)

function f2(du, u, p, t)
    du[1] = 0
    du[2] = u[4]
    du[3] = 2*u[1]*u[2]
    du[4] = -u[2] - u[1]^2 + u[2]^2 
end

tspan = (0.0, 0.1)

u0 = [0.55, 0.12, 0.03, 0.89]

prob1 = SplitODEProblem(f1, f2, u0, tspan);
sol1 = solve(prob1, ETDRK4(), dt=0.001);
nothing # hide
```
With our method we need to give the value of `epsilon`. In some case, you can get this value from `f1`.
You can pass a `DiffEqArrayOperator` as argument to the problem but the method is not always valid 
so we define a new type called `LinearHOODEOperator`:

```@example 4
using HOODESolver

f1 = LinearHOODEOperator(epsilon, A)
prob2 = SplitODEProblem(f1, f2, u0, tspan)
sol2 = solve(prob2, HOODEAB(), dt=0.01)

plot(sol1, vars=[3], label="EDTRK4")
plot!(sol2, vars=[3], label="HOODEAB")
plot!(sol2.t, getindex.(sol2.u, 3), m=:o)
```
# Quickstart

## Hénon-Heiles Example

We consider the system of Hénon-Heiles satisfied by $u(t)=(u_1, u_2, u_3, u_4)(t)$.

$$\frac{d u }{dt} = \frac{1}{\varepsilon} Au + f(u), \;\;\; u(t_{start})=u_{inn}\in\mathbb{R}^4,$$

where $A$ and $f$ are selected as follows

$$A=
\left(
\begin{array}{cccc}
0 & 0 & 1 & 0  \\
0 & 0 & 0 & 0  \\
-1 & 0 & 0 & 0  \\
0 & 0 & 0 & 0  
\end{array}
\right) \;\;\;\; \text{ and } \;\;\;\;
f(u) = \left(
\begin{array}{cccc}
0 \\
u_4\\
-2 u_1 u_2\\
-u_2-u_1^2+u_2^2
\end{array}
\right).$$

one chooses for example, $\varepsilon=0.001$
and $u_{in} = (0.12, 0.12, 0.12, 0.12)$

## Input parameters

The input arguments use the same format as the ODE package. 

Thus, first of all, we must define the arguments necessary to construct the problem (1), namely

- the $f$ function (in the form *Julia*) 
- the initial condition $u_{in}$.
- the initial time $t_{start}$ and final time $t_{end}$. 
- the second parameter of the 
- the $A$ matrix 
-  $\varepsilon \in ]0, 1]$

```@example 1
using HOODESolver, Plots

epsilon= 0.0001

A = [ 0 0 1 0 ; 
      0 0 0 0 ; 
     -1 0 0 0 ; 
      0 0 0 0 ]

f1 = LinearHOODEOperator( epsilon, A)

f2 = (u,p,t) ->  [ 0, u[4], 2*u[1]*u[2], -u[2] - u[1]^2 + u[2]^2 ] 

tspan = (0.0, 3.0)

u0 = [0.55, 0.12, 0.03, 0.89]

prob = SplitODEProblem(f1, f2, u0, tspan);
nothing # hide
```

From the `prob` problem, we can now switch to its numerical resolution. 

To do this, the numerical parameters are defined 
- the number of time slots $N_t$ which defines the time step $\Delta t = \frac{t_{\text{end}}-t_{start}}{N_t}$ but you can set the value of time step `dt` in the `solve` call.
- the $r$ order of the method 
- the number of $N_\tau$ points in the $\tau$ direction... 
- the order of preparation $q$ of the initial condition 

The default settings are : $N_t=100$, $r=4$, $N_\tau=32$ and $q=r+2=6$
To solve the problem with the default parameters, just call the `solve` command with the problem already defined as parameter

```@example 1     
sol = solve(prob, HOODEAB(), dt=0.1);
nothing # hide
```
Which is equivalent to `HOODEAB( order=4, nb_tau=32 )`


### Exhaustive definition of the parameters

- `prob` : problem defined by `SplitODEProblem` 
- `nb_tau=32` : $N_{\tau}$
- `order=4` : order $r$ of the method
- `order_prep=order+2` : order of preparation of initial data
- `dense=true` : indicates whether or not to keep the data from the fourier transform, if `dense=false`, processing is faster but interpolation can no longer be done.
- `nb_t=100` : $N_t$
- `getprecision=dense` : indicates whether the accuracy is calculated, the method used to calculate the accuracy multiplies the processing time by 2.
- `verbose=100` : trace level, if `verbose=0` then nothing is displayed.
- `par_u0` : If we have to make several calls to `solve` with the same initial data and in the same order, we can pass in parameter the already calculated data.
- `p_coef` : table with the coefficients of the Adams-Bashforth method. This array can be used to optimize several calls with the same parameters.

## Exit arguments

As an output, a structure of type `ODESolution` from [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/basics/solution/).
This structure can be seen as a function of t, it can also be seen as an array. Example:

```@repl 1     
sol.prob
sol.alg
t=2.541451547
sol(t)
sol[end]
sol(3.0)
```

To view the result, you can also use Plot, for example

```@example 1     
plot(sol) 
```

## Linear non-homogeneous case 

The following non-homogeneous linear system is considered to be satisfied by $u(t)=(u_1, u_2, u_3, u_4)(t)$

```math
\frac{d u }{dt} = \frac{1}{\varepsilon} Au + f(t, u), \;\;\; u(0)=u_in\in\mathbb{R}^4,
```

where ``A`` and ``f`` are selected as follows

```math
A=
\left(
\begin{array}{cccc}
0 & 0 & 1 & 0  \\
0 & 0 & 0 & 0  \\
-1 & 0 & 0 & 0  \\
0 & 0 & 0 & 0  
\end{array}
\right) \;\;\; \text{ and } \;\;\;
f(t, u) = Bu +\alpha t +\beta \;\; \text{ with  } \;\;
B\in {\mathcal M}_{4, 4}(\mathbb{R}), \alpha, \beta \in \mathbb{R}^4,
```

``B, \alpha, \beta`` are chosen randomly.

We wish to obtain a high precision, so we will use BigFloat real
numbers, they are encoded on 256 bits by default which gives a
precision bound of about ``2^{-256}. \approx 10^{-77}``.  At the
end, we compare a calculated result with an exact result. For this, you must
use the dedicated [`HOODEProblem`](@ref) type:

```@example 2
using HOODESolver
using Plots
using Random
```

```@example 2

rng = MersenneTwister(1111)

A=[0 0 1 0 ; 0 0 0 0 ; -1 0 0 0 ; 0 0 0 0]
B = 2rand(rng, BigFloat, 4, 4) - ones(BigFloat, 4, 4)
alpha = 2rand(rng, BigFloat, 4) - ones(BigFloat, 4)
beta = 2rand(rng, BigFloat, 4) - ones(BigFloat, 4)
u0 = [big"0.5", big"-0.123", big"0.8", big"0.7"]
t_start=big"0.0"
t_end=big"1.0"
epsilon=big"0.017"

fct = (u,p,t)-> B*u + t*p[1] +p[2]

prob = HOODEProblem(fct, u0, (t_start,t_end), (alpha, beta), A, epsilon, B)

sol = solve(prob, nb_t=10000, order=8)
sol.absprec
```

```@example 2
t=big"0.9756534187771"
sol(t)-getexactsol(sol.par_u0.parphi, u0, t)
```

```@example 2
plot(sol.t,sol.u_tr)
```

### Calculation of the exact solution

This involves calculating the exact solution $u(t)$ of the following equation at the instant $t$

```math
\frac{d u }{dt} = \frac{1}{\varepsilon} Au + Bu +\alpha t +\beta, \;\;\; u(t_{start})=u_{in}\in\mathbb{R}^4\text{, } A \text{ and }B \text{ are defined above }
```

Let

```math
\begin{aligned}
    M &= \frac{1}{\varepsilon} A + B\\
    C &= e^{-t_{start} M}u_{in} +M^{-1} e^{-t_{start} M} (t_{start}\alpha+\beta)+ M^{-2} e^{-t_{start} M} \alpha\\
    C_t &= -M^{-1} e^{-t M} (t\alpha+\beta)-M^{-2} e^{-t M} \alpha\\
    u(t) &= e^{t M} ( C + C_t)
\end{aligned}
```

Which, translated into Julia language, gives the code of the function `getexactsol` : 

```julia
function getexactsol(par::PreparePhi, u0, t)
    @assert !ismissing(par.matrix_B) "The debug matrix is not defined"
    sparse_A = par.sparse_Ap[1:(end-1),1:(end-1)]
    m = (1/par.epsilon)*sparse_A+par.matrix_B
    t_start = par.t_0
    if ismissing(par.paramfct)
        return exp((t-t0)*m)*u0
    end
    a, b = par.paramfct
    mm1 = m^(-1)
    mm2 = mm1^2
    e_t0 = exp(-t_start*m)
    C = e_t0*u0 + mm1*e_t0*(t_start*a+b)+mm2*e_t0*a
    e_inv = exp(-t*m)
    e = exp(t*m)
    C_t = -mm1*e_inv*(t*a+b)-mm2*e_inv*a
    return e*C+e*C_t
end
```

## Accuracy of the result according to the time interval

### Linear problem

From a problem of the previous type, as long as we can calculate the exact solution, it is possible to know exactly what the error is.
The initialization data being

```jl
using HOODESolver

A = [0 0 1 0 ; 0 0 0 0 ; -1 0 0 0 ; 0 0 0 0]

u0 = BigFloat.([-34//100, 78//100, 67//100, -56//10])

B = BigFloat.([12//100 -78//100 91//100 34//100
    -45//100 56//100 3//100 54//100
    -67//100 09//100 18//100 89//100
    -91//100 -56//100 11//100 -56//100])

alpha =  BigFloat.([12//100, -98//100, 45//100, 26//100])

beta =  BigFloat.([-4//100, 48//100, 23//100, -87//100])

epsilon = 0.015

t_end = big"1.0"

fct = (u,p,t)-> B*u + t*p[1] +p[2]

prob = HOODEProblem(fct,u0, (big"0.0",t_end), (alpha, beta), A, epsilon, B)
```

Note that the floats are coded on 512 bits.\
By varying $\Delta t$ from $10^{-2}$ to $5.10^{-6}$ (i.e. `nb_t` from `100` to `204800`) on a logarithmic scale, for odd orders from 3 to 17 we get these errors

#### Precision of the result with ε = 0.015

![](img/error_order.png)

Now with the same initial data, order being setted to 6, and $\varepsilon = 0.15, 0.015, \ldots, 1.5\times 10^{-7}$.\
Here floats are coded on 256 bits.

#### Precision of the result with order = 6

![](img/error_epsilon.png)

### Problem with Hénon-Heiles function

```jl
u0=BigFloat.([90, -44, 83, 13]//100)
t_end = big"1.0"
epsilon=big"0.0017"
fct = u -> [0, u[4], -2u[1]*u[2], -u[2]-u[1]^2+u[2]^2]
A = [0 0 1 0; 0 0 0 0;-1 0 0 0; 0 0 0 0]
prob = HOODEProblem(fct, u0, (big"0.0",t_end), missing, A, epsilon)
```
The float are coded on 512 bits.

#### Precision of the result with ε = 0.0017

![](img/error_order_hh.png)

Now with the same initial data, order being setted to 6, and $\varepsilon = 0.19, 0.019, \ldots, 1.9\times 10^{-8}$.\
Here floats are coded on 256 bits.

#### Precision of the result with ordre = 6

![](img/error_epsilon_hh.png)
