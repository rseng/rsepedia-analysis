# DCISolver - Dynamic Control of Infeasibility Solver

| **Documentation** | **CI** | **Coverage** | **Release** | **DOI** |
|:-----------------:|:------:|:------------:|:-----------:|:-------:|
| [![docs-stable][docs-stable-img]][docs-stable-url] [![docs-dev][docs-dev-img]][docs-dev-url] | [![build-ci][build-ci-img]][build-ci-url] | [![codecov][codecov-img]][codecov-url] | [![release][release-img]][release-url] | [![doi][doi-img]][doi-url] |

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://JuliaSmoothOptimizers.github.io/DCISolver.jl/stable
[docs-dev-img]: https://img.shields.io/badge/docs-dev-purple.svg
[docs-dev-url]: https://JuliaSmoothOptimizers.github.io/DCISolver.jl/dev
[build-ci-img]: https://github.com/JuliaSmoothOptimizers/DCISolver.jl/workflows/CI/badge.svg?branch=main
[build-ci-url]: https://github.com/JuliaSmoothOptimizers/DCISolver.jl/actions
[codecov-img]: https://codecov.io/gh/JuliaSmoothOptimizers/DCISolver.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/JuliaSmoothOptimizers/DCISolver.jl
[release-img]: https://img.shields.io/github/v/release/JuliaSmoothOptimizers/DCISolver.jl.svg?style=flat-square
[release-url]: https://github.com/JuliaSmoothOptimizers/DCISolver.jl/releases
[doi-img]: https://joss.theoj.org/papers/10.21105/joss.03991/status.svg
[doi-url]: https://doi.org/10.21105/joss.03991

DCI is a solver for equality-constrained nonlinear problems, i.e.,
optimization problems of the form

    min f(x)     s.t.     c(x) = 0.

It uses other JuliaSmoothOptimizers packages for development.
In particular, [NLPModels.jl](https://github.com/JuliaSmoothOptimizers/NLPModels.jl) is used for defining the problem, and [SolverCore](https://github.com/JuliaSmoothOptimizers/SolverCore.jl) for the output.
It uses [LDLFactorizations.jl](https://github.com/JuliaSmoothOptimizers/LDLFactorizations.jl) by default to compute the factorization in the tangent step. Follow [HSL.jl](https://github.com/JuliaSmoothOptimizers/HSL.jl)'s `MA57` installation for an alternative.
The feasibility steps are factorization-free and use iterative methods from [Krylov.jl](https://github.com/JuliaSmoothOptimizers/Krylov.jl)

## References

> Bielschowsky, R. H., & Gomes, F. A.
> Dynamic control of infeasibility in equality constrained optimization.
> SIAM Journal on Optimization, 19(3), 1299-1325 (2008).
> [10.1007/s10589-020-00201-2](https://doi.org/10.1007/s10589-020-00201-2)

## How to Cite

If you use DCISolver.jl in your work, please cite using the format given in [CITATION.bib](https://github.com/JuliaSmoothOptimizers/DCISolver.jl/blob/main/CITATION.bib).

## Installation

1. [LDLFactorizations.jl](https://github.com/JuliaSmoothOptimizers/LDLFactorizations.jl) is used by default. Follow [HSL.jl](https://github.com/JuliaSmoothOptimizers/HSL.jl)'s `MA57` installation for an alternative.
2. `pkg> add DCISolver`

## Example

```julia
using DCISolver, ADNLPModels

# Rosenbrock
nlp = ADNLPModel(x -> 100 * (x[2] - x[1]^2)^2 + (x[1] - 1)^2, [-1.2; 1.0])
stats = dci(nlp)

# Constrained
nlp = ADNLPModel(x -> 100 * (x[2] - x[1]^2)^2 + (x[1] - 1)^2, [-1.2; 1.0],
                 x->[x[1] * x[2] - 1], [0.0], [0.0])
stats = dci(nlp)
```

# Bug reports and discussions

If you think you found a bug, feel free to open an [issue](https://github.com/JuliaSmoothOptimizers/DCISolver.jl/issues).
Focused suggestions and requests can also be opened as issues. Before opening a pull request, start an issue or a discussion on the topic, please.

If you want to ask a question not suited for a bug report, feel free to start a discussion [here](https://github.com/JuliaSmoothOptimizers/Organization/discussions). This forum is for general discussion about this repository and the [JuliaSmoothOptimizers](https://github.com/JuliaSmoothOptimizers), so questions about any of our packages are welcome.
# Reference
​
## Contents
​
```@contents
Pages = ["reference.md"]
```
​
## Index
​
```@index
Pages = ["reference.md"]
```
​
```@autodocs
Modules = [DCISolver]
```# Advanced-usage of DCI

## Contents

```@contents
Pages = ["fine-tuneDCI.md"]
```

The main function exported by this package is the function `dci` whose basic usage has been illustrated previously.
It is also possible to fine-tune the parameters used in the implementation in two different ways.

## Examples

DCISolver.jl exports the function `dci`:
```
   dci(nlp :: AbstractNLPModel)
   dci(nlp :: AbstractNLPModel, x :: AbstractVector)
   dci(nlp :: AbstractNLPModel, meta :: MetaDCI, x :: AbstractVector)
   dci(nlp :: AbstractNLPModel, meta :: MetaDCI, workspace :: DCIWorkspace)
```
where `MetaDCI` is a structure handling all the parameters used in the algorithm, and `DCIWorkspace` pre-allocates all the memory used during the iterative process.

It is therefore possible to either call `dci(nlp, x, kwargs...)` and the keywords arguments are passed to the `MetaDCI` constructor or build an instance of `MetaDCI` directly.

```@example ex1
using ADNLPModels, DCISolver

nlp = ADNLPModel(
  x -> 100 * (x[2] - x[1]^2)^2 + (x[1] - 1)^2, 
  [-1.2; 1.0],
  x->[x[1] * x[2] - 1], 
  [0.0], [0.0],
  name = "Rosenbrock with x₁x₂=1"
)

#The alternative would be:
stats = dci(
  nlp, nlp.meta.x0, 
  max_time = 600., 
  linear_solver = :ldlfact, 
  TR_compute_step = :TR_lsmr
)
```
The alternative would be:
```@example ex1
meta = DCISolver.MetaDCI(
  nlp.meta.x0, nlp.meta.y0, 
  max_time = 600., 
  linear_solver = :ldlfact, 
  TR_compute_step = :TR_lsmr
)
stats = dci(nlp, meta, nlp.meta.x0)
```

The `DCIWorkspace` allows to reuse the same memory if one would re-solve a problem of the same dimension.

```@example ex1
workspace = DCISolver.DCIWorkspace(nlp, meta, nlp.meta.x0)
stats = dci(nlp, meta, workspace)
worspace.x0 .= ones(2) # change the initial guess, and resolve
stats = dci(nlp, meta, workspace)
```

## List of possible options

Find below a list of the main options of `dci`.

### Tolerances on the problem
We use `ϵ = atol + rtol * dualnorm`.
```
| Parameters           | Type          | Default      | Description                                    |
| -------------------- | ------------- | ------------ | ---------------------------------------------- |
| atol                 | AbstractFloat | 1e-5         | absolute tolerance.                            |
| rtol                 | AbstractFloat | 1e-5         | relative tolerance.                            |
| ctol                 | AbstractFloat | 1e-5         | feasibility tolerance.                         |
| unbounded_threshold  | AbstractFloat | -1e5         | below this threshold the problem is unbounded. |
| max_eval             | Integer       | 50000        | maximum number of cons + obj evaluations.      |
| max_time             | AbstractFloat | 120.         | maximum number of seconds.                     |
| max_iter             | Integer       | 500          | maximum number of iterations.                  |
| max_iter_normal_step | Integer       | typemax(Int) | maximum number of iterations in normal step.   |
```

### Compute Lagrange multipliers
```
| Parameters  | Type        | Default                                         | Description                                           |
| ----------- | ----------- | ----------------------------------------------- | ----------------------------------------------------- |
| comp_λ      | Symbol      | :cgls                                           | eval(comp_λ) is used to compute Lagrange multipliers. |
| λ_struct    | comp_λ_cgls | comp_λ_cgls(length(x0), length(y0), typeof(x0)) | companion structure of `comp_λ`.                      |
```

### Tangent step
```
| Parameters    | Type          | Default  | Description                                                                                               |
| ------------- | ------------- | -------- | --------------------------------------------------------------------------------------------------------- |
| linear_solver | Symbol        | :ldlfact | Solver for the factorization. options: :ma57.                                                             | 
| decrease_γ    | AbstractFloat | 0.1      | Regularization for the factorization: reduce γ if possible, > √eps(T), between tangent steps.             |
| increase_γ    | AbstractFloat | 100.0    | Regularization for the factorization: up γ if possible, < 1/√eps(T), during the factorization.            |
| δmin          | AbstractFloat | √eps(T)  | Regularization for the factorization: smallest value of δ used for the regularization.                    |
| tan_Δ         | AbstractFloat | 1.0      | Tangent step trust-region parameters: initial trust-region radius.                                        |
| tan_η₁        | AbstractFloat | 1e-2     | Tangent step trust-region parameters: decrease the trust-region radius when Ared/Pred < η₁.               |
| tan_η₂        | AbstractFloat | 0.75     | Tangent step trust-region parameters: increase the trust-region radius when Ared/Pred > η₂.               |
| tan_σ₁        | AbstractFloat | 0.25     | Tangent step trust-region parameters: decrease coefficient of the trust-region radius.                    |
| tan_σ₂        | AbstractFloat | 2.0      | Tangent step trust-region parameters: increase coefficient of the trust-region radius.                    |
| tan_small_d   | AbstractFloat | eps(T)   | Tangent step trust-region parameters: ||d|| is too small.                                                 |
| increase_Δtg  | AbstractFloat | 10.0     | Tangent step trust-region parameters: increase if possible, < 1 / √eps(T), the Δtg between tangent steps. |
```
### Normal step
```
| Parameters             | Type                                    | Default                                            | Description                                                                                               |
| ---------------------- | --------------------------------------- | -------------------------------------------------- | --------------------------------------------------------------------------------------------------------- |
| feas_step              | Symbol                                  | :feasibility_step                                  | Normal step                                                                                               |
| feas_η₁                | AbstractFloat                           | 1e-3                                               | Feasibility step: decrease the trust-region radius when Ared/Pred < η₁.                                   |
| feas_η₂                | AbstractFloat                           | 0.66                                               | Feasibility step: increase the trust-region radius when Ared/Pred > η₂.                                   |
| feas_σ₁                | AbstractFloat                           | 0.25                                               | Feasibility step: decrease coefficient of the trust-region radius.                                        |
| feas_σ₂                | AbstractFloat                           | 2.0                                                | Feasibility step: increase coefficient of the trust-region radius.                                        |
| feas_Δ₀                | AbstractFloat                           | 1.0                                                | Feasibility step: initial radius.                                                                         |
| feas_expected_decrease | AbstractFloat                           | 0.95                                               | Feasibility step: bad steps are when ‖c(z)‖ / ‖c(x)‖ >feas_expected_decrease.                             |
| bad_steps_lim          | Integer                                 | 3                                                  | Feasibility step: consecutive bad steps before using a second order step.                                 |
| TR_compute_step        | Symbol                                  | :TR_lsmr                                           | Compute the direction in feasibility step: options: :TR_dogleg.                                           |
| TR_compute_step_struct | Union{TR_lsmr_struct, TR_dogleg_struct} | TR_lsmr_struct(length(x0), length(y0), typeof(x0)) | Compute the direction in feasibility step: options: TR_dogleg_struct(length(x0), length(y0), typeof(x0)). |
```
### Parameters updating ρ (or redefine the function `compute_ρ`)
```
| Parameters  | Type          | Default | Description                                                     |
| ----------- | ------------- | ------- | --------------------------------------------------------------- |
| compρ_p1    | AbstractFloat | 0.75    | update ρ as `ρ = max(min(ngp, p1) * ρmax, ϵ)`.                  |
| compρ_p2    | AbstractFloat | 0.90    | update ρ as `ρ = primalnorm * p2` if not sufficiently feasible. |
| ρbar        | AbstractFloat | 2.0     | radius of the larger cylinder is `ρbar * ρ`.                    |
```
The computation of ρ can also be modified by importing `compute_ρ(dualnorm, primalnorm, norm∇fx, ρmax, ϵ, iter, meta::MetaDCI)`
# DCISolver - Dynamic Control of Infeasibility Solver

DCI is a solver for equality-constrained nonlinear problems, i.e.,
optimization problems of the form

```math
    \min_x \ f(x) \quad \text{s.t.} \quad  c(x) = 0,
```

based on the paper

> Bielschowsky, R. H., & Gomes, F. A.
> Dynamic control of infeasibility in equality constrained optimization.
> SIAM Journal on Optimization, 19(3), 1299-1325 (2008).
> [10.1007/s10589-020-00201-2](https://doi.org/10.1007/s10589-020-00201-2)

`DCISolver` is a JuliaSmoothOptimizers-compliant solver. It takes an [`AbstractNLPModel`](https://github.com/JuliaSmoothOptimizers/NLPModels.jl) as an input and returns a [`GenericExecutionStats`](https://github.com/JuliaSmoothOptimizers/SolverCore.jl/blob/16fc349908f46634f2c9acdddddb009b23634b71/src/stats.jl#L60).

We refer to [juliasmoothoptimizers.github.io](https://juliasmoothoptimizers.github.io) for tutorials on the NLPModel API. This framework allows the usage of models from Ampl (using [AmplNLReader.jl](https://github.com/JuliaSmoothOptimizers/AmplNLReader.jl)), CUTEst (using [CUTEst.jl](https://github.com/JuliaSmoothOptimizers/CUTEst.jl)), JuMP (using [NLPModelsJuMP.jl](https://github.com/JuliaSmoothOptimizers/NLPModelsJuMP.jl)), PDE-constrained optimization problems (using [PDENLPModels.jl](https://github.com/JuliaSmoothOptimizers/PDENLPModels.jl)) and models defined with automatic differentiation (using [ADNLPModels.jl](https://github.com/JuliaSmoothOptimizers/ADNLPModels.jl)).

## Installation

`DCISolver` is a registered package. To install this package, open the Julia REPL (i.e., execute the julia binary), type `]` to enter package mode, and install `DCISolver` as follows
```
add DCISolver
```

The DCI algorithm is an iterative method that has the flavor of a projected gradient algorithm and could be characterized as
a relaxed feasible point method with dynamic control of infeasibility. It is a combination of two steps: a tangent step and a feasibility step.
It uses [LDLFactorizations.jl](https://github.com/JuliaSmoothOptimizers/LDLFactorizations.jl) by default to compute the factorization in the tangent step. Follow [HSL.jl](https://github.com/JuliaSmoothOptimizers/HSL.jl)'s `MA57` installation for an alternative.
The feasibility steps are factorization-free and use iterative methods from [Krylov.jl](https://github.com/JuliaSmoothOptimizers/Krylov.jl).

## Example

We consider in this example the minization of the Rosenbrock function over an equality constraint.
```math
    \min_x \ 100 * (x₂ - x₁²)² + (x₁ - 1)² \quad \text{s.t.} \quad  x₁x₂=1,
```
The problem is modeled using `ADNLPModels.jl` with `[-1.2; 1.0]` as default initial point, and then solved using `dci`.
```@example
using DCISolver, ADNLPModels, Logging
nlp = ADNLPModel(
  x -> 100 * (x[2] - x[1]^2)^2 + (x[1] - 1)^2, 
  [-1.2; 1.0],
  x -> [x[1] * x[2] - 1], 
  [0.0], [0.0],
  name = "Rosenbrock with x₁x₂=1"
)
stats = with_logger(NullLogger()) do
  dci(nlp)
end

println(stats)
```

# Bug reports and discussions

If you think you found a bug, feel free to open an [issue](https://github.com/JuliaSmoothOptimizers/DCISolver.jl/issues).
Focused suggestions and requests can also be opened as issues. Before opening a pull request, start an issue or a discussion on the topic, please.

If you want to ask a question not suited for a bug report, feel free to start a discussion [here](https://github.com/JuliaSmoothOptimizers/Organization/discussions). This forum is for general discussion about this repository and the [JuliaSmoothOptimizers](https://github.com/JuliaSmoothOptimizers), so questions about any of our packages are welcome.
# Benchmarks

## CUTEst benchmark

With a JSO-compliant solver, such as DCI, we can run the solver on a set of problems, explore the results, and compare to other JSO-compliant solvers using specialized benchmark tools. 
We are following here the tutorial in [SolverBenchmark.jl](https://juliasmoothoptimizers.github.io/SolverBenchmark.jl/v0.3/tutorial/) to run benchmarks on JSO-compliant solvers.
``` @example ex1
using CUTEst
```

To test the implementation of DCI, we use the package [CUTEst.jl](https://github.com/JuliaSmoothOptimizers/CUTEst.jl), which implements `CUTEstModel` an instance of `AbstractNLPModel`. 

``` @example ex1
using SolverBenchmark
```

Let us select equality-constrained problems from CUTEst with a maximum of 100 variables or constraints. After removing problems with fixed variables, examples with a constant objective, and infeasibility residuals.

``` @example ex1
_pnames = CUTEst.select(
  max_var = 100, 
  min_con = 1, 
  max_con = 100, 
  only_free_var = true, 
  only_equ_con = true, 
  objtype = 3:6
)

#Remove all the problems ending by NE as Ipopt cannot handle them.
pnamesNE = _pnames[findall(x->occursin(r"NE\b", x), _pnames)]
pnames = setdiff(_pnames, pnamesNE)
cutest_problems = (CUTEstModel(p) for p in pnames)

length(cutest_problems) # number of problems
```

We compare here DCISolver with [Ipopt](https://link.springer.com/article/10.1007/s10107-004-0559-y) (Wächter, A., & Biegler, L. T. (2006). On the implementation of an interior-point filter line-search algorithm for large-scale nonlinear programming. Mathematical programming, 106(1), 25-57.), via the [NLPModelsIpopt.jl](https://github.com/JuliaSmoothOptimizers/NLPModelsIpopt.jl) thin wrapper, with DCISolver on a subset of CUTEst problems.

``` @example ex1
using DCISolver, NLPModelsIpopt
```
 To make stopping conditions comparable, we set `Ipopt`'s parameters `dual_inf_tol=Inf`, `constr_viol_tol=Inf` and `compl_inf_tol=Inf` to disable additional stopping conditions related to those tolerances, `acceptable_iter=0` to disable the search for an acceptable point.

``` @example ex1
#Same time limit for all the solvers
max_time = 1200. #20 minutes
tol = 1e-5

solvers = Dict(
  :ipopt => nlp -> ipopt(
    nlp,
    print_level = 0,
    dual_inf_tol = Inf,
    constr_viol_tol = Inf,
    compl_inf_tol = Inf,
    acceptable_iter = 0,
    max_cpu_time = max_time,
    x0 = nlp.meta.x0,
    tol = tol,
  ),
  :dcildl => nlp -> dci(
    nlp,
    nlp.meta.x0,
    linear_solver = :ldlfact,
    max_time = max_time,
    max_iter = typemax(Int64),
    max_eval = typemax(Int64),
    atol = tol,
    ctol = tol,
    rtol = tol,
  ),
)

stats = bmark_solvers(solvers, cutest_problems)
```
The function `bmark_solvers` return a `Dict` of `DataFrames` with detailed information on the execution. This output can be saved in a data file.
``` @example ex1
using JLD2
@save "ipopt_dcildl_$(string(length(pnames))).jld2" stats
```
The result of the benchmark can be explored via tables,
``` @example ex1
pretty_stats(stats[:dcildl])
```
or it can also be used to make performance profiles.
``` @example ex1
using Plots
gr()

legend = Dict(
  :neval_obj => "number of f evals", 
  :neval_cons => "number of c evals", 
  :neval_grad => "number of ∇f evals", 
  :neval_jac => "number of ∇c evals", 
  :neval_jprod => "number of ∇c*v evals", 
  :neval_jtprod  => "number of ∇cᵀ*v evals", 
  :neval_hess  => "number of ∇²f evals", 
  :elapsed_time => "elapsed time"
)
perf_title(col) = "Performance profile on CUTEst w.r.t. $(string(legend[col]))"

styles = [:solid,:dash,:dot,:dashdot] #[:auto, :solid, :dash, :dot, :dashdot, :dashdotdot]

function print_pp_column(col::Symbol, stats)
  
  ϵ = minimum(minimum(filter(x -> x > 0, df[!, col])) for df in values(stats))
  first_order(df) = df.status .== :first_order
  unbounded(df) = df.status .== :unbounded
  solved(df) = first_order(df) .| unbounded(df)
  cost(df) = (max.(df[!, col], ϵ) + .!solved(df) .* Inf)

  p = performance_profile(
    stats, 
    cost, 
    title=perf_title(col), 
    legend=:bottomright, 
    linestyles=styles
  )
end

print_pp_column(:elapsed_time, stats) # with respect to time
```

``` @example ex1
print_pp_column(:neval_jac, stats) # with respect to number of jacobian evaluations
```

## CUTEst benchmark with Knitro

In this second part, we present the result of a similar benchmark with a maximum of 10000 variables and constraints (82 problems), and including the solver [`KNITRO`](https://link.springer.com/chapter/10.1007/0-387-30065-1_4) (Byrd, R. H., Nocedal, J., & Waltz, R. A. (2006). K nitro: An integrated package for nonlinear optimization. In Large-scale nonlinear optimization (pp. 35-59). Springer, Boston, MA.) via [`NLPModelsKnitro.jl`](https://github.com/JuliaSmoothOptimizers/NLPModelsKnitro.jl). The script is included in [/benchmark/script10000_knitro.jl)](https://github.com/JuliaSmoothOptimizers/DCISolver.jl/benchmark/script10000_knitro.jl). We report here a performance profile with respect
to the elapsed time to solve the problems and to the sum of evaluations of objective and constrain functions, see [/benchmark/figures.jl)](https://github.com/JuliaSmoothOptimizers/DCISolver.jl/benchmark/figures.jl) for the code generating the profile wall.

![](./assets/ipopt_knitro_dcildl_82.png)
