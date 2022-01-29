# POMDPStressTesting.jl

[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://sisl.github.io/POMDPStressTesting.jl/dev/)
[![Build Status](https://github.com/sisl/POMDPStressTesting.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/sisl/POMDPStressTesting.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/sisl/POMDPStressTesting.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/sisl/POMDPStressTesting.jl)


Adaptive stress testing of black-box systems, implemented within the [POMDPs.jl](https://github.com/JuliaPOMDP/POMDPs.jl) ecosystem.

See the [documentation](https://sisl.github.io/POMDPStressTesting.jl/dev/) for more details.

# Citation

If you use this package for research purposes, please cite the following:

[![status](https://joss.theoj.org/papers/04dc39ea89e90938727d789a2e402b0b/status.svg)](https://joss.theoj.org/papers/04dc39ea89e90938727d789a2e402b0b)

```
@article{moss2021pomdpstresstesting,
  title = {{POMDPStressTesting.jl}: Adaptive Stress Testing for Black-Box Systems},
  author = {Robert J. Moss},
  journal = {Journal of Open Source Software},
  year = {2021},
  volume = {6},
  number = {60},
  pages = {2749},
  doi = {10.21105/joss.02749}
}
```

# Interface
To stress test a new system, the user has to define the `GrayBox` and `BlackBox` interface outlined in [`src/GrayBox.jl`](https://github.com/mossr/POMDPStressTesting.jl/blob/master/src/GrayBox.jl) and [`src/BlackBox.jl`](https://github.com/mossr/POMDPStressTesting.jl/blob/master/src/BlackBox.jl).

### GrayBox Interface
The `GrayBox` simulator and environment interface includes:
* `GrayBox.Simulation` type to hold simulation variables
* `GrayBox.environment(sim::Simulation)` to return the collection of environment distributions
* `GrayBox.transition!(sim::Simulation)` to transition the simulator, returning the log-likelihood

### BlackBox  Interface
The `BlackBox` system interface includes:
* `BlackBox.initialize!(sim::Simulation)` to initialize/reset the system under test
* `BlackBox.evaluate!(sim::Simulation)` to evaluate/execute the system under test
* `BlackBox.distance(sim::Simulation)` to return how close we are to an event
* `BlackBox.isevent(sim::Simulation)` to indicate if a failure event occurred
* `BlackBox.isterminal(sim::Simulation)` to indicate the simulation is in a terminal state

Functions ending with `!` may modify the `Simulation` object in place.


# Solvers
Several solvers are implemented.

#### Reinforcement learning solvers
* [`MCTSPWSolver`](https://github.com/mossr/POMDPStressTesting.jl/blob/master/src/solvers/mcts.jl)

#### Deep reinforcement learning solvers<sup>1</sup>
* [`TRPOSolver`](https://github.com/mossr/POMDPStressTesting.jl/blob/master/src/solvers/drl/trpo.jl)
* [`PPOSolver`](https://github.com/mossr/POMDPStressTesting.jl/blob/master/src/solvers/drl/ppo.jl)

#### Stochastic optimization solvers
* [`CEMSolver`](https://github.com/mossr/POMDPStressTesting.jl/blob/master/src/solvers/cem.jl)

#### Baseline solvers
* [`RandomSearchSolver`](https://github.com/mossr/POMDPStressTesting.jl/blob/master/src/solvers/random_search.jl)


# Example
[![Example Notebook](https://img.shields.io/badge/example-notebook-blue)](https://nbviewer.jupyter.org/github/sisl/POMDPStressTesting.jl/blob/master/notebooks/Walk1D.ipynb)

An example implementation of the AST interface is provided for the Walk1D problem:
* **Julia source**: [`test/Walk1D.jl`](https://github.com/mossr/POMDPStressTesting.jl/blob/master/test/Walk1D.jl)
* **Jupyter notebook**: [`Walk1D.ipynb`](https://nbviewer.jupyter.org/github/sisl/POMDPStressTesting.jl/blob/master/notebooks/Walk1D.ipynb)
* **Descriptive tutorial-style write-up**: [`walk1d.pdf`](./test/pdf/walk1d.pdf) (created using [TeX.jl](https://github.com/mossr/TeX.jl))

<!-- (https://github.com/mossr/POMDPStressTesting.jl/blob/master/test/walk1d.pdf) -->

<kbd>
<p align="center">
  <a href="./test/pdf/walk1d.pdf">
    <img src="./test/svg/walk1d.svg">
  </a>
</p>
</kbd>

<!-- With an accompanying notebook: [`Walk1D.ipynb`](https://github.com/mossr/POMDPStressTesting.jl/blob/master/notebooks/Walk1D.ipynb) -->

# Installation

Install the `POMDPStressTesting.jl` package via:
```julia
] add POMDPStressTesting
```

### Testing
To run the test suite, you can use the Julia package manager.
```julia
] test POMDPStressTesting
```

# Contributing
We welcome contributions! Please fork the repository and submit a new Pull Request.

---
Package maintained by Robert Moss: mossr@cs.stanford.edu

<sup>1</sup> TRPO and PPO thanks to [Shreyas Kowshik's](https://github.com/shreyas-kowshik/RL-baselines.jl) initial implementation.
---
title: 'POMDPStressTesting.jl: Adaptive Stress Testing for Black-Box Systems'
tags:
  - Julia
  - stress testing
  - black-box systems
  - POMDPs.jl
authors:
  - name: Robert J. Moss
    orcid: 0000-0003-2403-454X
    affiliation: 1
affiliations:
 - name: Stanford University
   index: 1
date: 26 August 2020
bibliography: paper.bib
header-includes: |
    \usepackage{listings}
---
\lstdefinelanguage{Julia}{
    keywords=[3]{initialize!, transition!, evaluate!, distance, isevent, isterminal, environment},
    keywords=[2]{Nothing, Tuple, Real, Bool, Simulation, BlackBox, GrayBox, Sampleable, Environment},
    keywords=[1]{function, abstract, type, end},
    sensitive=true,
    morecomment=[l]{\#},
    morecomment=[n]{\#=}{=\#},
    morestring=[s]{"}{"},
    morestring=[m]{'}{'},
    alsoletter=!?,
    literate={,}{{\color[HTML]{0F6FA3},}}1
             {\{}{{\color[HTML]{0F6FA3}\{}}1
             {\}}{{\color[HTML]{0F6FA3}\}}}1
}

\lstset{
    language         = Julia,
    backgroundcolor  = \color[HTML]{F2F2F2},
    basicstyle       = \small\ttfamily\color[HTML]{19177C},
    numberstyle      = \ttfamily\scriptsize\color[HTML]{7F7F7F},
    keywordstyle     = [1]{\bfseries\color[HTML]{1BA1EA}},
    keywordstyle     = [2]{\color[HTML]{0F6FA3}},
    keywordstyle     = [3]{\color[HTML]{0000FF}},
    stringstyle      = \color[HTML]{F5615C},
    commentstyle     = \color[HTML]{AAAAAA},
    rulecolor        = \color[HTML]{000000},
    frame=lines,
    xleftmargin=10pt,
    framexleftmargin=10pt,
    framextopmargin=4pt,
    framexbottommargin=4pt,
    tabsize=4,
    captionpos=b,
    breaklines=true,
    breakatwhitespace=false,
    showstringspaces=false,
    showspaces=false,
    showtabs=false,
    columns=fullflexible,
    keepspaces=true,
    numbers=none,
}


# Summary

\href{https://github.com/sisl/POMDPStressTesting.jl}{POMDPStressTesting.jl} is a package that uses reinforcement learning and stochastic optimization to find likely failures in black-box systems through a technique called adaptive stress testing [@ast].
Adaptive stress testing (AST) has been used to find failures in safety-critical systems such as aircraft collision avoidance systems [@ast_acasx], flight management systems [@ast_fms], and autonomous vehicles [@ast_av].
The POMDPStressTesting.jl package is written in Julia [@julia] and is part of the wider POMDPs.jl ecosystem [@pomdps_jl], which provides access to simulation tools, policies, visualizations, and---most importantly---solvers.
We provide different solver variants including online planning algorithms such as Monte Carlo tree search [@mcts] and deep reinforcement learning algorithms such as trust region policy optimization (TRPO) [@trpo] and proximal policy optimization (PPO) [@ppo].
Stochastic optimization solvers such as the cross-entropy method [@cem] are also available and random search is provided as a baseline.
Additional solvers can easily be added by adhering to the POMDPs.jl interface.

The AST formulation treats the falsification problem (i.e., finding failures) as a Markov decision process (MDP) with a reward function that uses a measure of distance to a failure event to guide the search towards failure.
The reward function also uses the state transition probabilities to guide towards \textit{likely} failures.
Reinforcement learning aims to maximize the discounted sum of expected rewards, therefore maximizing the sum of log-likelihoods is equivalent to maximizing the likelihood of a trajectory.
A gray-box simulation environment steps the simulation and outputs the state transition probabilities, and the black-box system under test is evaluated in the simulator and outputs an event indication and the real-valued distance metric (i.e., how close we are to failure).
To apply AST to a general black-box system, a user has to implement the following Julia interface:

\begin{lstlisting}[language=Julia]
# GrayBox simulator and environment
abstract type GrayBox.Simulation end
function GrayBox.environment(sim::Simulation)::GrayBox.Environment end
function GrayBox.transition!(sim::Simulation)::Real end

# BlackBox.interface(input::InputType)::OutputType
function BlackBox.initialize!(sim::Simulation)::Nothing end
function BlackBox.evaluate!(sim::Simulation)::Tuple{Real, Real, Bool} end
function BlackBox.distance(sim::Simulation)::Real end
function BlackBox.isevent(sim::Simulation)::Bool end
function BlackBox.isterminal(sim::Simulation)::Bool end
\end{lstlisting}

Our package builds off work originally done in the AdaptiveStressTesting.jl package [@ast], but POMDPStressTesting.jl adheres to the interface defined by POMDPs.jl and provides different action modes and solver types.
Related falsification tools (i.e. tools that do not include most-likely failure analysis) are \textsc{S-TaLiRo} [@staliro], Breach [@breach], and \textsc{FalStar} [@falstar].
These packages use a combination of optimization, path planning, and reinforcement learning techniques to solve the falsification problem.
The tool most closely related to POMDPStressTesting.jl is the AST Toolbox in Python [@ast_av], which wraps around the gym reinforcement learning environment [@gym].
The author has contributed to the AST Toolbox and found the need to create a similar package in pure Julia for better performance and to interface with the POMDPs.jl ecosystem.

# Statement of Need

Validating autonomous systems is a crucial requirement before their deployment into real-world environments.
Searching for likely failures using automated tools enable engineers to address potential problems during development.
Because many autonomous systems are in environments with rare failure events, it is especially important to incorporate likelihood of failure within the search to help inform the potential problem mitigation.
This tool provides a simple interface for general black-box systems to fit into the adaptive stress testing problem formulation and gain access to solvers.
Due to varying simulation environment complexities, random seeds can be used as the AST action when the user does not have direct access to the environmental probability distributions or when the environment is complex.
Alternatively, directly sampling from the distributions allows for finer control over the search.
The interface is designed to easily extend to other autonomous system applications and explicitly separating the simulation environment from the system under test allows for wider validation of complex black-box systems.



# Research and Industrial Usage

POMDPStressTesting.jl has been used to find likely failures in aircraft trajectory prediction systems [@ast_fms], which are flight-critical subsystems used to aid in-flight automation.
A developmental commercial flight management system was stress tested so the system engineers could mitigate potential issues before system deployment [@ast_fms].
In addition to traditional requirements-based testing for avionics certification [@do178c], this work is being used to find potential problems during development.
There is also ongoing research on the use of POMDPStressTesting.jl for assessing the risk of autonomous vehicles and determining failure scenarios of autonomous lunar rovers. 


# Acknowledgments

We acknowledge Ritchie Lee for his guidance and original work on adaptive stress testing and the AdaptiveStressTesting.jl package and Mark Koren, Xiaobai Ma, and Anthony Corso for their work on the AST Toolbox Python package and the CrossEntropyMethod.jl package.
We also acknowledge Shreyas Kowshik for his initial implementation of the TRPO and PPO algorithms in Julia.
We want to thank the Stanford Intelligent Systems Laboratory for their development of the POMDPs.jl ecosystem and the MCTS.jl package; particular thanks to Zachary Sunberg.
We also want to thank Mykel J. Kochenderfer for his support and research input and for advancing the Julia community.


# References
# POMDPStressTesting
[![Build Status](https://github.com/sisl/POMDPStressTesting.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/sisl/POMDPStressTesting.jl/actions/workflows/CI.yml) [![codecov](https://codecov.io/gh/sisl/POMDPStressTesting.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/sisl/POMDPStressTesting.jl)

This package is used to find likely failures in a black-box software system.
The package is integrated with the [POMDPs.jl](https://github.com/JuliaPOMDP/POMDPs.jl) ecosystem; giving access to solvers, policies, and visualizations (although no prior knowledge of the POMDPs.jl package is needed to use POMDPStressTesting.jl—see the [Guide](@ref guide)). It uses a technique called adaptive stress testing (AST)[^1] to find likely failures using a distance metric to a failure event and the likelihood of an environment sample to guide the search.

A POMDP is a [partially observable Markov decision process](https://en.wikipedia.org/wiki/Partially_observable_Markov_decision_process), which is a framework to define a sequential decision making problem where the true state is unobservable. In the context of this package, we use the POMDP acronym mainly to tie the package to the POMDPs.jl package, but the system that is stress tested can also be defined as a POMDP.

This package is intended to help developers stress test their systems before deployment into the real-world (see existing use cases for aircraft collision avoidance systems[^1] and aircraft trajectory prediction systems[^2]). It is also used for research purposes to expand on the AST concept by allowing additional solution methods to be explored and tested.

---

[^1] Ritchie Lee et al., *"Adaptive Stress Testing: Finding Likely Failure Events with Reinforcement Learning
"*, 2020. [https://arxiv.org/abs/1811.02188](https://arxiv.org/abs/1811.02188)

[^2] Robert J. Moss, Ritchie Lee, Nicholas Visser, Joachim Hochwarth, James G. Lopez, Mykel J. Kochenderfer, *"Adaptive Stress Testing of Trajectory Predictions in Flight Management Systems"*, DASC 2020. [https://arxiv.org/abs/2011.02559](https://arxiv.org/abs/2011.02559)

---

## Package Features
- Search for failures in a black-box [system](@ref system)
- Define probability distributions of your [simulation environment](@ref sim_env)
- Find likely system failures using a variety of [solvers](@ref solvers)
- Calculate and visualize [failure metrics](@ref metrics_visualizations)
- Replay found failures


## Contents
```@contents
Pages = [
    "installation.md",
    "guide.md",
    "solvers.md",
    "example.md",
    "contrib.md"
]
```

## [`BlackBox` System Definition](@id system)
A black-box system could be an external software executable, code written in Julia, or code written in another language.
The system is generally a sequential decision making system than can be stepped foward in time.
It is termed "black-box" because all we need is to be able to initialize it (using `BlackBox.initialize!`), evaluate or step the system forward in time (using `BlackBox.evaluate!`), and parse the output of the system to determine the distance metric (using `BlackBox.distance`), the failure event indication (using `BlackBox.isevent`), and whether the system is in a terminal state (using `BlackBox.isterminal`).
- See the [`BlackBox` interface](@ref blackbox_interface) for implementation details

The `BlackBox` system interface includes:
- `BlackBox.initialize!(sim::Simulation)` to initialize/reset the system under test
- `BlackBox.evaluate!(sim::Simulation)` to evaluate/execute the system under test
- `BlackBox.distance(sim::Simulation)` to return how close we are to an event
- `BlackBox.isevent(sim::Simulation)` to indicate if a failure event occurred
- `BlackBox.isterminal(sim::Simulation)` to indicate the simulation is in a terminal state


## [`GrayBox` Simulator/Environment Definition](@id sim_env)
The gray-box simulator and environment define the parameters of your simulation and the probability distributions governing your simulation environment. It is termed "gray-box" because we need access to the probability distributions of the environment in order to get the log-likelihood of a sample used by the simulator (which is ulimately used by the black-box system).
- See the [`GrayBox` interface](@ref graybox_interface) for implementation details.

The `GrayBox` simulator and environment interface includes:
- `GrayBox.Simulation` type to hold simulation variables
- `GrayBox.environment(sim::Simulation)` to return the collection of environment distributions
- `GrayBox.transition!(sim::Simulation)` to transition the simulator, returning the log-likelihood

## Failure and Distance Definition
A *failure* event of the system under test is defined be the user. The user defines the function `BlackBox.isevent` to return an boolean indicating a failure or not given the current state of the simulation. An example failure used in the context of AST would be a collision when stress testing autonomous vehicles or aircraft collision avoidance systems.

The real-valued *distance* metric is used to indicate "how close are we to a failure?" and is defined by the user in the `BlackBox.distance` function. This metric is used to guide the search process towards failures by receiving a signal of the distance to a failure. An example distance metric for the autonomous vehicle problem would be the distance between the autonomous vehicle and a pedestrian, where if a failure is a collision with a pedestrian then we'd like to minimize this distance metric to find failures.

## Citation

If you use this package for research purposes, please cite the following:

[![status](https://joss.theoj.org/papers/04dc39ea89e90938727d789a2e402b0b/status.svg)](https://joss.theoj.org/papers/04dc39ea89e90938727d789a2e402b0b)

```
@article{moss2021pomdpstresstesting,
  title = {{POMDPStressTesting.jl}: Adaptive Stress Testing for Black-Box Systems},
  author = {Robert J. Moss},
  journal = {Journal of Open Source Software},
  year = {2021},
  volume = {6},
  number = {60},
  pages = {2749},
  doi = {10.21105/joss.02749}
}
```
# Installation

To install the package, run:

```julia
] add POMDPStressTesting
```

## Testing
[![Build Status](https://github.com/sisl/POMDPStressTesting.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/sisl/POMDPStressTesting.jl/actions/workflows/CI.yml) [![codecov](https://codecov.io/gh/sisl/POMDPStressTesting.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/sisl/POMDPStressTesting.jl)


To run the test suite, open the Julia Pkg mode using `]` and then run:

```julia
test POMDPStressTesting
```

Testing is automated using Travis CI, which runs the [test/runtests.jl](https://github.com/mossr/POMDPStressTesting.jl/blob/master/test/runtests.jl) file.# [Example Problem](@id example)
[![Example Notebook](https://img.shields.io/badge/example-notebook-blue)](https://nbviewer.jupyter.org/github/sisl/POMDPStressTesting.jl/blob/master/notebooks/Walk1D.ipynb)

This section walks through a full example to illustrate how to use the POMDPStressTesting package. See the [Jupyter notebook](https://nbviewer.jupyter.org/github/sisl/POMDPStressTesting.jl/blob/master/notebooks/Walk1D.ipynb) for the full implementation.

Some definitions to note for this example problem:
- **System**: a one-dimensional walking agent
- **Environment**: distribution of random walking actions, sampled from a standard normal distribution $$\mathcal N(0.1)$$
- **Failure event**: agent walks outside of the ±10 region
- **Distance metric**: how close to the ±10 edge is the agent?

#### Problem Description (Walk1D)
This problem, called Walk1D, samples random walking distances from a standard normal distribution $$\mathcal N(0,1)$$ and defines failures as walking past a certain threshold (which is set to ±10 in this example). AST will either select the seed which deterministically controls the sampled value from the distribution (i.e. from the transition model) or will directly sample the provided environmental distributions. These action modes are determined by the seed-action or sample-action options (i.e. [`ASTSeedAction`](@ref ast_action_type) or [`ASTSampleAction`](@ref ast_action_type), respectively) . AST will guide the simulation to failure events using a notion of distance to failure, while simultaneously trying to find the set of actions that maximizes the log-likelihood of the samples.


> Refer to the [notebook](https://nbviewer.jupyter.org/github/sisl/POMDPStressTesting.jl/blob/master/notebooks/Walk1D.ipynb) for the full implementation.

> For the non-notebook version, see the Julia file [test/Walk1D.jl](https://github.com/mossr/POMDPStressTesting.jl/blob/master/test/Walk1D.jl)# Library/Interface

This section details the interface and functions provided by POMDPStressTesting.jl.

```@meta
CurrentModule = POMDPStressTesting
```

## Contents

```@contents
Pages = ["interface.md"]
```

## Index

```@index
Pages = ["interface.md"]
```

# Modules
```@docs
POMDPStressTesting
```

```@docs
AST
```

```@docs
GrayBox
```

```@docs
BlackBox
```

## [`GrayBox`](@id graybox_interface)
```@docs
GrayBox.Simulation
GrayBox.Sample
GrayBox.Environment
GrayBox.EnvironmentSample
GrayBox.environment
GrayBox.transition!
```

## [`BlackBox`](@id blackbox_interface)
```@docs
BlackBox.initialize!
BlackBox.evaluate!
BlackBox.distance
BlackBox.isevent
BlackBox.isterminal
```

```@meta
CurrentModule = AST
```

## `AST`
```@docs
search!
initialstate
reward
gen
isterminal
discount
random_action
action
actions
convert_s
go_to_state
record
record_trace
get_top_path
rollout
rollout_end
playback
online_path
state_trace
action_trace
q_trace
action_q_trace
reset_metrics!
```

# AST Types
```@docs
ASTParams
ASTAction
ASTSeedAction
ASTSampleAction
ASTState
ASTMetrics
ASTMDP
```

```@meta
CurrentModule = POMDPStressTesting
```

# [Metrics](@id metrics)
```@docs
print_metrics
```

# [Visualizations/Figures](@id visualizations)
```@docs
visualize
episodic_figures
distribution_figures
```# [Guide](@id guide)

To use this package to stress test your own system, the user has to provide the following:
- Implementations of the `BlackBox` system interface functions (outlined [here](@ref blackbox_interface))
- Implementations of the `GrayBox` simulator/environment interface functions (outlined [here](@ref graybox_interface))

## Problem Setup
Once the `GrayBox` and `BlackBox` interfaces are defined (see [Example](@ref example) for a full example), the user has access to the solvers and simply needs to set up the AST problem.

First, set up your simulation structure, where `YourSimulation <: GrayBox.Simulation` (note, you'll need to replace `YourSimulation` with your own structure you've defined as part of the `GrayBox` interface):
```julia
sim::GrayBox.Simulator = YourSimulation()
```

Then, set up the adaptive stress testing (AST) Markov decision process (MDP) structure, given your `sim` object:
```julia
mdp::ASTMDP = ASTMDP{ASTSeedAction}(sim) # could use `ASTSeedAction` or `ASTSampleAction`
```

#### [AST Action Type](@id ast_action_type)
- **`ASTSeedAction`**
    - This action type samples seeds for a random number generator (RNG), which means the `GrayBox.transition!` function must sample from the environment themselves and apply the transition.
        - Useful when it's difficult to access the individual sampled environment outputs.
- **`ASTSampleAction`**
    - This action type samples directory from the `GrayBox.Environment` and will pass the sample(s) (as `GrayBox.EnvironmentSample`) to the `GrayBox.transition!` function to be directly applied
        - Useful when you have full access to the simulation environment and can apply each sample directly in the transition function.

## Solving the AST MDP

Now you can choose your solver (see [Solvers](@ref solvers)) and run `solve` given the AST `mdp` (Markov decision process) to produce an online `planner` (no search has been performed yet):
```julia
solver = MCTSPWSolver()
planner = solve(solver, mdp)
```

## Searching for Failures
Once the problem is set up, you can search for failures using `search!` given your `planner`. This will return the best action trace it found.
```julia
action_trace = search!(planner)
```

## Playback Failures
Given either the `action_trace` or the top `k` performing action traces, using [`get_top_path(mdp)`](@ref), you can playback the particular action:
```julia
final_state = playback(planner, action_trace)
```

## [Metrics and Visualizations](@id metrics_visualizations)
Afterwards, you can look at performance metrics and visualizations (see [Metrics](@ref metrics) and [Visualization](@ref visualizations)):
```julia
print_metrics(planner)
```

To plot, first install the `PyPlot` and `Seaborn` packages and load them. We use [`Requires.jl`](https://github.com/JuliaPackaging/Requires.jl) to handle these dependencies.
```julia
using PyPlot
using Seaborn
```

You can plot episodic metrics, including running miss distance mean, minimum miss distance, and cumulative failures all over episode (i.e. iteration):
```julia
episodic_figures(planner)
```

You can also plot the miss distance distribution and log-likelihood distribution:
```julia
distribution_figures(planner)
```
# Contributing

We welcome all contributions!

- Please [fork the repository](https://github.com/sisl/POMDPStressTesting.jl) and submit a new [Pull Request](https://github.com/sisl/POMDPStressTesting.jl/pulls)
- Report issues through our [GitHub issue tracker](https://github.com/sisl/POMDPStressTesting.jl/issues)
- For further support, either file an [issue](https://github.com/sisl/POMDPStressTesting.jl/issues) or email Robert Moss at [mossr@cs.stanford.edu](mailto:mossr@cs.stanford.edu)

## Style Guide
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

We follow the [Blue](https://github.com/invenia/BlueStyle) style guide for Julia.# [Solvers](@id solvers)
This section describes the solution methods (i.e. solvers) used to in the AST problem formulation to search for likely failure events.

Several solvers are implemented:

## Reinforcement learning
Monte Carlo tree search (MCTS) is used as the search algorithm in the standard AST formulation.
We piggy-back on the [MCTS.jl](https://github.com/JuliaPOMDP/MCTS.jl) package with modification specific to AST.
Modifications to MCTS are described here: [^1]
* [`MCTSPWSolver`](https://github.com/mossr/POMDPStressTesting.jl/blob/master/src/solvers/mcts.jl)

[^1] Robert J. Moss, Ritchie Lee, Nicholas Visser, Joachim Hochwarth, James G. Lopez, Mykel J. Kochenderfer, *"Adaptive Stress Testing of Trajectory Predictions in Flight Management Systems"*, DASC 2020. [https://arxiv.org/abs/2011.02559](https://arxiv.org/abs/2011.02559)

## Deep reinforcement learning
Deep reinforcement learning solvers include *trust region policy optimization* (TRPO) [^2] and *proximal policy optimization* (PPO) [^3].
We'd like to thank [Shreyas Kowshik's](https://github.com/shreyas-kowshik/RL-baselines.jl) for their initial Julia implementation of these methods.
* [`TRPOSolver`](https://github.com/mossr/POMDPStressTesting.jl/blob/master/src/solvers/drl/trpo.jl)
* [`PPOSolver`](https://github.com/mossr/POMDPStressTesting.jl/blob/master/src/solvers/drl/ppo.jl)

[^2] John Schulman, Sergey Levine, Philipp Moritz, Michael I. Jordan, and Pieter Abbeel, *"Trust Region Policy Optimization"*, ICML 2015. [https://arxiv.org/abs/1502.05477](https://arxiv.org/abs/1502.05477)

[^3] John Schulman, Filip Wolski, Prafulla Dhariwal, Alec Radford, and Oleg Klimov, *"Proximal Policy Optimization"*, 2017. [https://arxiv.org/abs/1707.06347](https://arxiv.org/abs/1707.06347)

## Stochastic optimization
Solvers that use stochastic optimization include the [cross-entropy method](https://en.wikipedia.org/wiki/Cross-entropy_method) solver `CEMSolver`. This solution method is adapted from the [CrossEntropyMethod.jl](https://github.com/sisl/CrossEntropyMethod.jl) package by Anthony Corso.
* [`CEMSolver`](https://github.com/mossr/POMDPStressTesting.jl/blob/master/src/solvers/cem.jl)


## Baselines
Baseline solvers are used for comparison to more sophisticated search methods. Currently, the only baseline solver is the `RandomSearchSolver` which uses Monte Carlo rollouts of a random policy to search for failures.
* [`RandomSearchSolver`](https://github.com/mossr/POMDPStressTesting.jl/blob/master/src/solvers/random_search.jl)


