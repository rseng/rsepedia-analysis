# Code Structure in `src/`

## Patterns Module

The `Patterns.jl` file contains type definitions for pattern matching building blocks 
called `Pattern`s, shared between pattern matching backends.
This module provides the type hierarchy required to build patterns, the
left hand side of rules.

## Rules 

The `Rules` folder contains 
- `rules.jl`: definitions for rule types used in various rewriting backends.
- `matchers.jl`: Classical rewriting pattern matcher.

# `Syntax.jl`
Contains the frontend to Rules and Patterns (`@rule` macro and `Pattern` function), using the compatible SymbolicUtils.jl syntax.

# EGraphs Module 
Contains code for the e-graphs rewriting backend. See [egg paper](https://dl.acm.org/doi/pdf/10.1145/3434304) for an high level overview.

- `egraphs.jl`: Definition of `ENode`, `EClass` and `EGraph` types, EClass unioning, metadata access, defintion of EGraphs, adding, merging, rebuilding.
- `ematch.jl`: E-Graph Pattern matching virtual machine interpreter.
- `analysis.jl`: Core algorithms for analyzing egraphs and extracting terms from egraphs.
- `saturation.jl`: Core algorithm for equality saturation, rewriting on e-graphs, e-graphs search.  Search phase of equality saturation. Uses multiple-dispatch on rules, Write phase of equality saturation. Application and instantiation of `Patterns` from matching/search results. Definition of `SaturationParams` type, parameters for equality saturation, Definition of equality saturation execution reports. Utility functions and macros to check equality of terms in egraphs.
- `Schedulers.jl`: Module containing definition of Schedulers for equality saturation. 


## `Library.jl`
Contains utility functions and examples of ready-to-use theories of rules. Macros that generate single rules corresponding to common algebraic properties and macros for generating theories from common algebraic structures.  
<p align="center">
<img width="400px" src="https://raw.githubusercontent.com/juliasymbolics/Metatheory.jl/master/docs/src/assets/dragon.jpg"/>
</p>

# Metatheory.jl

[![Docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://juliasymbolics.github.io/Metatheory.jl/dev/)
[![Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliasymbolics.github.io/Metatheory.jl/stable/)
![CI](https://github.com/juliasymbolics/Metatheory.jl/workflows/CI/badge.svg)
[![codecov](https://codecov.io/gh/juliasymbolics/Metatheory.jl/branch/master/graph/badge.svg?token=EWNYPD7ASX)](https://codecov.io/gh/juliasymbolics/Metatheory.jl)
[![arXiv](https://img.shields.io/badge/arXiv-2102.07888-b31b1b.svg)](https://arxiv.org/abs/2102.07888)
[![status](https://joss.theoj.org/papers/3266e8a08a75b9be2f194126a9c6f0e9/status.svg)](https://joss.theoj.org/papers/3266e8a08a75b9be2f194126a9c6f0e9)
[![Zulip](https://img.shields.io/badge/Chat-Zulip-blue)](https://julialang.zulipchat.com/#narrow/stream/277860-metatheory.2Ejl)

**Metatheory.jl** is a general purpose term rewriting, metaprogramming and algebraic computation library for the Julia programming language, designed to take advantage of the powerful reflection capabilities to bridge the gap between symbolic mathematics, abstract interpretation, equational reasoning, optimization, composable compiler transforms, and advanced
homoiconic pattern matching features. The core features of Metatheory.jl are a powerful rewrite rule definition language, a vast library of functional combinators for classical term rewriting and an *e-graph rewriting*, a fresh approach to term rewriting achieved through an equality saturation algorithm. Metatheory.jl can manipulate any kind of
Julia symbolic expression type, as long as it satisfies the [TermInterface.jl](https://github.com/JuliaSymbolics/TermInterface.jl).

Metatheory.jl provides:
- An eDSL (domain specific language) to define different kinds of symbolic rewrite rules.
- A classical rewriting backend, derived from the [SymbolicUtils.jl](https://github.com/JuliaSymbolics/SymbolicUtils.jl) pattern matcher, supporting associative-commutative rules. It is based on the pattern matcher in the [SICM book](https://mitpress.mit.edu/sites/default/files/titles/content/sicm_edition_2/book.html).
- A flexible library of rewriter combinators.
- An e-graph rewriting (equality saturation) backend and pattern matcher, based on the [egg](https://egraphs-good.github.io/) library, supporting backtracking and non-deterministic term rewriting by using a data structure called *e-graph*, efficiently incorporating the notion of equivalence in order to reduce the amount of user effort required to achieve optimization tasks and equational reasoning.
- `@capture` macro for flexible metaprogramming.

Intuitively, Metatheory.jl transforms Julia expressions
in other Julia expressions and can achieve such at both compile and run time. This allows Metatheory.jl users to perform customized and composable compiler optimizations specifically tailored to single, arbitrary Julia packages.
Our library provides a simple, algebraically composable interface to help scientists in implementing and reasoning about semantics and all kinds of formal systems, by defining concise rewriting rules in pure, syntactically valid Julia on a high level of abstraction. Our implementation of equality saturation on e-graphs is based on the excellent, state-of-the-art technique implemented in the [egg](https://egraphs-good.github.io/) library, reimplemented in pure Julia.

## 1.0 is out!

The first stable version of Metatheory.jl is out! The goal of this release is to unify the symbolic manipulation ecosystem of Julia packages. Many features have been ported from SymbolicUtils.jl. Now, Metatheory.jl can be used in place of SymbolicUtils.jl when you have no need of manipulating mathematical expressions. SymbolicUtils.jl can now completely leverage on the generic stack of rewriting features provided by Metatheory.jl, highly decoupled from the symbolic term representation thanks to [TermInterface.jl](https://github.com/JuliaSymbolics/TermInterface.jl). Read more in [NEWS.md](https://github.com/JuliaSymbolics/Metatheory.jl/blob/master/NEWS.md).

The introduction of  [TermInterface.jl](https://github.com/JuliaSymbolics/TermInterface.jl) has allowed for large potential in generalization of term rewriting and symbolic analysis and manipulation features. It’s been a few months we’ve been talking about the integration between Metatheory.jl with Symbolics.jl, as it has been shown in the ["High-performance symbolic-numerics via multiple dispatch"](https://arxiv.org/abs/2105.03949) paper.

## Recommended Readings - Selected Publications

- The [Metatheory.jl manual](https://juliasymbolics.github.io/Metatheory.jl/stable/) 
- The [Metatheory.jl introductory paper](https://joss.theoj.org/papers/10.21105/joss.03078#) gives a brief high level overview on the library and its functionalities.
- The Julia Manual [metaprogramming section](https://docs.julialang.org/en/v1/manual/metaprogramming/) is fundamental to understand what homoiconic expression manipulation is and how it happens in Julia.
- An [introductory blog post on SIGPLAN](https://blog.sigplan.org/2021/04/06/equality-saturation-with-egg/) about `egg` and e-graphs rewriting.
- [egg: Fast and Extensible Equality Saturation](https://dl.acm.org/doi/pdf/10.1145/3434304) contains the definition of *E-Graphs* on which Metatheory.jl's equality saturation rewriting backend is based. This is a strongly recommended reading.
- [High-performance symbolic-numerics via multiple dispatch](https://arxiv.org/abs/2105.03949): a paper about how we used Metatheory.jl to optimize code generation in [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl)

## Contributing

If you'd like to give us a hand and contribute to this repository you can:
- Find a high level description of the project architecture in [ARCHITECTURE.md](https://github.com/juliasymbolics/Metatheory.jl/blob/master/ARCHITECTURE.md)
- Read the contribution guidelines in [CONTRIBUTING.md](https://github.com/juliasymbolics/Metatheory.jl/blob/master/CONTRIBUTING.md)

## Installation

You can install the stable version:
```julia
julia> using Pkg; Pkg.add("Metatheory")
```

Or you can install the developer version (recommended by now for latest bugfixes)
```julia
julia> using Pkg; Pkg.add(url="https://github.com/JuliaSymbolics/Metatheory.jl")
```

## Documentation

Extensive Metatheory.jl is available [here](https://juliasymbolics.github.io/Metatheory.jl/dev)

## Citing

If you use Metatheory.jl in your research, please [cite](https://github.com/juliasymbolics/Metatheory.jl/blob/master/CITATION.bib) our works.

--- 

# Sponsors

If you enjoyed Metatheory.jl and would like to help, you can donate a coffee or choose place your logo and name in this page. [See 0x0f0f0f's Github Sponsors page](https://github.com/sponsors/0x0f0f0f/)!

<p align="center">
<a href="https://planting.space"> 
    <img width="300px" src="https://raw.githubusercontent.com/juliasymbolics/Metatheory.jl/master/.github/plantingspace.png"/>
</a>
</p>
## 1.2
- Fixes when printing patterns
- Can pass custom `similarterm` to `SaturationParams` by using `SaturationParams.simterm`.

## 1.1
- EGraph pattern matcher can now match against both symbols and function objects
- Fixes for Symbolics.jl integration


## 1.0

Metatheory.jl + SymbolicUtils.jl = ❤️

- Metatheory.jl now supports the same syntax as [SymbolicUtils.jl](https://github.com/JuliaSymbolics/SymbolicUtils.jl/) for the rule definition DSL!
- The classical pattern matcher has been redesigned, and it is a port of [SymbolicUtils.jl](https://github.com/JuliaSymbolics/SymbolicUtils.jl/)'s pattern matcher. Now Metatheory.jl can be used in place of SU's rewriting backend.
- Performance improvements: caching of ground terms when doing e-matching in equality saturation.
- Dynamic Rules do not use RuntimeGeneratedFunctions when not needed.
- Removed `@metatheory_init`
- Rules now support type and function predicates as in SymbolicUtils.jl
- Redesigned the library
- Introduced `@timerewrite` to time the execution of classical rewriting systems.# Notes for Metatheory.jl Contributors

Welcome, and thanks for considering Metatheory.jl! Please be sure to respect our [community standards](https://julialang.org/community/standards) in all interactions.

We gratefully acknowledge the general [Julia CONTRIBUTING.md document](https://github.com/JuliaLang/julia/blob/master/CONTRIBUTING.md), from which much of this was adapted.

## Learning Julia

A pre-requisite for using Metatheory.jl is to know at least a little about Julia itself. [The learning page](https://julialang.org/learning) has a great list of resources for new and experienced users alike. [This tutorial video](https://www.youtube.com/watch?v=vWkgEddb4-A) is one recommended starting point, as is the "[Invitation to Julia](https://www.youtube.com/watch?v=gQ1y5NUD_RI)" workshop video from JuliaCon 2015  ([slide materials here](https://github.com/dpsanders/invitation_to_julia)). The [Julia documentation](https://docs.julialang.org) covers the language and core library features, and is searchable.

## Learning Metatheory.jl

Our [main documentaion](https://github.com/JuliaSymbolics/Metatheory.jl/) provides an overview and some examples of using Metatheory.jl.
The core package is hosted at [Metatheory.jl](https://github.com/JuliaSymbolics/Metatheory.jl/).

## Before filing an issue

Julia's own "[How to file a bug report](https://github.com/JuliaLang/julia/blob/master/CONTRIBUTING.md#how-to-file-a-bug-report)" has many useful tips to help make sure that all necessary information is included.

Try to report the issue in the package responsible for the error.
You can often make good guesses by examining the backtrace (in cases where an
error is thrown), using `@which`, stepping in with the debugger, or just
using the search bar at the top left of [Metatheory.jl](https://github.com/JuliaSymbolics/Metatheory.jl/).

## Contributing documentation

*By contributing you agree to be bound by Metatheory.jl' MIT license*

Many documentation issues are easy! For small changes, you can just click on one of the files in the `docs/src` directory, click on the "pencil icon," and [edit it in your browser](https://help.github.com/en/github/managing-files-in-a-repository/editing-files-in-another-users-repository). Any changes you suggest will first be vetted by an experienced developer, so there is no need to worry that you'll mess something up.

Changes to the "docstrings" (the string preceding a method in source code) should be made in the package in which they appear.

For bigger documentation changes, it is probably best to clone the package and submit the changes as an ordinary pull request, as described below under "Contributing code." You can build the package locally if you install [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl), and run `include("make.jl")` in the `docs/` folder. To see the completed documentation, open the `build/index.md` file in your browser.

## Contributing code

*By contributing you agree to be bound by Metatheory.jl' MIT license*

If you've never submitted a pull request before, it can take a little while to become familiar with the process. In addition to the steps below, [GitHub has a tutorial and exercises](https://try.github.io/). See also the excellent [Git book](https://git-scm.com/book/en/v2). There are also many good external tutorials on this subject, like [this one](https://yangsu.github.io/pull-request-tutorial/).

### Contributor Checklist

* Create a [GitHub account](https://github.com/signup/free).

* If you plan to fix a bug, feel free to first report the bug as an issue on its own.
  In the text, you can mention whether you're planning on addressing it yourself.
  *Pro tip*: if you do submit a pull request to fix it, put "Fixes #<issue number>" in the commit message and it will close automatically when your pull request is merged.

  If you're concerned your change might be controversial, you can also use an issue to propose your change in general terms and discuss it before implementation.

* Fork whatever repository you plan to commit to by clicking on the "Fork" button at the upper-right of the home page.

* If you haven't already implemented your changes, check the package out for development: hit `]` in the Julia REPL and then type (for example) `dev Images`.
You'll get a copy of the full repository in your `~/.julia/dev` folder. See the [package manager documentation](https://julialang.github.io/Pkg.jl/v1/) for further details.

* Make your changes. Generally you should be working on a branch, so your work doesn't conflict with ongoing development in the `master` branch. Ensure you follow the [Julia style guide](https://docs.julialang.org/en/v1/manual/style-guide/index.html) for your contribution.

* Test your changes. We aspire to have test coverage for every bit of "user visible" functionality. Tests are stored, appropriately, in the `test/` folder of each package. You can run existing tests yourself and add new ones. Sometimes testing is more work than the actual change itself, but having tests ensures that no well-meaning future developer will accidentally mess up your functionality---it's worth it!  *"A fix is for today. A test is forever."*

* Submit your changes up to your fork and then submit a pull request-!

* See what happens to the automated tests that run on Github Actions. If there are errors, check the logs and see whether they look like they are related to your changes; if so, try to fix the problem by adding new commits to your pull request. Once the tests pass, hooray! :tada:

* Relax and wait for feedback. We try to review contributions quickly and courteously. But we are human, and sometimes we get busy with other things or fail to notice an email; if it's been a while since you submitted your pull request, try posting a polite reminder about the existence of your pull request.

* Discuss any feedback you receive as necessary. It's fine to defend your approach, but also be open to making changes based on suggestions you receive.

* Sooner or later, the fate of your pull request will become clear. If it gets approved, an established contributor will merge it. It's not officially released into the wild until a contributor releases a new version of the package; if that doesn't happen quickly, don't hesitate to make an inquiry in case it's simply been overlooked.

From the whole team, thanks in advance for your contribution!

### Contribution tips

* [Revise](https://github.com/timholy/Revise.jl) is a package that
tracks changes in source files and automatically updates function
definitions in your running Julia session. Using it, you can make
extensive changes without needing to rebuild the package in order to test
your changes.

* Debuggers can help you get to the root of a problem. There are many choices and interfaces:
  + [Juno](https://github.com/JunoLab/Juno.jl) has a polished GUI for debugging
  + [Debugger](https://github.com/JuliaDebug/Debugger.jl) has a polished command-line interface
  + [Rebugger](https://github.com/timholy/Rebugger.jl) has an innovative but somewhat less-polished command-line interface
  + [Infiltrator](https://github.com/JuliaDebug/Infiltrator.jl) offers more limited debugging, but often it's precisely what you need while avoiding the performance penalties that some of the other options suffer from.
---
title: 'Metatheory.jl: Fast and Elegant Algebraic Computation in Julia with Extensible Equality Saturation'
tags:
  - Julia
  - compiler
  - symbolic
  - algebra
  - rewriting
  - optimization
authors:
  - name: Alessandro Cheli #^[Custom footnotes for e.g. denoting who the corresponding author is can be included like this.]
    orcid: 0000-0002-8122-9469
    affiliation: 1 # (Multiple affiliations must be quoted)
affiliations:
 - name: University of Pisa, Pisa, Italy
   index: 1
date: 11 February 2021
bibliography: paper.bib

---

# Statement of Need

<!-- The Julia programming language is a fresh approach to technical computing [@bezanson2017julia], disrupting the popular conviction that a programming language cannot be very high level, easy to learn, and performant at the same time. One of the most practical features of Julia is the excellent metaprogramming and macro system, allowing for programmatic generation and manipulation of Julia expressions as first-class values in the core language, with a well-known paradigm similar to LISP idioms such as Scheme,
a programming language property colloquially referred to as *homoiconicity*. -->

The Julia programming language is a fresh approach to technical computing [@bezanson2017julia], disrupting the popular conviction that a programming language cannot be high-level, easy to learn, and performant at the same time. One of the most practical features of Julia is the excellent metaprogramming and macro system, allowing for *homoiconicity*: programmatic generation and manipulation of expressions as first-class values, a well-known paradigm found in LISP dialects such as Scheme.

Metatheory.jl is a general-purpose metaprogramming and algebraic computation library for the Julia programming language, designed to take advantage of its powerful reflection capabilities to bridge the gap between symbolic mathematics,
abstract interpretation, equational reasoning, optimization, composable compiler transforms, and advanced homoiconic pattern-matching features. Intuitively, Metatheory.jl transforms Julia expressions into other Julia expressions at both compile time and run time. This allows users to perform customized and composable compiler optimizations that are specifically tailored to single, arbitrary Julia packages. The library provides a simple, algebraically composable interface to help scientists to implement and reason about all kinds of formal systems, by defining concise rewriting rules as syntactically-valid Julia code. The primary benefit of using Metatheory.jl is the algebraic nature of the specification of the rewriting system. Composable blocks of rewrite rules bear a strong resemblance to algebraic
structures encountered in everyday scientific literature.

<!-- Rewrite rules are defined as regular Julia expressions, manipulating other syntactically valid Julia expressions: since Julia supports LaTeX-like abbreviations of UTF8 mathematical symbols as valid operators and symbols,
rewrite theories in Metatheory.jl can bear a strong structural and visual resemblance to mathematical formalisms encountered in paper literature. -->

<!-- Theories can then be executed through two, highly composable, rewriting backends. The first backend relies on a *classic* fixed-point recursive iteration of AST, with a match-and-replace algorithm built on top of the [@matchcore] pattern matcher. This backend is suitable for deterministic recursive algorithms that intensively use pattern matching on syntax trees, for example, defining an interpreter from operational or denotational semantics. Nevertheless, when using this classical approach, even trivial equational rules such as commutativity and associativity may cause the rewriting algorithm to loop indefinitely, or to return unexpected results. This is known as *rewrite order* and is notoriously recognized for requiring extensive user reasoning about the ordering and structuring of rules to ensure termination. -->

# Summary

Metatheory.jl offers a concise macro system to define *theories*: composable blocks of rewriting rules that can be executed through two, highly composable, rewriting backends. The first is based on standard rewriting, built on top of the pattern matcher developed in @matchcore.
This approach, however, suffers from the usual problems of rewriting systems. For example, even trivial equational rules such as commutativity may lead to non-terminating systems and thus need to be adjusted by some sort of structuring or rewriting order, which is known to require extensive user reasoning.


The other back-end for Metatheory.jl, the core of our contribution, is designed so that it does not require the user to reason about rewriting order. To do so it relies on equality saturation on *e-graphs*, the state-of-the-art technique adapted from the `egg` Rust library [@egg].

*E-graphs* can compactly represent many equivalent expressions and programs. Provided with a theory of rewriting rules, defined in pure Julia, the *equality saturation* process iteratively executes an e-graph-specific pattern matcher and inserts the matched substitutions. Since e-graphs can contain loops, infinite derivations can be represented compactly and it is not required that the described rewrite system be terminating or confluent.

The saturation process relies on the definition of e-graphs to include *rebuilding*, i.e. the automatic process of propagation and maintenance of congruence closures.
One of the core contributions of @egg is a delayed e-graph rebuilding process that is executed at the end of each saturation step, whereas previous definitions of e-graphs in the literature included rebuilding after every rewrite operation.
Provided with *equality saturation*, users can efficiently derive (and analyze) all possible equivalent expressions contained in an e-graph. The saturation process can be required to stop prematurely as soon as chosen properties about the e-graph and its expressions are proved. This latter back-end based on *e-graphs* is suitable for partial evaluators, symbolic mathematics, static analysis, theorem proving and superoptimizers.

<!-- The other back-end for Metatheory.jl, the core of our contribution, is designed to not require the user to reason about rewriting order by employing equality saturation on e-graphs. This backend allows programmers to define equational theories in pure Julia without worrying about rule ordering and structuring, by relying on state-of-the-art techniques for equality saturation over *e-graphs* adapted from the `egg` Rust library [@egg].
Provided with a theory of equational rewriting rules, *e-graphs* compactly represent many equivalent programs. Saturation iteratively executes an e-graph specific pattern matcher to efficiently compute (and analyze) all possible equivalent expressions contained in the e-graph congruence closure. This latter back-end is suitable for partial evaluators, symbolic mathematics, static analysis, theorem proving and superoptimizers. -->

![These four e-graphs represent the process of equality saturation, adding many equivalent ways to write $a * (2 * 3) / 6$ after each iteration. \label{fig:egraph}](egraphs.png)


The original `egg` library [@egg] is
the first implementation of generic and extensible e-graphs [@nelson1980fast]; the contributions of `egg` include novel amortized algorithms for fast and efficient equivalence saturation and analysis.
Differently from the original Rust implementation of `egg`, which handles expressions defined as Rust strings and data structures, our system directly manipulates homoiconic Julia expressions, and can therefore fully leverage the Julia subtyping mechanism [@zappa2018julia], allowing programmers to build expressions containing not only symbols but all kinds of Julia values.
This permits rewriting and analyses to be efficiently based on runtime data contained in expressions. Most importantly, users can -- and are encouraged to -- include type assertions in the left-hand side of rewriting rules in theories.

One of the project goals of Metatheory.jl, beyond being easy to use and composable, is to be fast and efficient. Both the first-class pattern matching system and the generation of e-graph analyses from theories rely on RuntimeGeneratedFunctions.jl [@rgf], generating callable functions at runtime that efficiently bypass Julia's world age problem (explained and formalized in @belyakova2020world) with the full performance of a standard Julia anonymous function.


## Analyses and Extraction

With Metatheory.jl, modeling analyses and conditional/dynamic rewrites are straightforward. It is possible to check conditions on runtime values or to read and write from external data structures during rewriting. The analysis mechanism described in `egg` [@egg] and re-implemented in our contribution lets users define ways to compute additional analysis metadata from an arbitrary semi-lattice domain, such as costs of nodes or logical statements attached to terms. Other than for inspection, analysis data can be used to modify expressions in the e-graph both during rewriting steps and after e-graph saturation.

Therefore using the equality saturation (e-graph) backend, extraction can be performed as an on-the-fly e-graph analysis or after saturation. Users
can define their own cost function, or choose between a variety of predefined cost functions for automatically extracting the best-fitting expressions from an equivalence class represented in an e-graph.

# Example Usage

In this example, we build rewrite systems, called `theories` in Metatheory.jl, for simplifying expressions
in the usual commutative monoid of multiplication and the commutative group of addition, and we compose
the `theories` together with a *constant folding* theory. The pattern matcher for the e-graphs backend
allows us to use the existing Julia type hierarchy for integers and floating-point numbers with a high level
of abstraction. As a contribution over the original egg [@egg] implementation, left-hand sides of rules in Metatheory.jl can contain type assertions on pattern variables, to give rules that depend on consistent type hierarchies and  to seamlessly access literal Julia values in the right-hand side of dynamic rules.

We finally introduce two simple rules for simplifying fractions, that
for the sake of simplicity, do not check any additional analysis data.
\autoref{fig:egraph} contains a friendly visualization of a consistent fragment of the equality saturation process in this example.
You can see how loops evidently appear in the definition of the rewriting rules.
While the classic rewriting backend would loop indefinitely or stop early when repeatedly matching these rules,
the e-graph backend natively supports this level of abstraction and allows the
programmer to completely forget about the ordering and looping of rules.
Efficient scheduling heuristics are applied automatically to prevent instantaneous
combinatorial explosion of the e-graph, thus preventing substantial slowdown of the equality saturation
process.

```julia
using Metatheory
using Metatheory.EGraphs

comm_monoid = @theory begin
  # commutativity
  a * b => b * a
  # identity
  a * 1 => a
  # associativity
  a * (b * c) => (a * b) * c
  (a * b) * c => a * (b * c)
end;

comm_group = @theory begin
  # commutativity
  a + b => b + a
  # identity
  a + 0 => a
  # associativity
  a + (b + c) => (a + b) + c
  (a + b) + c => a + (b + c)
  # inverse
  a + (-a) => 0
end;

# dynamic rules are defined with the `|>` operator
folder = @theory begin
  a::Real + b::Real |> a+b
  a::Real * b::Real |> a*b
end;

div_sim = @theory begin
  (a * b) / c => a * (b / c)
  a::Real / a::Real  |>  (a != 0 ? 1 : error("division by 0"))
end;

t = union(comm_monoid, comm_group, folder, div_sim) ;

g = EGraph(:(a * (2*3) / 6)) ;
saturate!(g, t) ;
ex = extract!(g, astsize)
# :a

```

# Conclusion

Many applications of equality saturation to advanced optimization tasks have been recently published. Herbie [@panchekha2015automatically]
is a tool for automatically improving the precision of floating point expressions, which recently switched to `egg` as the core rewriting backend. However, Herbie requires interoperation and conversion of expressions between different languages and libraries. In @yang2021equality, the authors used `egg` to superoptimize tensor signal flow graphs describing neural networks.  Implementing similar case studies in pure Julia would make valid research contributions on their own. We are confident that a well-integrated and homoiconic equality saturation engine in pure Julia will permit exploration of many new metaprogramming applications, and allow them to be implemented in an elegant, performant and concise way. Code for Metatheory.jl is available in @metatheory, or at [https://github.com/0x0f0f0f/Metatheory.jl](https://github.com/0x0f0f0f/Metatheory.jl).

# Acknowledgements

We acknowledge Max Willsey and contributors for their work on the original `egg` library [@egg], Christopher Rackauckas and Christopher Foster for their efforts in developing RuntimeGeneratedFunctions [@rgf], Taine Zhao for developing MLStyle [@mlstyle] and MatchCore [@matchcore], and Philip Zucker for his original idea of implementing E-Graphs in Julia [@philzuck1; @philzuck2] and support during the development of the project. Special thanks to Filippo Bonchi for a friendly review of a preliminary version of this article.

# References
This is a simple script to convert Metatheory.jl <https://github.com/0x0f0f0f/Metatheory.jl> theories into an Egg <https://egraphs-good.github.io/> query for comparison.

Get a rust toolchain <https://rustup.rs/>

Make a new project 

```
cargo new my_project
cd my_project
```

Add egg as a dependency to the Cargo.toml. Add the last line shown here.

```
[package]
name = "autoegg"
version = "0.1.0"
authors = ["Philip Zucker <philzook58@gmail.com>"]
edition = "2018"

# See more keys and their definitions at https://doc.rust-lang.org/cargo/reference/manifest.html

[dependencies]
egg = "0.6.0"
```

Copy and paste the Julia script in the project folder. Replace the example theory and query with yours in the script

Run it

```
julia gen_egg.jl
```

Now you can run it in Egg

```
cargo run --release
```

Profit.
# Classical Term Rewriting

## Rule-based rewriting

Rewrite rules match and transform an expression. A rule is written using either
the `@rule` or `@theory` macros. It creates a callable `Rule` object.

### Basics of rule-based term rewriting in Metatheory.jl

**NOTE:** for a real world use case using mathematical constructs, please refer
to [SymbolicUtils.jl](https://github.com/JuliaSymbolics/SymbolicUtils.jl). SU
provides optimized types for mathematical expressions, code generation and a
polished set of rules for simplification.

Here is a simple symbolic rewrite rule, that uses formula for the double angle of the sine function:

```julia:rewrite1
using Metatheory

r1 = @rule sin(2(~x)) --> 2sin(~x)*cos(~x)

expr = :(sin(2z))
r1(expr)
```

The `@rule` macro takes a pair of patterns  -- the _matcher_ and the _consequent_ (`@rule matcher OPERATOR consequent`). If an expression matches the matcher pattern, it is rewritten to the consequent pattern. `@rule` returns a callable object that applies the rule to an expression. There are different kinds of rule in Metatheory.jl:

**Rule operators**:
- `LHS => RHS`: create a `DynamicRule`. The RHS is *evaluated* on rewrite.
- `LHS --> RHS`: create a `RewriteRule`. The RHS is **not** evaluated but *symbolically substituted* on rewrite.
- `LHS == RHS`: create a `EqualityRule`. In e-graph rewriting, this rule behaves like `RewriteRule` but can go in both directions. Doesn't work in classical rewriting.
- `LHS ≠ RHS`: create a `UnequalRule`. Can only be used in e-graphs, and is used to eagerly stop the process of rewriting if LHS is found to be equal to RHS.


You can use **dynamic rules**, defined with the `=>`
operator, to dynamically compute values in the right hand of expressions. This is the default behaviour of rules in [SymbolicUtils.jl](https://github.com/JuliaSymbolics/SymbolicUtils.jl)
Dynamic rules, are similar to anonymous functions. Instead of a symbolic
substitution, the right hand of a dynamic `=>` rule is evaluated during
rewriting: the values that produced a match are bound to the pattern variables.

`~x` in the example is what is a **slot variable** (or *pattern* variable) named `x`. In a matcher pattern, slot variables are placeholders that match exactly one expression. When used on the consequent side, they stand in for the matched expression. If a slot variable appears twice in a matcher pattern, **in classical rewriting** all corresponding matches must be equal (as tested by `Base.isequal` function). Hence this rule says: if you see something added to itself, make it twice of that thing, and works as such.

If you try to apply this rule to an expression with triple angle, it will return `nothing` -- this is the way a rule signifies failure to match.
```julia:rewrite2
r1(:(sin(3z))) === nothing
```

Slot variable (matcher) is not necessary a single variable

```julia:rewrite3
r1(:(sin(2*(w-z))))
```

but it must be a single expression

```julia:rewrite4
r1(:(sin(2*(w+z)*(α+β)))) === nothing
```

Rules are of course not limited to single slot variable

```julia:rewrite5
r2 = @rule sin(~x + ~y) --> sin(~x)*cos(~y) + cos(~x)*sin(~y);

r2(:(sin(α+β)))
```

If you want to match a variable number of subexpressions at once, you will need a **segment variable**. `~xs...` in the following example is a segment variable:

```julia:rewrite6
@rule(+(~xs...) => xs)(:(x + y + z))
```

`~xs` is a vector of subexpressions matched. You can use it to construct something more useful:

```julia:rewrite7
r3 = @rule *(~ys...)^~x => :((*)($(map(y-> :($y^$x), ys)...)));

r3(:((w*w*α*β)^2))
```

### Predicates for matching

Matcher pattern may contain slot variables with attached predicates, written as `~x::p` where `p` is either

- A function that takes a matched expression and returns a boolean value. Such a slot will be considered a match only if `p` returns true.
- A Julia type. Will be considered a match if and only if the value matching against `x` has a type that is a subtype of `p` (`typeof(x) <: p`)

Similarly `~x::g...` is a way of attaching a predicate `g` to a segment variable. In the case of segment variables `g` gets a vector of 0 or more expressions and must return a boolean value. If the same slot or segment variable appears twice in the matcher pattern, then at most one of the occurance should have a predicate.

For example,

```julia:pred1
r = @rule +(~x, ~y::(ys->iseven(length(ys)))...) => "odd terms";

@show r(:(a + b + c + d))
@show r(:(b + c + d))
@show r(:(b + c + b))
@show r(:(a + b))
```


### Declaring Slots

Slot variables can be declared without the `~` using the `@slots` macro

```julia:slots1
@slots x y @rule sin(x + y) => sin(x)*cos(y) + cos(x)*sin(y);
```

This works for segments as well:

```julia:slots2
@slots xs @rule(+(~xs...) => xs);
```

The `@slots` macro is superfluous for the `@rule`, `@capture` and `@theory` macros.
Slot variables may be declared directly as the first arguments to those macros:

```julia:slots3
@rule x y sin(x + y) => sin(x)*cos(y) + cos(x)*sin(y);
```

### Theories

In almost all use cases, it is practical to define many rules grouped together.
A set of rewrite rules and equalities is called a *theory*, and can be defined with the
`@theory` macro. This macro is just syntax sugar to define vectors of rules in a nice and readable way. 


```julia
t = @theory x y z begin 
    x * (y + z) --> (x * y) + (x * z)
    x + y       ==  (y + x)
    #...
end;
```

Is the same thing as writing

```julia
v = [
    @rule x y z  x * (y + z) --> (x * y) + (x * z)
    @rule x y x + y == (y + x)
    #...
];
```

Theories are just collections and
can be composed as regular Julia collections. The most
useful way of composing theories is unioning
them with the '∪' operator.
You are not limited to composing theories, you can
manipulate and create them at both runtime and compile time
as regular vectors.

```julia
using Metatheory
using Metatheory.Library

comm_monoid = @commutative_monoid (*) 1
comm_group = @theory a b c begin
    a + 0 --> a
    a + b --> b + a
    a + inv(a) --> 0 # inverse
    a + (b + c) --> (a + b) + c
end
distrib = @theory a b c begin
    a * (b + c) => (a * b) + (a * c)
end
t = comm_monoid ∪ comm_group ∪ distrib
```

## Composing rewriters

Rules may be *chained together* into more
sophisticated rewirters to avoid manual application of the rules. A rewriter is
any callable object which takes an expression and returns an expression or
`nothing`. If `nothing` is returned that means there was no changes applicable
to the input expression. The Rules we created above are rewriters.

The `Metatheory.Rewriters` module contains some types which create and transform
rewriters.

- `Empty()` is a rewriter which always returns `nothing`
- `Chain(itr)` chain an iterator of rewriters into a single rewriter which applies
   each chained rewriter in the given order.
   If a rewriter returns `nothing` this is treated as a no-change.
- `RestartedChain(itr)` like `Chain(itr)` but restarts from the first rewriter once on the
   first successful application of one of the chained rewriters.
- `IfElse(cond, rw1, rw2)` runs the `cond` function on the input, applies `rw1` if cond
   returns true, `rw2` if it retuns false
- `If(cond, rw)` is the same as `IfElse(cond, rw, Empty())`
- `Prewalk(rw; threaded=false, thread_cutoff=100)` returns a rewriter which does a pre-order 
   (*from top to bottom and from left to right*) traversal of a given expression and applies 
   the rewriter `rw`. `threaded=true` will use multi threading for traversal.
   Note that if `rw` returns `nothing` when a match is not found, then `Prewalk(rw)` will
   also return nothing unless a match is found at every level of the walk. If you are
   applying multiple rules, then `Chain` already has the appropriate passthrough behavior.
   If you only want to apply one rule, then consider using `PassThrough`.
   `thread_cutoff` 
   is the minimum number of nodes in a subtree which should be walked in a threaded spawn.
- `Postwalk(rw; threaded=false, thread_cutoff=100)` similarly does post-order 
   (*from left to right and from bottom to top*) traversal.
- `Fixpoint(rw)` returns a rewriter which applies `rw` repeatedly until there are no changes to be made.
- `FixpointNoCycle` behaves like `Fixpoint` but instead it applies `rw` repeatedly only while it is returning new results.
- `PassThrough(rw)` returns a rewriter which if `rw(x)` returns `nothing` will instead
   return `x` otherwise will return `rw(x)`.

### Chaining rewriters

Several rules may be chained to give chain of rules. Chain is an array of rules which are subsequently applied to the expression.
Important feature of `Chain` is that it returns the expression instead of `nothing` if it doesn't change the expression
It is important to notice, that chain is ordered, so if rules are in different order it wouldn't work the same as in earlier example


One way to circumvent the problem of order of applying rules in chain is to use
`RestartedChain`, it restarts the chain after each successful application of a
rule, so after a rule is hit it (re)starts again and it can apply all the other
rules to the resulting expression. You can also use `Fixpoint` to apply the
rules until there are no changes.

# Metatheory.jl 1.0

```@raw html
<p align="center">
<img width="400px" src="https://raw.githubusercontent.com/juliasymbolics/Metatheory.jl/master/docs/src/assets/dragon.jpg"/>
</p>
```

[![Docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://juliasymbolics.github.io/Metatheory.jl/dev/)
[![Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliasymbolics.github.io/Metatheory.jl/stable/)
![CI](https://github.com/juliasymbolics/Metatheory.jl/workflows/CI/badge.svg)
[![codecov](https://codecov.io/gh/juliasymbolics/Metatheory.jl/branch/master/graph/badge.svg?token=EWNYPD7ASX)](https://codecov.io/gh/juliasymbolics/Metatheory.jl)
[![arXiv](https://img.shields.io/badge/arXiv-2102.07888-b31b1b.svg)](https://arxiv.org/abs/2102.07888)
[![status](https://joss.theoj.org/papers/3266e8a08a75b9be2f194126a9c6f0e9/status.svg)](https://joss.theoj.org/papers/3266e8a08a75b9be2f194126a9c6f0e9)
[![Zulip](https://img.shields.io/badge/Chat-Zulip-blue)](https://julialang.zulipchat.com/#narrow/stream/277860-metatheory.2Ejl)

**Metatheory.jl** is a general purpose term rewriting, metaprogramming and algebraic computation library for the Julia programming language, designed to take advantage of the powerful reflection capabilities to bridge the gap between symbolic mathematics, abstract interpretation, equational reasoning, optimization, composable compiler transforms, and advanced
homoiconic pattern matching features. The core features of Metatheory.jl are a powerful rewrite rule definition language, a vast library of functional combinators for classical term rewriting and an *e-graph rewriting*, a fresh approach to term rewriting achieved through an equality saturation algorithm. Metatheory.jl can manipulate any kind of
Julia symbolic expression type, as long as it satisfies the [TermInterface.jl](https://github.com/JuliaSymbolics/TermInterface.jl).

Metatheory.jl provides:
- An eDSL (domain specific language) to define different kinds of symbolic rewrite rules.
- A classical rewriting backend, derived from the [SymbolicUtils.jl](https://github.com/JuliaSymbolics/SymbolicUtils.jl) pattern matcher, supporting associative-commutative rules. It is based on the pattern matcher in the [SICM book](https://mitpress.mit.edu/sites/default/files/titles/content/sicm_edition_2/book.html).
- A flexible library of rewriter combinators.
- An e-graph rewriting (equality saturation) backend and pattern matcher, based on the [egg](https://egraphs-good.github.io/) library, supporting backtracking and non-deterministic term rewriting by using a data structure called *e-graph*, efficiently incorporating the notion of equivalence in order to reduce the amount of user effort required to achieve optimization tasks and equational reasoning.
- `@capture` macro for flexible metaprogramming.

Intuitively, Metatheory.jl transforms Julia expressions
in other Julia expressions and can achieve such at both compile and run time. This allows Metatheory.jl users to perform customized and composable compiler optimizations specifically tailored to single, arbitrary Julia packages.
Our library provides a simple, algebraically composable interface to help scientists in implementing and reasoning about semantics and all kinds of formal systems, by defining concise rewriting rules in pure, syntactically valid Julia on a high level of abstraction. Our implementation of equality saturation on e-graphs is based on the excellent, state-of-the-art technique implemented in the [egg](https://egraphs-good.github.io/) library, reimplemented in pure Julia.

## 1.0 is out!

The first stable version of Metatheory.jl is out! The goal of this release is to unify the symbolic manipulation ecosystem of Julia packages. Many features have been ported from SymbolicUtils.jl. Now, Metatheory.jl can be used in place of SymbolicUtils.jl when you have no need of manipulating mathematical expressions. SymbolicUtils.jl can now completely leverage on the generic stack of rewriting features provided by Metatheory.jl, highly decoupled from the symbolic term representation thanks to [TermInterface.jl](https://github.com/JuliaSymbolics/TermInterface.jl). Read more in [NEWS.md](https://github.com/JuliaSymbolics/Metatheory.jl/blob/master/NEWS.md).

## Recommended Readings - Selected Publications

- The [Metatheory.jl manual](https://juliasymbolics.github.io/Metatheory.jl/stable/) 
- The [Metatheory.jl introductory paper](https://joss.theoj.org/papers/10.21105/joss.03078#) gives a brief high level overview on the library and its functionalities.
- The Julia Manual [metaprogramming section](https://docs.julialang.org/en/v1/manual/metaprogramming/) is fundamental to understand what homoiconic expression manipulation is and how it happens in Julia.
- An [introductory blog post on SIGPLAN](https://blog.sigplan.org/2021/04/06/equality-saturation-with-egg/) about `egg` and e-graphs rewriting.
- [egg: Fast and Extensible Equality Saturation](https://dl.acm.org/doi/pdf/10.1145/3434304) contains the definition of *E-Graphs* on which Metatheory.jl's equality saturation rewriting backend is based. This is a strongly recommended reading.
- [High-performance symbolic-numerics via multiple dispatch](https://arxiv.org/abs/2105.03949): a paper about how we used Metatheory.jl to optimize code generation in [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl)

## Contributing

If you'd like to give us a hand and contribute to this repository you can:
- Find a high level description of the project architecture in [ARCHITECTURE.md](https://github.com/juliasymbolics/Metatheory.jl/blob/master/ARCHITECTURE.md)
- Read the contribution guidelines in [CONTRIBUTING.md](https://github.com/juliasymbolics/Metatheory.jl/blob/master/CONTRIBUTING.md)

If you enjoyed Metatheory.jl and would like to help, please also consider a [tiny donation 💕](https://github.com/sponsors/0x0f0f0f/)!

## Installation

You can install the stable version:
```julia
julia> using Pkg; Pkg.add("Metatheory")
```

Or you can install the developer version (recommended by now for latest bugfixes)
```julia
julia> using Pkg; Pkg.add(url="https://github.com/JuliaSymbolics/Metatheory.jl")
```

## Documentation

Extensive Metatheory.jl is available [here](https://juliasymbolics.github.io/Metatheory.jl/dev)

## Citing

If you use Metatheory.jl in your research, please [cite](https://github.com/juliasymbolics/Metatheory.jl/blob/master/CITATION.bib) our works.

--- 

```@raw html
<p align="center">
<a href="https://planting.space"> 
    <img width="300px" src="https://raw.githubusercontent.com/juliasymbolics/Metatheory.jl/master/.github/plantingspace.png"/>
</a>
</p>
```# API Documentation


## Syntax

```@autodocs
Modules = [Metatheory.Syntax]
```

---

## Patterns

```@autodocs
Modules = [Metatheory.Patterns]
```

---

## Rules 

```@autodocs
Modules = [Metatheory.Rules]
```

---

## Rules 

```@autodocs
Modules = [Metatheory.Rules]
```

---

## Rewriters

```@autodocs
Modules = [Metatheory.Rewriters]
```

---

## EGraphs

```@autodocs
Modules = [Metatheory.EGraphs]
```

---

## EGraph Schedulers

```@autodocs
Modules = [Metatheory.EGraphs.Schedulers]
```# Interfacing with Metatheory.jl

This section is for Julia package developers who may want to use the rule rewriting systems on their own expression types.

## Defining the interface

Metatheory.jl matchers can match any Julia object that implements an interface to traverse it as a tree. The interface in question, is defined in the [TermInterface.jl](https://github.com/JuliaSymbolics/TermInterface.jl) package. Its purpose is to provide a shared interface between various symbolic programming Julia packages. 

In particular, you should define methods from TermInterface.jl for an expression tree type `T` with symbol types `S` to  work
with SymbolicUtils.jl

You can read the documentation of [TermInterface.jl](https://github.com/JuliaSymbolics/TermInterface.jl) on the [Github repository](https://github.com/JuliaSymbolics/TermInterface.jl).
# EGraphs and Equality Saturation

An *EGraph* is an efficient data structure for representing congruence relations.
EGraphs are data structures originating from theorem provers. Several projects
have very recently repurposed EGraphs to implement state-of-the-art,
rewrite-driven compiler optimizations and program synthesizers using a technique
known as equality saturation. Metatheory.jl provides a general purpose,
customizable implementation of EGraphs and equality saturation, inspired from
the [egg](https://egraphs-good.github.io/) library for Rust. You can read more
about the design of the EGraph data structure and equality saturation algorithm
in the [egg paper](https://dl.acm.org/doi/pdf/10.1145/3434304).

See [Alessandro Cheli](https://0x0f0f0f.github.io/) and [Philip Zucker](https://www.philipzucker.com/)'s 
[talk at JuliaCon 2021](https://www.youtube.com/watch?v=tdXfsTliRJk) for an overview of the concepts introduced in this chapter of the manual (**NOTE**: Syntax in the talk slideshow is out of date).

## What can I do with EGraphs in Metatheory.jl?

In classical term rewriting, rewrites are typically destructive and forget the
matched left-hand side. Therefore, rules are applied in an arbitrary or
controlled order - this often results in local minima and looping. For decades,
programmers and scientists using term rewriting systems have spent their time
trying to find confluent and terminating systems of rules. This requires a lot
of effort and time. When studying any computational, mathematical or scientific
system governed by equational rules, about non obviously oriented equations, such as `(a + b) + c = a + (b + c
)`?

E-Graphs come to our help. 
EGraphs are bipartite graphs of [ENode](@ref)s and [EClass](@ref)es:
a data structure for efficiently represent and rewrite on many equivalent expressions at the same time. A sort of fast data structure for sets of trees. Subtrees and parents are shared if possible. This makes EGraphs similar to DAGs.
Most importantly, with EGraph rewriting you can use **bidirectional rewrite rules**, such as **equalities** without worrying about
the ordering and confluence of your rewrite system!
Therefore, rule application in EGraphs is non-destructive - everything is
copied! This allows users to run non-deterministic rewrite systems. Many rules
can match at the same time and the previous state of expressions will not be
lost.

The EGraph backend for Metatheory.jl allows you to create an
EGraph from a starting expression, to add more expressions to the EGraph with
`addexpr!`, and then to effectively fill the EGraph with all possible equivalent
expressions resulting from applying rewrite rules from a [theory](../rewrite#Theories), by using the
`saturate!` function. You can then easily extract expressions from an e-graph by calling `extract!` with a cost
function.

A killer feature of [egg](https://egraphs-good.github.io/) and Metatheory.jl
are **EGraph Analyses**. They allow you to annotate expressions and equivalence classes in an EGraph with values from a semilattice domain, and then to:
* Automatically extract optimal expressions from an EGraph deciding from analysis data.
* Have conditional rules that are executed if some criteria is met on analysis data
* Have dynamic rules that compute the right hand side based on analysis data.

## Library

The `Metatheory.Library` module contains utility functions and macros for creating
rules and theories from commonly used algebraic structures and
properties, to be used with the e-graph backend.
```julia
using Metatheory.Library

comm_monoid = @commutative_monoid (*) 1
```


#### Theories and Algebraic Structures

**The e-graphs backend can directly handle associativity, equalities
commutativity and distributivity**, rules that are
otherwise known of causing loops and require extensive user reasoning 
in classical rewriting.

```julia
t = @theory a b c begin
    a * b == b * a
    a * 1 == a
    a * (b * c) == (a * b) * c
end
```


## Equality Saturation

We can programmatically build and saturate an EGraph.
The function `saturate!` takes an `EGraph` and a theory, and executes
equality saturation. Returns a report
of the equality saturation process.
`saturate!` is configurable, customizable parameters include
a `timeout` on the number of iterations, a `eclasslimit` on the number of e-classes in the EGraph, a `stopwhen` functions that stops saturation when it evaluates to true.
```julia
g = EGraph(:((a * b) * (1 * (b + c))));
report = saturate!(g, t);
# access the saturated EGraph
report.egraph

# show some fancy stats
report
```

```
Equality Saturation Report
=================
        Stop Reason: saturated
        Iterations: 1
        EGraph Size: 9 eclasses, 51 nodes
 ───────────────────────────────────────────────────────────────────────────────────────
                                                Time                   Allocations      
                                        ──────────────────────   ───────────────────────
            Tot / % measured:                1.18s / 0.45%            955KiB / 68.1%    

 Section                        ncalls     time   %tot     avg     alloc   %tot      avg
 ───────────────────────────────────────────────────────────────────────────────────────
 Apply                               1   4.63ms  87.5%  4.63ms    512KiB  78.7%   512KiB
 Search                              1    656μs  12.4%   656μs    139KiB  21.3%   139KiB
   a * (b * c) == (a * b) * c        1    242μs  4.58%   242μs   79.2KiB  12.2%  79.2KiB
   a * b == b * a                    1    153μs  2.89%   153μs   34.2KiB  5.26%  34.2KiB
   a * 1 == a                        1    115μs  2.17%   115μs   14.4KiB  2.21%  14.4KiB
   appending matches                 3   4.06μs  0.08%  1.35μs      544B  0.08%     181B
 Rebuild                             1   3.75μs  0.07%  3.75μs     0.00B  0.00%    0.00B
 ───────────────────────────────────────────────────────────────────────────────────────
```

With the EGraph equality saturation backend, Metatheory.jl can prove **simple**
equalities very efficiently. The `@areequal` macro takes a theory and some
expressions and returns true iff the expressions are equal according to the
theory. The following example may return true with an appropriate example theory. 

```julia 
julia> @areequal some_theory (x+y)*(a+b) ((a*(x+y))+b*(x+y)) ((x*(a+b))+y*(a+b)) 
```


## Configurable Parameters

[`EGraphs.saturate!`](@ref) can accept an additional parameter of type
[`EGraphs.SaturationParams`](@ref) to configure the equality saturation algorithm.
Extensive documentation for the configurable parameters is available in the [`EGraphs.SaturationParams`](@ref) API docstring.

```julia
# create the saturation params
params = SaturationParams(timeout=10, eclasslimit=4000)
saturate!(egraph, theory, params)
```


```@meta
CurrentModule = Base
```

## Outline of the Equality Saturation Algorithm

The `saturate!` function behaves as following.
Given a starting e-graph `g`, a set of rewrite rules `t` and some parameters `p` (including an iteration limit `n`):
* For each rule in `t`, search through the e-graph for l.h.s.
* For each match produced, apply the rewrite
* Do a bottom-up traversal of the e-graph to rebuild the congruence closure
* If the e-graph hasn’t changed from last iteration, it has saturated. If so, halt saturation.
* Loop at most n times.

Note that knowing if an expression with a set of rules saturates an e-graph or never terminates
is still an open research problem



## Extracting from an EGraph

Since e-graphs non-deterministically represent many equivalent symbolic terms,
extracting an expression from an EGraph is the process of selecting and
extracting a single symbolic expression from the set of all the possible
expressions contained in the EGraph. Extraction is done through the `extract!`
function, and the theoretical background behind this procedure is an [EGraph
Analysis](https://dl.acm.org/doi/pdf/10.1145/3434304); A cost function is
provided as a parameter to the `extract!` function. This cost function will
examine mostly every e-node in the e-graph and will determine which e-nodes will
be chosen from each e-class through an automated, recursive algorithm.

Metatheory.jl already provides some simple cost functions, such as `astsize`,
which expresses preference for the smallest expressions contained in equivalence
classes.

Here's an example
Given the theory:

```@example extraction
using Metatheory
using Metatheory.Library

comm_monoid = @commutative_monoid (*) 1;
t = @theory a b c begin
    a + 0 --> a
    a + b --> b + a
    a + inv(a) --> 0 # inverse
    a + (b + c) --> (a + b) + c
	a * (b + c) --> (a * b) + (a * c)
	(a * b) + (a * c) --> a * (b + c)
	a * a --> a^2
	a --> a^1
	a^b * a^c --> a^(b+c)
	log(a^b) --> b * log(a)
	log(a * b) --> log(a) + log(b)
	log(1) --> 0
	log(:e) --> 1
	:e^(log(a)) --> a
	a::Number + b::Number => a + b
	a::Number * b::Number => a * b
end
t = comm_monoid ∪ t ;
nothing # hide
```

We can extract an expression by using

```@example extraction

expr = :((log(e) * log(e)) * (log(a^3 * a^2)))
g = EGraph(expr)
saturate!(g, t)
ex = extract!(g, astsize)
```

The second argument to `extract!` is a **cost function**. [astsize](@ref) is 
a cost function provided by default, which computes the size of expressions.


## Defining custom cost functions for extraction.

A *cost function* for *EGraph extraction* is a function used to determine
which *e-node* will be extracted from an *e-class*. 

It must return a positive, non-complex number value and, must accept 3 arguments.
1) The current [ENode](@ref) `n` that is being inspected. 
2) The current [EGraph](@ref) `g`.
3) The current analysis type `an`.

From those 3 parameters, one can access all the data needed to compute
the cost of an e-node recursively.

* One can use [TermInterface.jl](https://github.com/JuliaSymbolics/TermInterface.jl) methods to access the operation and child arguments of an e-node: `operation(n)`, `arity(n)` and `arguments(n)`
* Since e-node children always point to e-classes in the same e-graph, one can retrieve the [EClass](@ref) object for each child of the currently visited enode with `g[id] for id in arguments(n)`
* One can inspect the analysis data for a given eclass and a given analysis type `an`, by using [hasdata](@ref) and [getdata](@ref).
* Extraction analyses always associate a tuple of 2 values to a single e-class: which e-node is the one that minimizes the cost
and its cost. More details can be found in the [egg paper](https://dl.acm.org/doi/pdf/10.1145/3434304) in the *Analyses* section. 

Here's an example:

```julia
# This is a cost function that behaves like `astsize` but increments the cost 
# of nodes containing the `^` operation. This results in a tendency to avoid 
# extraction of expressions containing '^'.
function cost_function(n::ENodeTerm, g::EGraph, an::Type{<:AbstractAnalysis})
    cost = 1 + arity(n)

    operation(n) == :^ && (cost += 2)

    for id in arguments(n)
        eclass = g[id]
        # if the child e-class has not yet been analyzed, return +Inf
        !hasdata(eclass, an) && (cost += Inf; break)
        cost += last(getdata(eclass, an))
    end
    return cost
end

# All literal expressions (e.g `a`, 123, 0.42, "hello") have cost 1
cost_function(n::ENodeLiteral, g::EGraph, an::Type{<:AbstractAnalysis}) = 1
```

## EGraph Analyses

An *EGraph Analysis* is an efficient and automated way of analyzing all the possible
terms contained in an e-graph. Metatheory.jl provides a toolkit to ease and 
automate the process of EGraph Analysis. An *EGraph Analysis* defines a domain
of values and associates a value from the domain to each [EClass](@ref) in the graph.
Theoretically, the domain should form a [join semilattice](https://en.wikipedia.org/wiki/Semilattice).
Rewrites can cooperate with e-class analyses by depending on analysis facts and adding
equivalences that in turn establish additional facts. 

In Metatheory.jl, EGraph Analyses are identified by a *type* that is subtype of `AbstractAnalysis`.
An [`EGraph`](@ref) can only contain one analysis per type.
The following functions define an interface for analyses based on multiple dispatch 
on `AbstractAnalysis` types: 
* [islazy](@ref) should return true if the analysis should NOT be computed on-the-fly during egraphs operation, only when required.  
* [make](@ref) should take an ENode and return a value from the analysis domain.
* [join](@ref) should return the semilattice join of two values in the analysis domain (e.g. *given two analyses value from ENodes in the same EClass, which one should I choose?*)
* [modify!](@ref) Can be optionally implemented. Can be used modify an EClass on-the-fly given its analysis value.

### Defining a custom analysis

In this example, we will provide a custom analysis that tags each EClass in an EGraph
with `:even` if it contains an even number or with `:odd` if it represents an odd number,
or `nothing` if it does not contain a number at all. Let's suppose that the language of the symbolic expressions
that we are considering will contain *only integer numbers, variable symbols and the `*` and `+` operations.*

Since we are in a symbolic computation context, we are not interested in the
the actual numeric result of the expressions in the EGraph, but we only care to analyze and identify
the symbolic expressions that will result in an even or an odd number.

Defining an EGraph Analysis is similar to the process [Mathematical Induction](https://en.wikipedia.org/wiki/Mathematical_induction).
To define a custom EGraph Analysis, one should start by defining a type that 
subtypes `AbstractAnalysis` that will be used to identify this specific analysis and 
to dispatch against the required methods.

```julia
using Metatheory
using Metatheory.EGraphs
abstract type OddEvenAnalysis <: AbstractAnalysis end
```

The next step, the base case of induction, is to define a method for
[make](@ref) dispatching against our `OddEvenAnalysis`. First, we want to
associate an analysis value only to the *literals* contained in the EGraph. To do this we
take advantage of multiple dispatch against `ENodeLiteral`.

```julia
function EGraphs.make(an::Type{OddEvenAnalysis}, g::EGraph, n::ENodeLiteral)
    if n.value isa Integer
        return iseven(n.value) ? :even : :odd
    else 
        return nothing
    end
end
```

Now we have to consider the *induction step*. 
Knowing that our language contains only `*` and `+` operations, and knowing that:
* odd * odd = odd
* odd * even = even
* even * even = even

And we know that 
* odd + odd = even 
* odd + even = odd 
* even + even = even

We can now define a method for `make` dispatching against 
`OddEvenAnalysis` and `ENodeTerm`s to compute the analysis value for *nested* symbolic terms. 
We take advantage of the methods in [TermInterface](https://github.com/JuliaSymbolics/TermInterface.jl) 
to inspect the content of an `ENodeTerm`.
From the definition of an [ENode](@ref), we know that children of ENodes are always IDs pointing
to EClasses in the EGraph.

```julia
function EGraphs.make(an::Type{OddEvenAnalysis}, g::EGraph, n::ENodeTerm)
    # Let's consider only binary function call terms.
    if exprhead(n) == :call && arity(n) == 2
        op = operation(n)
        # Get the left and right child eclasses
        child_eclasses = arguments(n)
        l = g[child_eclasses[1]]
        r = g[child_eclasses[2]]

        # Get the corresponding OddEvenAnalysis value of the children
        # defaulting to nothing 
        ldata = getdata(l, an, nothing)
        rdata = getdata(r, an, nothing)

        if ldata isa Symbol && rdata isa Symbol
            if op == :*
                return (ldata == :even || rdata == :even) ? :even : :odd
            elseif op == :+
                return (ldata == rdata) ? :even : :odd
            end
        elseif isnothing(ldata) && rdata isa Symbol && op == :*
            return rdata
        elseif ldata isa Symbol && isnothing(rdata) && op == :*
            return ldata
        end
    end

    return nothing
end
```

We have now defined a way of tagging each ENode in the EGraph with `:odd` or `:even`, reasoning 
inductively on the analyses values. The [analyze!](@ref) function will do the dirty job of doing 
a recursive walk over the EGraph. The missing piece, is now telling Metatheory.jl how to merge together
analysis values. Since EClasses represent many equal ENodes, we have to inform the automated analysis
how to extract a single value out of the many analyses values contained in an EGraph.
We do this by defining a method for [join](@ref).

```julia
function EGraphs.join(an::Type{OddEvenAnalysis}, a, b)
    if a == b 
        return a 
    else
        # an expression cannot be odd and even at the same time!
        # this is contradictory, so we ignore the analysis value
        return nothing 
    end
end
```

We do not care to modify the content of EClasses in consequence of our analysis.
Therefore, we can skip the definition of [modify!](@ref).
We are now ready to test our analysis.

```julia
t = @theory a b c begin 
    a * (b * c) == (a * b) * c
    a + (b + c) == (a + b) + c
    a * b == b * a
    a + b == b + a
    a * (b + c) == (a * b) + (a * c)
end

function custom_analysis(expr)
    g = EGraph(expr)
    saturate!(g, t)
    analyze!(g, OddEvenAnalysis)
    return getdata(g[g.root], OddEvenAnalysis)
end

custom_analysis(:(3*a)) # :odd
custom_analysis(:(3*(2+a)*2)) # :even
custom_analysis(:(3y * (2x*y))) # :even
```
