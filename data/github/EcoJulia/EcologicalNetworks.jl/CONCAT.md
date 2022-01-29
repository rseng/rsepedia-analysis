# EcologicalNetworks.jl

This `julia` package provides a common interface to analyze all types of data
on ecological networks. It is designed to be general, easy to expand, and work
on bipartite/unipartite as well as deterministic/quantitative/probabilistic
networks. The current version is compatible with `julia` version 1.5 and up.

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://ecojulia.github.io/EcologicalNetworks.jl/stable)
[![](https://img.shields.io/badge/docs-dev-orange.svg)](https://ecojulia.github.io/EcologicalNetworks.jl/dev)

## Getting started

Install (from the package mode):

~~~
add EcologicalNetworks
~~~

Optionally, if you want to visualize networks, install

```
add EcologicalNetworksPlots
```

To access network data stored in [mangal.io](http://mangal.io), install

```
add Mangal
```

That's it. Now head over to the
[documentation](http://EcoJulia.github.io/EcologicalNetworks.jl/stable/).

## How's the code doing?

### Released version

[![DOI](https://zenodo.org/badge/25148478.svg)](https://zenodo.org/badge/latestdoi/25148478)
[![license](https://img.shields.io/badge/license-MIT%20%22Expat%22-yellowgreen.svg)](https://github.com/EcoJulia/EcologicalNetworks.jl/blob/master/LICENSE.md)

[![GitHub tag](https://img.shields.io/github/tag/EcoJulia/EcologicalNetworks.jl.svg)]()
[![GitHub issues](https://img.shields.io/github/issues/EcoJulia/EcologicalNetworks.jl.svg)]()
[![GitHub pull requests](https://img.shields.io/github/issues-pr/EcoJulia/EcologicalNetworks.jl.svg)]()

### On `master`

![CI](https://github.com/EcoJulia/EcologicalNetworks.jl/workflows/CI/badge.svg?branch=master)
![CI](https://github.com/EcoJulia/EcologicalNetworks.jl/workflows/TagBot/badge.svg?branch=master)
![CI](https://github.com/EcoJulia/EcologicalNetworks.jl/workflows/CompatHelper/badge.svg?branch=master)
[![codecov.io](http://codecov.io/github/EcoJulia/EcologicalNetworks.jl/coverage.svg?branch=master)](http://codecov.io/github/EcoJulia/EcologicalNetworks.jl?branch=master)

### On `develop`

![CI](https://github.com/EcoJulia/EcologicalNetworks.jl/workflows/CI/badge.svg?branch=develop)
![CI](https://github.com/EcoJulia/EcologicalNetworks.jl/workflows/TagBot/badge.svg?branch=develop)
![CI](https://github.com/EcoJulia/EcologicalNetworks.jl/workflows/CompatHelper/badge.svg?branch=develop)
[![codecov.io](http://codecov.io/github/EcoJulia/EcologicalNetworks.jl/coverage.svg?branch=develop)](http://codecov.io/github/EcoJulia/EcologicalNetworks.jl?branch=develop)
# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as contributors and maintainers pledge to making participation in our project and our community a harassment-free experience for everyone, regardless of age, body size, disability, ethnicity, gender identity and expression, level of experience, nationality, personal appearance, race, religion, or sexual identity and orientation.

## Our Standards

Examples of behavior that contributes to creating a positive environment include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a professional setting

## Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable behavior and are expected to take appropriate and fair corrective action in response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or reject comments, commits, code, wiki edits, issues, and other contributions that are not aligned to this Code of Conduct, or to ban temporarily or permanently any contributor for other behaviors that they deem inappropriate, threatening, offensive, or harmful.

## Scope

This Code of Conduct applies both within project spaces and in public spaces when an individual is representing the project or its community. Examples of representing a project or community include using an official project e-mail address, posting via an official social media account, or acting as an appointed representative at an online or offline event. Representation of a project may be further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by contacting the project team at timothee.poisot@umontreal.ca. The project team will review and investigate all complaints, and will respond in a way that it deems appropriate to the circumstances. The project team is obligated to maintain confidentiality with regard to the reporter of an incident. Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good faith may face temporary or permanent repercussions as determined by other members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4, available at [http://contributor-covenant.org/version/1/4][version]

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/4/
# Guidelines for contribution

This document described guidelines for contributions for projects led by members
of the [Poisot lab -- Quantitative and Computational Ecology][pl], Université de
Montréal. It should be followed by lab members (at all times), as well as people
with suggested added features, or willing to contribute enhancement or bug
fixes.

[pl]: http://poisotlab.io/

Edits to this document must be suggested as pull requests on the github project
located at <https://github.com/PoisotLab/PLCG>, and not made on the document
itself.

It is suggested to include this file to your projects via a `Makefile` rule:

``` make
CONTRIBUTING.md:
  wget -O $@ https://raw.githubusercontent.com/PoisotLab/PLCG/master/README.md
```

This work is licensed under a Creative Commons Attribution 4.0 International
License.

You should have received a copy of the license along with this
work. If not, see <http://creativecommons.org/licenses/by/4.0/>.

## General code of behaviour

We follow the [Contributor Covenant][cc] and (more pragmatically) [The No
Asshole Rule][nar]. No amazing technical or scientific contribution is an excuse
for creating a toxic environment.

[cc]: http://contributor-covenant.org/version/1/2/0/
[nar]: https://en.wikipedia.org/wiki/The_No_Asshole_Rule

## Issues

Reporting issues is the simplest and most efficient way to contribute. A useful
issue includes:

1. The version of all relevant software
2. Your operating system
3. A minimal reproducible example (the shortest code we can copy/paste to reproduce the issue)
4. A description of what was expected to happen
5. A description of what happened instead

Some repositories have issues templates enabled -- if so, they should be used.

## Authorship

Any contribution to the code of a software grants authorship on the next paper
describing this software, and immediate authorship to the next release (which
will receive a DOI). Any contribution to a paper hosted as a public
version-controlled repository grants authorship of the paper. For this reason,
it is important to correctly identify yourself. For `R` packages, this is done
through the `DESCRIPTION` file. Whenever a `CONTRIBUTORS` file is found, add
your name to it as part of the workflow described next. For authorship of
papers, we go for alphabetical order except for the first and last authors.

## Workflow

This section describes the general steps to make sure that your contribution is
integrated rapidly. The general workflow is as follows:

1. Fork the repository (see *Branches, etc.* below)
2. Create an *explicitly named branch*
3. Create a pull request *even if you haven't pushed code yet*
4. Be as explicit as possible on your goals
5. Do not squash / rebase commits while you work -- we will do so when merging

### Pull requests

Creating a pull request *before* you push any code will signal that you are
interested in contributing to the project. Once this is done, push often, and be
explicit about what the commits do (see commits, below). This gives the
opportunity for feedback during your work, and allow for tweaks in what you are
doing.

A *good* pull request (in addition to satisfying to all of the criteria below)
is:

1. Single purpose - it should do one thing, and one thing only
2. Short - it should ideally involve less than 250 lines of code
3. Limited in scope - it should ideally not span more than a handful of files
4. Well tested and well documented
5. Written in a style similar to the rest of the codebase

This will ensure that your contribution is rapidly reviewed and evaluated.

### Branches, etc.

The *tagged* versions of anything on `master` are stable releases. The `master`
branch itself is the latest version, which implies that it might be broken at
times. For more intensive development, there is usually a `dev` branch, or
feature-specific branches. All significant branches are under continuous
integration *and* code coverage analysis.

### Of emojis

[Atom's][atom] guideline suggest the use of emojis to easily identify what is
the purpose of each commit. This is a good idea and should be followed, and it
also saves a few characters from the commit first line. Specifically, prepend
your commits as follow (adapted from [these guidelines]).

| If the commit is about... | ...then use        | Example                                        |
|:--------------------------|:-------------------|:-----------------------------------------------|
| Work in progress          | `:construction:`   | :construction: new graphics                    |
| Bug fix                   | `:bug:`            | :bug: mean fails if NA                         |
| Code maintenance          | `:wrench:`         | :wrench: fix variable names                    |
| New test                  | `:rotating_light:` | :rotating_light: wget JSON resource            |
| New data                  | `:bar_chart:`      | :bar_chart: example pollination network        |
| New feature               | `:sparkles:`       | :sparkles: (anything amazing)                  |
| Documentation             | `:books:`          | :books: null models wrapper                    |
| Performance improvement   | `:racehorse:`      | :racehorse: parallelizes null model by default |
| Upcoming release          | `:package:`        | :package: v1.0.2                               |

[atom]: https://github.com/atom/atom/blob/master/CONTRIBUTING.md
[these guidelines]: https://github.com/dannyfritz/commit-message-emoji

## Tests and coverage

Most of our repositories undergo continuous integration, and code coverage
analysis.

It is expected that:

1. Your commits will not break a passing repo
2. Your commits *can* make a breaking repo pass
3. Your additions or changes will be adequately tested

This information will be given in the thread of your pull request. Most of our
repositories are set-up in a way that new commits that decrease coverage by more
than 10% won't be merged.

Good tests make sure that:

1. The code works as it should
2. The code breaks as it should (see *Program defensively* below)
3. Functions play nicely together

Before merging the content of *any* pull request, the following tests are done:

1. The branch to be merged *from* passes the build
2. The future branch (*i.e.* with the changeset) passes the build
3. The code coverage does not decrease by more than (usually) 10%

## Coding practices

### Use functions

Don't repeat yourself. If you have to do something more than once, write a
function for it. Functions should be as short as possible, and should be single
purpose. Think of every piece of code you write as an element of a library.

### Program defensively

Check user input, check data structure, and fail often but explicitly.

Bad:

``` julia
function add(x, y)
  return x + y
end
```

Acceptable:

``` julia
function add(x, y)
  @assert typeof(x) == Int64
  @assert typeof(y) == Int64
  return x + y
end
```

Good:

``` julia
function add(x, y)
  if typeof(x) != Int64
    throw(TypeError("The first argument should be an Int64", "add", Int64, typeof(x)))
  end
  if typeof(y) != Int64
    throw(TypeError("The second argument should be an Int64", "add", Int64, typeof(y)))
  end
  return x + y
end
```

(Of course all three solutions are actually bad because they are not the
idiomatic way of doing this, but you get the general idea.)

### Use type-static functions

Any given function must *always* return an object of the same type.

A function like `y = (x) -> x>3?true:x` (returns `true` if `x` is larger than 3,
else returns `x`) is not acceptable, as it is hard to predict, hard to debug,
and hard to use.

### Comment

Comment your code -- comments should not explain what the code *does* (this is
for the documentation), but how it does it.

### Variables names

Be as explicit as possible. Have a look at the rest of the codebase before you
start. Using `i`, `j`, `k` for loops is, as usual, fine.

### Docstrings and documentation

Write some.

`R` packages must be compiled with `roxygen`, `python` code must have docstrings
for all documentable objects, and `Julia` functions must be documented using
base docstrings, and `jldoctests` wherever possible.

There are three levels of documentation: the API documentation (which will be
generated from the docstrings), the documentation for fellow developers (which
resides ideally entirely in the comments), and documentation for the end users
(which include use cases, *i.e.* how is the feature use in the real world). Your
pull request must include relevant changes and additions to the documentation.

### Versioning

We use [semantic versioning][sv] (`major`.`minor`.`patch`). Anything that adds
no new feature should increase the `patch` number, new non-API-breaking changes
should increase `minor`, and major changes should increase `major`. Any increase
of a higher level resets the number after it (*e.g*, `0.3.1` becomes `1.0.0` for
the first major release). It is highly recommended that you do not start working
on API-breaking changes without having received a go from us first.

[sv]: http://semver.org/
# {describe feature here}

<!-- Thank you for this pull request! Please replace any information between curly brackets {like so} with the relevant information. -->

{short description of what the PR does}

## Related issues

<!-- If relevant, please add links to the relevant issues, using `#00` -->

## Checklist

- [ ] The code is covered by unit tests
  - [ ] Additional cases or module in `test/...`
  - [ ] Relevant lines in `test/runtests.jl`
- [ ] The code is *documented*
  - [ ] Docstrings written for every function
  - [ ] If needed, additional lines in `docs/src/...`
- [ ] All contributors have added their name and affiliation to `.zenodo.json`

## Pinging

Pinging @tpoisot
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: bug, need triage
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
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: enhancement, need triage
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
---
name: To-do
about: Short notes on development tasks
title: ''
labels: need triage, to do
assignees: ''

---

## What to do?

...

## Why?

...

## Any ideas how?

...
# EcologicalNetworks

This package provides a common interface for the analysis of ecological
networks, using `julia`. It is *very* opinionated about the "right" way
to do things, but we have documented our opinions in several publications
(see the references at the bottom of this page, and in the documentation of
all functions).

The package is built around a type system for networks, which is intended to
capture the different types of data and communities ecologists need to handle.
This makes the package extensible, both by writing additional methods with
a very fine-tuned dispatch, or by adding additional types that should work
out of the box (or be very close to).

This package is a *library* for the analysis of ecological networks. On purpose,
we do not provide "wrapper"-type functions that would perform an entire
analysis. We experimented with this idea during development, and rapidly
realized that even for the most simple research project, we needed to make small
tweaks that made the wrappers a nuisance. We decided to give you lego blocks,
and it's your job to build the kick-ass spaceship.

We tried to avoid making the package into yet another Domain Specific Language.
This means that when an operation should be expressed using the julian syntax,
we made it this way. Transforming networks from a type to another is done with
`convert`. Random networks are drawn with `rand`. Swapping of interactions
is done with `shuffle`. There is support for slicing of networks, as well
as the entire operations on sets. A lot of methods from `Base` have been
overloaded, and this *should* make the code easy to write and read, since
it looks almost exactly like any other *Julia* code on arrays.

### Why should I use this package?

It offers a single interface to analyse almost all type of networks for ecology.
It's somewhat fast (very specialized packages are likely to be faster). It's
built around the very best practices in network analysis. We think the type
system is very cool. It's very well tested and adequately documented. We used it
for research and teaching for months before releasing it. It's actively
maintained and we will keep adding functionalities.

You don't have to use it if you don't want to.

### But it doesn't even make figures!

The code for network visualization is in a companion package named
`EcologicalNetworksPlots`. There are two reasons for this decision.

First, network visualization, although attractive, is not necessary for network
analysis. It can help, but given the wrong network layout technique, it can also
introduce biases. When the volume of networks increased, we found that
visualization became less and less informative. Because it is not strictly
speaking a tool for analysis, it is not part of this package.

Second, it helps to keep software dependency small. Most of our work using this
package is done on clusters of one sort of the other, and having fewer
dependencies means that installation is easier. `EcologicalNetworksPlots` can be
installed like any other Julia package. It is also documented on [its own
website][ENP].

[ENP]: https://poisotlab.github.io/EcologicalNetworksPlots.jl/stable/

### And worse, you forgot my favorite method!

Yeah, about that. We probably didn't.

A lot of methods were considered for inclusion in the package, but ultimately
discarded because we were not 100% confident in their robustness, reliability,
validity, or interpretation. As we said, the package is *very* opinionated about
the right way to do things, and new functions require more time for maintenance
and testing; it makes sense for us to focus on things we trust.

If your favorite measure or method is missing, there are two solutions. First,
this package is essentially a library of functions to build network analyses, so
you can use this to create a function that does what you want. For example, if
you want to take the square root of a quantitative network, you can overload the
`√` method from base this way:

~~~ julia
import Base: √

function √(N::T) where {T <: QuantitativeNetwork}
   @assert all(N.edges .> zero(eltype(N.edges)))
   # Take the square root of the interaction strength
   sqrt_matrix = sqrt.(N.edges)
   # Return a new network with the correct types
   return T(sqrt_matrix, EcologicalNetworks._species_objects(N)...)
end
~~~

The second solution (which is actually a second *step* after you have been
writing your own method), is to submit a pull request to the package, to have
your new methods available in the next release. Currently, we will be very
selective about which methods are added (because every line of code needs to be
maintained).

## References

About the analysis of ecological networks in general, the package covers (or
will cover over time) most of the measures we identified as robust in the
following publication:

Delmas, Eva, Mathilde Besson, Marie-Hélène Brice, Laura A. Burkle, Giulio V.
Dalla Riva, Marie-Josée Fortin, Dominique Gravel, et al. « Analysing Ecological
Networks of Species Interactions ». Biological Reviews (2018), 112540.
https://doi.org/10.1111/brv.12433.


We highly recommend we keep it nearby when using the package. A lot of
decisions taken during development are grounded in the analysis of the
literature we conducted over a few years. Anything else is now documented
in the functions themselves.

## How can I contribute?

Good question!

The easiest way to contribute is to use the package, and [open an issue][issue]
whenever you can't manage to do something, think the syntax is not clear, or
the documentation is confusing. This is in fact one of the best ways to help.

[issue]: https://github.com/PoisotLab/EcologicalNetworks.jl/issues

If you want to contribute code, you can fork this repository, and start adding
the functions you want, or changing the code. Please work from the `develop`
branch (`master` does not accept pull requests except from maintainers, and
cannot be pushed to unless a series of conditions are met). It's better if all
of your code is tested and documented, but we will work with you when receiving
the pull request anyways.
# Public methods

```@autodocs
Modules = [EcologicalNetworks]
Private = false
Public = true
Order = [:function]
```
# Internals

```@autodocs
Modules = [EcologicalNetworks]
Private = true
Public = false
Order = [:function]
```
# Network β-diversity

Measures of β-diversity work by first calculating the unique/shared items (using
the `βs`, `βos`, and `βwn` functions), then passing on these arguments to one of
the `KGLXX` functions to return a (dis)similarity score. The `KGL` functions are
named for Koleff, Gaston, and Lennon -- the number of each function matches the
number in Table 1.

## β-diversity components

The package implements functions for the βs, βos, and βwn components of network
dissimilarity. In the original publication, we also described βst, which was the
proprotion of dissimilarity due to species turnover, and defined as βst = βwn -
βos *for measures of dissimilarity bounded between 0 and 1*. After discussing
with colleagues and considering our own use-cases, it appears that the
interpretation of βst is not always straightforward, and so we have decided to
exclude it form the available functions.

```@docs
βs
βos
βwn
```

## β-diversity measures

```@docs
KGL01
```

## Basic operations on networks

Internally, the functions for β-diversity rely on the usual operations on sets.
The act of combining two networks, for example, is a `union` operation.

```@docs
setdiff
union
intersect
```
# Modularity

The analysis of network modularity is done in three steps:

1. generate a starting point, using one of the starter functions
2. optimize modularity
3. analyse the output

All starter functions take a network as input, and return a tuple of this
network and a dictionary where every species maps onto its module. This forms
the input of all other modularity related functions.

## Measures

```@docs
Q
Qr
```

## Starters

```@docs
n_random_modules
each_species_its_module
lp
```

## Optimizers

```@docs
brim
salp
```

## Functional roles

```@docs
functional_cartography
```
# Nestedness

```@docs
η
nodf
ρ
```
# Links, degree, connectance

## Connectance and number of links

```@docs
sum
links
connectance
linkage_density
```

## Degree

```@docs
degree
```

```@docs
degree_var
```

## Species without interactions

```@docs
isdegenerate
simplify!
simplify
```

## Species-level specificity

```@docs
specificity
```
# Measures of overlap

```@docs
overlap
AJS
EAJS
```
# Motif enumeration

## List of canonical motifs

```@docs
unipartitemotifs
```

## Motif counting

```@docs
find_motif
```

## Probabilistic case

```@docs
expected_motif_count
```
# Paths and centrality

## Number of paths and shortest path

```@docs
number_of_paths
shortest_path
bellman_ford
dijkstra
```

## Centrality measures

```@docs
centrality_degree
centrality_closeness
centrality_katz
centrality_harmonic
```
# Resilience

We provide the metrics proposed by Gao et al (2016) which summarize the global
behaviour of complex unipartite networks. The dynamics of a system of N components
(nodes/species) can follow the coupled nonlinear differential equation

$$\frac{\text{d}x_i}{\text{d}t} = F(x_i) + \sum_{j=1}^N A_{ij}G(x_i, x_j)$$

where the adjacency matrix $A$ captures the interaction between the
components. This system can be described in 1-D using an effective term

$$\frac{\text{d}x_\text{eff}}{\text{d}t} = F(x_\text{eff}) + \beta_\text{eff}G(x_\text{eff}, x_\text{eff})$$

with $\beta_\text{eff}$ a single resilience parameter which can capture
the effect of perturbing the system (node/link removal, weight change...).
This resilience parameter can be computed from an `AbstractUnipartiteNetwork`
using the functions `βeff` or `resilience`.

It can be shown that

$$\beta_\text{eff} = \langle s \rangle + \mathcal{S}  \mathcal{H}\,,$$

with

- $\langle s \rangle$ the average weighted degree (computed using `s_mean`),
- $\mathcal{S}$ the symmetry(computed using `symmetry`),
- $\mathcal{H}$ the heterogeneity (computed using `heterogeneity`).


> Goa, J., Barzael, B. and Barabási 2016. Universal resilience patterns in complex networks.
> Nature 530(7590), 307-312. doi:10.1038/nature16948

## Available functions

```@docs
s
σ_in
σ_out
symmetry
heterogeneity
resilience
```
# Indices based on information theory

Indices based on information theory, such as entropy, mutual information etc, can easily be computed. To this end, the ecological network is transformed in a bivariate distribution. This is done by normalizing the adjacency or incidence matrix to obtain a doubly stochastic matrix. The information theoretic indices are computed either from this matrix or directly from the ecological network. Note that when using an array is input, the functions do not perform any checks whether the matrix is normalized and nonnegative. When the input is an ecological network, the functions automatically convert the network to a normalized probability matrix.

One can compute individual indices or use the function `information_decomposition` which performs the entire decomposition at once. This decomposition yields for a given network the deviation of the marginal distributions of the species with the uniform distribution (quantifying the evenness), the mutual information (quantifying the specialisation) and the variance of information (quantifying the freedom and stability of the interactions). These indices satisfy the following balance equation for the top ($T$) and bottom ($B$) throphic level:

$$\log(nm) = D(B,T) + 2 I(B;T) + V(B;T)$$

$$\log(n) = D(B) + I(B;T) + H(B|T)$$

$$\log(m) = D(T) + I(B;T) + H(T|B)$$

Here, $n$ and $m$ are number of bottom and top species, respectively.

Indices can be calculated for the joint distribution, as well as for the marginal distributions of the two trophic levels (if applicable), by changing an optional argument `dim=1` of the function.

## Network conversion

```@docs
make_joint_distribution
```

## Indices

```@docs
entropy
conditional_entropy
mutual_information
variation_information
potential_information
diff_entropy_uniform
```

## Decomposition

```@docs
information_decomposition
```

## Effective interactions

```@docs
convert2effective
```

## References

Stock, M.; Hoebeke, L.; De Baets, B. « Disentangling the Information in Species Interaction Networks ».
Entropy 2021, 23, 703. https://doi.org/10.3390/e23060703
Networks are following the `AbstractArray` interface rather closely, with some
additional functionalities to get at interactions by using species names rather
than networks positions.

## Accessing elements

```@docs
getindex
```

## Changing elements

The ecological networks types are *all* mutable.

```@docs
setindex!
```

## Network size

```@docs
size
```This page presents the core functions to manipulate networks. Whenever possible,
the approach of `EcologicalNetworks` is to overload functions from `Base`.

## Accessing species

```@docs
species
```

## Accessing interactions

### Presence of an interaction

```@docs
has_interaction
interactions
```
### Random network samples

```@docs
sample
```

## Network utilities

### Species richness

```@docs
richness
```

### Changing network shape

```@docs
permutedims
nodiagonal
nodiagonal!
```

### Invert interactions

```@docs
Base.:!
```
Ecological networks types follow the iteration interface, allowing them to be
used within a loop. Iteration is done on *interactions* (species are already
arrays), and return a named tuple with various information.

## Network length

```@docs
length
```

## Iteration

```@docs
iterate
```Conversions between types are used to perform two usual operations: make a
bipartite network unipartite, and remove quantitative information. There are two
high-level functions which work by using the union types, and a series of
type-to-type functions (the later should be avoided, and exists only to make the
high-level functions work).

```@docs
convert
```
One feature of `EcologicalNetwork` which makes the rest of the package works is
the type system to represent networks. This is not the most enthralling reading,
but this pacge will walk you through the different options, and discuss how and
when to use them.

## Network representation

All networks types have a field `edges` to store the *sparse* adjacency matrix,
and fields `S`, or `T` and `B`, for species in unipartite and bipartite networks
respectively. `edges` is always a two-dimensional array (see below for more
information), where interactions go *from the rows*, *to the columns*. Network
types are **mutable**. Operations that will modify the network end with a `!`,
as is the julian convention.

Fields `S`, `T`, and `B` are one-dimensional arrays of `AllowedSpeciesTypes` --
they currently can be `String` or `Symbol`, and represent the species/nodes
names. Future allowed types will be added in later releases.

You should *never* have to manipulate the network by calling its fields
directly. The `species` function will give you access to the species, and the
network slicing operations (see later sections) will let you access subset of
the network / individual interactions / set of neighbours.

Network types are iterable: this is equivalent to calling the `interactions`
function on a network. On small networks, `interactions` is faster (but
allocates the whole memory at once). On large networks, it can be less true, and
using the iteration approach can save some time. The iteration protocol is the
same as for all other Julia collections:

~~~
for (int_number, interaction) in N
  @info "Interaction $(int_number) -- $(interaction)"
end
~~~

The objects returned by the iteration protocol are named tuples with fields `to`
and `from` (always), and can have additional fields `probability` and
`strength`.

### Partiteness

In unipartite networks, the adjacency matrix `edges` is square, and has as many
rows/columns as there are elements in `S`. This is always checked and enforced
upon construction of the object, so you *cannot* have a mismatch.

In bipartite networks, the matrix `edges` is not necessarily square, and has
dimensions equal to the lengths of `T` (rows) and `B` (columns). This too is
checked upon construction.

All elements in `S` *must* be unique (no duplicate node names). In addition, all
names in the union of `T` and `B` must be unique too (so that when a bipartite
network is cast to a unipartite one, the constraint on unique names in `S` is
respected). These constraints are enforced when constructing the object, and
will return explicit error messages if not met.

### Type of information

At all points, you can have a look at the types of the interactions and the
species objects -- the next entries in this documentation give additional
information about the types allowed.

```@docs
eltype
```

## Union types

All networks are grouped upon the `AbstractEcologicalNetwork` type:

```@docs
AbstractEcologicalNetwork
```

The type of nodes that are allowed is determined by the *non-exported*
`EcologicalNetworks._check_species_validity` function. To allow an additional type of
node, you can write the following:

~~~ julia

struct Foo
  name::AbstractString
  bar::AbstractFloat
end

import EcologicalNetworks
function EcologicalNetworks._check_species_validity(::Type{Foo})
end
~~~

Note that **integers are never valid species identifiers**. By default, `String`
and `Symbol` are used. The function `_check_species_validity` should do *nothing*
for an accepted type (and it will throw an error for any other type).

### By partiteness

```@docs
AbstractBipartiteNetwork
AbstractUnipartiteNetwork
```

### By interaction type

```@docs
BinaryNetwork
QuantitativeNetwork
ProbabilisticNetwork
DeterministicNetwork
```


## List of available types

These are the types that you *actually* declare and use. They are presented last
because it is easier to understand what they are when you get a sense for the
different union types.

```@docs
UnipartiteNetwork
BipartiteNetwork
UnipartiteQuantitativeNetwork
BipartiteQuantitativeNetwork
UnipartiteProbabilisticNetwork
BipartiteProbabilisticNetwork
```
In this section, we will measure the dissimilarity between bipartite
host-parasite networks.

```@example betadiv
using EcologicalNetworks
using Plots
using Plots.PlotMeasures
```

We use networks that span the entirety of Eurasia. Because these networks are
originally quantitative, we will remove the information on interaction strength
using `convert`. Note that we convert to an union type (`BinaryNetwork`) -- the
`convert` function will select the appropriate network type to return based on
the partiteness. The core operations on sets (`union`, `diff`, and `intersect`)
are implemented for the `BinaryNetwork` type. As such, generating the "metaweb"
(*i.e.* the list of all species and all interactions in the complete dataset)
is:

```@example betadiv
all_hp_data = filter(x -> occursin("Hadfield", x.Reference), web_of_life());
ids = getfield.(all_hp_data, :ID);
networks = convert.(BinaryNetwork, web_of_life.(ids));
metaweb = reduce(union, networks)
```

From this metaweb, we can measure $\beta_{OS}'$, *i.e.* the dissimilarity of
every network to the expectation in the metaweb. Measuring the distance between
two networks is done in two steps. Dissimilarity is first partitioned into three
components (common elements, and elements unique to both samples), then the
value is measured based on the cardinality of these components. The functions to
generate the partitions are `βos` (dissimilarity of interactions between shared
species), `βs` (dissimilarity of species composition), and `βwn` (whole network
dissimilarity). The output of these functions is passed to one of the functions
to measure the actual $β$-diversity.

```@example betadiv
βcomponents = [βos(metaweb, n) for n in networks];
βosprime = KGL02.(βcomponents);
```

Finally, we measure the pairwise distance between all networks (because we use a
symmetric measure, we only need $n\times(n-1)$ distances):

```@example betadiv
S, OS, WN = Float64[], Float64[], Float64[]
for i in 1:(length(networks)-1)
  for j in (i+1):length(networks)
    push!(S, KGL02(βs(networks[i], networks[j])))
    push!(OS, KGL02(βos(networks[i], networks[j])))
    push!(WN, KGL02(βwn(networks[i], networks[j])))
  end
end
```

We can now visualize these data:

```@example betadiv
p1 = histogram(βosprime, frame=:origin, bins=20, c=:white, leg=false, grid=false, margin=10mm)
xaxis!(p1, "Difference to metaweb", (0,1))
yaxis!(p1, (0,10))

p2 = plot([0,1],[0,1], c=:grey, ls=:dash, frame=:origin, grid=false, lab="", legend=:bottomleft, margin=10mm)
scatter!(p2, S, OS, mc=:black, lab="shared sp.", msw=0.0)
scatter!(p2, S, WN, mc=:lightgrey, lab="all sp.", msw=0.0, m=:diamond)
xaxis!(p2, "Species dissimilarity", (0,1))
yaxis!(p2, "Network dissimilarity", (0,1))

plot(p1,p2, size=(700,300))
```In this example, we will show how the modular structure of an ecological network
can be optimized. Finding the optimal modular structure can be a time-consuming
process, as it relies on heuristic which are not guaranteed to converge to the
global maximum. There is no elegant alternative to trying multiple approaches,
repeating the process multiple times, and having some luck.

```@example modularity
using EcologicalNetworks
using EcologicalNetworksPlots
using Plots
using Plots.PlotMeasures
```

## Generating modular partitions

For the first approach, we will generate random partitions of the species across
3 to 12 modules, and evaluate 20 replicate attempts for each of these
combinations. The output we are interested in is the number of modules, and the
overall modularity.

```@example modularity
N = convert(BipartiteNetwork, web_of_life("M_PA_003"))

n = repeat(3:12, outer=20)
m = Array{Dict}(undef, length(n))

for i in eachindex(n)
  # Each run returns the network and its modules
  # We discard the network, and assign the modules to our object
  _, m[i] = n_random_modules(n[i])(N) |> x -> brim(x...)
end
```

```@example modularity
q = map(x -> Q(N,x), m);
c = (m .|> values |> collect) .|> unique .|> length;
```

```@example modularity
p1 = scatter(c, q, c=:grey, msw=0.0, leg=false, frame=:origin, grid=false, margin = 10mm)
xaxis!(p1, "Number of modules")
yaxis!(p1, "Modularity", (0, 0.5))
```

## Measuring modularity

Now that we have the modular partition for every attempt, we can count the
modules in it, and measure its modularity. Out of all attempts, we want to get
the most modular one, *i.e.* the one with highest modularity. In some simple
problems, there may be several partitions with the highest value, so we can
either take the first, or one at random:

```@example modularity
optimal = rand(findall(q.== maximum(q)));
best_m = m[optimal];
```

This can be plotted using `EcologicalNetworksPlots`:

```@example modularity
I = initial(RandomInitialLayout, N)
for step in 1:4000
  position!(SpringElectric(1.2; gravity=0.1), I, N)
end
p2 = plot(I, N, aspectratio=1)
scatter!(p2, I, N, bipartite=true, nodefill=best_m, markercolor=:isolum)
```

## Species functional roles

We can finally look at the functional roles of the species:

```@example modularity
roles = functional_cartography(N, best_m)
```

This function returns a tuple (an unmodifiable set of values) of coordinates for
every species, indicating its within-module contribution, and its participation
coefficient. These results can be plotted to separate species in module hubs,
network hubs, peripherals, and connectors. Note that in the context of
ecological networks, this classification is commonly used. It derives from
previous work on metabolic networks, which subdivides the plane in 7 (rather
than 4) regions. For the sake of completeness, we have added the 7 regions to
the plot as well.

```@example modularity
plot(Shape([-2, 2.5, 2.5, -2], [0, 0, 0.05, 0.05]), lab="", frame=:box, lc=:grey, opacity=0.3, c=:grey, lw=0.0, grid=false, margin = 10mm) #R1
plot!(Shape([-2, 2.5, 2.5, -2], [0.05, 0.05, 0.62, 0.62]), lab="", c=:transparent) #R2
plot!(Shape([-2, 2.5, 2.5, -2], [0.62, 0.62, 0.80, 0.80]), lab="", lc=:grey, opacity=0.3, c=:grey, lw=0.0) #R3
plot!(Shape([-2, 2.5, 2.5, -2], [0.80, 0.80, 1.0, 1.0]), lab="", c=:transparent) #R4

plot!(Shape([2.5, 3.0, 3.0, 2.5], [0, 0, 0.3, 0.3]), lab="", c=:transparent) #R5
plot!(Shape([2.5, 3.0, 3.0, 2.5], [0.3, 0.3, 0.75, 0.75]), lab="", lc=:grey, opacity=0.3, c=:grey, lw=0.0) #R6
plot!(Shape([2.5, 3.0, 3.0, 2.5], [0.75, 0.75, 1.0, 1.0]), lab="", c=:transparent) #R7

vline!([2.5], c=:black, ls=:dot, lw=2.0)
hline!([0.62], c=:black, ls=:dot, lw=2.0)
collect(values(roles)) |> x -> scatter!(x, leg=false, c=:white)
yaxis!("Among-module connectivity", (0,1))
xaxis!("Within-module degree", (-2, 3))
```
In this example, we will visualize different measures of centrality, using a
large food web.

```@example centr
using EcologicalNetworks
using EcologicalNetworksPlots
using Plots
using Plots.PlotMeasures
```

We will work on a quantitative intertidal food web:

```@example centr
N = convert(BinaryNetwork, pajek()[:ChesUpper])
```

We will then find a force-directed layout to visualize it:

```@example centr
I = initial(RandomInitialLayout, N)
for step in 1:200richness(N)
    position!(ForceAtlas2(1.2), I, N)
end
```

We can now show the different values of centrality on this plot:

## Degree centrality

```@example centr
plot(I, N)
scatter!(I, N, nodefill=centrality_degree(N), c=:YlGnBu, aspectratio=1)
```

## Katz centrality

```@example centr
plot(I, N)
scatter!(I, N, nodefill=centrality_katz(N), c=:YlGnBu, aspectratio=1)
```

## Closeness centrality

```@example centr
plot(I, N)
scatter!(I, N, nodefill=centrality_closeness(N), c=:YlGnBu, aspectratio=1)
```

## Harmonic centrality

```@example centr
plot(I, N)
scatter!(I, N, nodefill=centrality_harmonic(N), c=:YlGnBu, aspectratio=1)
```

## Eigenvector centrality

```@example centr
plot(I, N)
scatter!(I, N, nodefill=centrality_eigenvector(N), c=:YlGnBu, aspectratio=1)
```In this illustration, we will simulate extinctions of hosts, to show how the
package can be extended by using the core functions described in the "Interface"
section. Simply put, the goal of this example is to write a function to randomly
remove one host species, remove all parasite species that end up not connected
to a host, and measuring the effect of these extinctions on the remaining
network. Rather than measuring the network structure in the function, we will
return an array of networks to be manipulated later:

```@example ext
using EcologicalNetworks
using Plots
using Plots.PlotMeasures
using Statistics
```

```@example ext
function extinctions(N::T) where {T <: AbstractBipartiteNetwork}

  # We start by making a copy of the network to extinguish
  Y = [copy(N)]

  # While there is at least one species remaining...
  while richness(last(Y)) > 1
    # We remove one species randomly
    remain = sample(species(last(Y); dims=2), richness(last(Y); dims=2)-1, replace=false)

    # Remaining species
    R = last(Y)[:,remain]
    simplify!(R)

    # Then add the simplified network (without the extinct species) to our collection
    push!(Y, copy(R))
  end
  return Y
end
```

One classical analysis is to remove host species, and count the richness of
parasite species, to measure their robustness to host extinctions -- this is
usually done with multiple scenarios for order of extinction, but we will focus
on the random order here. Even though `EcologicalNetworks` has a built-in
function for richness, we can write a small wrapper around it:

```@example ext
function parasite_richness(N::T) where {T<:BinaryNetwork}
  return richness(N; dims=1)
end
```

Writing multiple functions that take a single argument allows to chain them in a
very expressive way: for example, measuring the richness on all timesteps in a
simulation is `N |> extinctions .|> parasite_richness`, or alternatively,
`parasite_richness.(extinctions(N))`. In @fig:extinctions, we illustrate the
output of this analysis on 100 simulations (average and standard deviation) for
one of the networks.

```@example ext
N = convert(BinaryNetwork, web_of_life("A_HP_050"))

X = Float64[]
Y = Float64[]
for i in 1:200
  timeseries = extinctions(N)
  path_l = parasite_richness.(timeseries)./richness(N; dims=1)
  prop_r = 1.0.-richness.(timeseries; dims=2)./richness(N; dims=2)
  append!(X, prop_r)
  append!(Y, path_l)
end
x = sort(unique(X))
y = zeros(Float64, length(x))
sy = zeros(Float64, length(x))
for (i, tx) in enumerate(x)
  y[i] = mean(Y[X.==tx])
  sy[i] = std(Y[X.==tx])
end

pl = plot(x, y, ribbon=sy, c=:black, fill=(:lightgrey), lw=2, ls=:dash, leg=false, margin = 10mm, grid=false, frame=:origin, xlim=(0,1), ylim=(0,1))
xaxis!(pl, "Proportion of hosts removed")
```
In this example, we will show how `EcologicalNetworks.jl` can be integrated with [`Mangal.jl`](https://github.com/EcoJulia/Mangal.jl) to analyse many ecological networks. Specifically, we will show how to analyse the association between meaningful network properties (i.e. connectance, nestedness, and modularity) using all food webs archived on the [`mangal.io`](https://mangal.io/#/) online database.

To conduct this analysis, we need to upload the following packages:

```@example mangal
using EcologicalNetworks
using Mangal
using DataFrames
using Plots
using Plots.PlotMeasures
```

We first retrieve relevant metadata for all 1,386 networks archived on `mangal.io` using the `Mangal.jl` package. We count the number of species $S$ and the total number of interactions $L$ in each network, as well as their number of trophic interactions (predation and herbivory). We store these information in a data frame along with the networks' ID numbers, and print the first 5 elements. Due to the high number of networks we handle, note that this step might take some time to run.

```@example mangal
number_of_networks = count(MangalNetwork)
count_per_page = 100
number_of_pages = convert(Int, ceil(number_of_networks/count_per_page))

mangal_networks = DataFrame(fill(Int64, 5),
                 [:id, :S, :L, :pred, :herb],
                 number_of_networks)

global cursor = 1
for page in 1:number_of_pages
    global cursor
    networks_in_page = Mangal.networks("count" => count_per_page, "page" => page-1)
    for current_network in networks_in_page
        S = count(MangalNode, current_network)
        L = count(MangalInteraction, current_network)
        pred = count(MangalInteraction, current_network, "type" => "predation")
        herb = count(MangalInteraction, current_network, "type" => "herbivory")
        mangal_networks[cursor,:] .= (current_network.id, S, L, pred, herb)
        cursor = cursor + 1
    end
end

first(mangal_networks, 5)
```

We now have all the information we need to identify all food webs archived on `mangal.io`. Here we consider as food webs any ecological networks mainly composed of trophic interactions. We find that 259 networks meet this condition.  

```@example mangal
foodwebs = mangal_networks[mangal_networks[!, :pred] .+ mangal_networks[!, :herb] ./ mangal_networks[!, :L] .> 0.5, :]

first(foodwebs, 5)
```

To analyse their properties, we first need to read all of these food webs using the `Mangal.jl` package, and then convert them to UnipartiteNetworks using the `EcologicalNetworks.jl` package. 

```@example mangal
mangal_foodwebs = network.(foodwebs.id)

unipartite_foodwebs = convert.(UnipartiteNetwork, mangal_foodwebs)

unipartite_foodwebs[1:5]
```

We can then compute any measure supported by `EcologicalNetworks.jl` for UnipartiteNetworks. In this example, we compute connectance, nestedness, and modularity. To compute network modularity, we use 100 random species assignments in 3 to 15 groups as our starters, the BRIM algorithm to optimize the modularity for each of these random partitions, and retain the maximum value for each food web. Readers are invited to take a look at the [documentation](https://ecojulia.github.io/EcologicalNetworks.jl/dev/properties/modularity/) for further details on how to compute modularity.

```@example mangal
foodweb_measures = DataFrame(fill(Float64, 3),
                 [:connect, :nested, :modul],
                 length(unipartite_foodwebs))

# connectance
foodweb_measures.connect = connectance.(unipartite_foodwebs)

# nestedness
foodweb_measures.nested = ρ.(unipartite_foodwebs)

# modularity (BRIM algorithm)
number_of_modules = repeat(3:15, outer=100)
modules = Array{Dict}(undef, length(number_of_modules))

for i in eachindex(unipartite_foodwebs)
    current_network = unipartite_foodwebs[i]
    for j in eachindex(number_of_modules)
        _, modules[j] = n_random_modules(number_of_modules[j])(current_network) |> x -> brim(x...)
    end
    partition_modularity = map(x -> Q(current_network,x), modules);
    foodweb_measures.modul[i] = maximum(partition_modularity)
end

first(foodweb_measures, 5)
```

The association between these food-web measures can then be plotted. We find that modularity is negatively associated with connectance and nestedness, whereas nestesdness and connectance are positively associated.

```@example mangal
pal=RGB(204/255,121/255,167/255)

scatter(foodweb_measures.connect, foodweb_measures.nested,
    alpha=0.6, color=pal,
    lab="", framestyle=:box,
    xlabel="Connectance",
    ylabel="Nestedness",
    margin = 10mm)
```

```@example mangal
scatter(foodweb_measures.connect, foodweb_measures.modul,
    alpha=0.6, color=pal,
    lab="",framestyle=:box,
    xlabel="Connectance",
    ylabel="Modularity",
    margin = 10mm)
```

```@example mangal
scatter(foodweb_measures.modul, foodweb_measures.nested,
    alpha=0.6, color=pal,
    lab="", framestyle=:box,
    xlabel="Modularity",
    ylabel="Nestedness",
    margin = 10mm)
```
# Null models

Randomization of networks is mostly used to perform null hypothesis significance
testing, or to draw random realizations of a probabilistic network. There are
two ways to perform networks randomization: either by shuffling interactions
within the networks while enforcing some constraints (`shuffle`) or by drawing
random samples from a probabilistic network (`rand`).

## Draw a network from a probabilistic network

```@docs
rand
```

## Generate probabilistic networks from deterministic networks

These functions generate a probabilistic network from a deterministic network,
where the probability of every interaction is determined by the degree
distribution (or connectance) of the network.

```@docs
null1
null2
null3
```

### Random Dot Product Graph model

The interaction probability in a Random Dot Product Graphs model are given by the dot product of species representations in two metric spaces of a given dimension (one describing species as consumers, one as producers -- or predators and preys).

```@docs
RDPG
```

The Random Dot Product Graph spaces are computed via a truncated Singular Value Decomposition of the food web adjacency matrix.

```@docs
svd_truncated
```

## Shuffle interactions

```@docs
shuffle!
shuffle
```
# Optimal transportation for species interaction networks

By solving an *optimal transportation* problem, one can estimate the interaction intensities given (1) a matrix of interaction utility values (i.e., the preferences for certain interactions between the species) and (2) the abundances of the top and bottom species. It hence allows predicting *how* species will interact. The interactions estimated intensities are given by the ecological coupling matrix $Q$, which has the optimal trade-off between utility (enriching for prefered interactions) and entropy (interactions are distributed as randomly as possible). The function `optimaltransportation` has the following inputs:
- a utility matrix `M`;
- the (relative) abundances of the top and bottom species `a` and `b`;
- the entropic regularization parameter `λ` (set default to 1).
For details, we refer to the paper presenting this framework.

```@docs
optimaltransportation
```

## References

Stock, M., Poisot, T., & De Baets, B. (2021). « Optimal transportation theory for
species interaction networks. » Ecology and Evolution, 00(1), ece3.7254.
https://doi.org/10.1002/ece3.7254# Structural network models

Structure of ecological networks is non-random. Network architecture can
have a strong effect on important ecosystem properties (Mougi and Kondoh
2012, Thébault and Fontaine 2010). Many of the structural features of
food-webs can be simulated using small number of simple rules. Despite this
simplicity these models can often accurately reproduce some of the second
order characteristics of empirical food-webs (Stouffer et al. 2005). These
characteristics of phenomenological stochastic models allow for their wide
applications e.g. to simulate biomass dynamics using dynamical models or
study extinction cascades.

> Mougi, A. and Kondoh, M. (2012) ‘Diversity of Interaction Types and
Ecological Community Stability’, Science, 337(6092), pp. 349–351. doi:
10.1126/science.1220529.

> Thébault, E. and Fontaine, C. (2010) ‘Stability of Ecological Communities
and the Architecture of Mutualistic and Trophic Networks’, Science,
329(5993), pp. 853–856. doi: 10.1126/science.1188321.

> Stouffer, D. B. et al. (2005) ‘Quantitative Patterns in the Structure
of Model and Empirical Food Webs’, Ecology, 86(5), pp. 1301–1311. doi:
10.1890/04-0957.

Many models with various interactions assignment algorithms have been
proposed. `EcologicalNetworks` provides functions to generate random ecological
networks of the `UnipartiteNetwork` type. Listed below are those most often
used in ecological studies.

## Cascade model

This model uses one abstract trophic trait. For any given consumer links can
be randomly assigned to a resource species with the trait value smaller than
that of a consumer.

```@docs
cascademodel
```

## Niche model

Niche model extended cascade model by introducing ranges for each consumer. In
this model consumers can predate on resources which trait values are within
the predators' 'niche' range.

```@docs
nichemodel
```

## Nested hierarchy model

In order to reproduce more faithfully properties of complex and
multidimensional natural nested hierarchy model tries to use simple rules
to incorporate also the phylogenetic similarity in resource composition
of predators.

```@docs
nestedhierarchymodel
```

## Minimum potential niche model

This model attempts to explicitly simulate forbidden links in empirical food webs.

```@docs
mpnmodel
```
