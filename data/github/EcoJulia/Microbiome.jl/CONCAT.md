![Microbiome.jl logo](logo.png)

[![status](https://joss.theoj.org/papers/450fa18f47932c5fd3b837edeac91440/status.svg)](https://joss.theoj.org/papers/450fa18f47932c5fd3b837edeac91440) ![EcoJulia maintainer: kescobo](https://img.shields.io/badge/EcoJulia%20Maintainer-kescobo-blue.svg)
## Latest Release 

[![Latest Release](https://img.shields.io/github/release/EcoJulia/Microbiome.jl.svg)](https://github.com/EcoJulia/Microbiome.jl/releases/latest)

[![Docs stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://docs.ecojulia.org/Microbiome.jl/stable/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/EcoJulia/Microbiome.jl/blob/master/LICENSE)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

## Development Status

[![Docs dev](https://img.shields.io/badge/docs-latest-blue.svg)](https://docs.ecojulia.org/Microbiome.jl/latest/)
[![CI](https://github.com/EcoJulia/Microbiome.jl/workflows/CI/badge.svg)](https://github.com/EcoJulia/Microbiome.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/EcoJulia/Microbiome.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/EcoJulia/Microbiome.jl)

## Description

Microbiome.jl is a package for manipulating and analyzing
microbiome and microbial community data.

## Installation

To use the latest version of `Microbiome.jl`,
you must be on `julia` v1.6 or greater.

Install `Microbiome.jl` from the Julia REPL:

```julia
julia> using Pkg

julia> pkg"add Microbiome"
```

If you are interested in the cutting edge of the development, please check out
the master branch to try new features before release.

```julia
julia> pkg"add Microbiome#main"
```

## Companion Packages

- You might also be interested in some functionality provided by
  [BiobakeryUtils.jl](https://github.com/EcoJulia/BiobakeryUtils.jl).
- Microbiome.jl uses [EcoBase](https://github.com/EcoJulia/EcoBase.jl) under the hood
  for many of its types and methods.
- Microbiome.jl has a `Tables.jl` interface, so you can convert `CommunityProfiles`
  to any `Tables.jl`-compatible sink.
  - `CommunityProfile` is not yet a sink itself
# News for Microbiome.jl

## v0.8.1

- Some changes for JOSS review ([#115](https://github.com/EcoJulia/Microbiome.jl/pull/115)),
  mostly in documentation
## v0.8.0

This release brings a host of changes,
including a manuscript for [submission to the Journal of Open Source Software](https://github.com/openjournals/joss-reviews/issues/3876).
Also, we're now part of [EcoJulia](http://github.com/EcoJulia)! (moved from [BioJulia](http://github.com/BioJulia)).

Notable changes include:

- Much better metadata handling in `CommunityProfile`s, (see esp issues [#63](https://github.com/EcoJulia/Microbiome.jl/issues/63) [#64](https://github.com/EcoJulia/Microbiome.jl/issues/64) [#65](https://github.com/EcoJulia/Microbiome.jl/issues/65) [#75](https://github.com/EcoJulia/Microbiome.jl/issues/75) [#76](https://github.com/EcoJulia/Microbiome.jl/issues/76))
- Add `commjoin` function to merge the contents of multiple `CommunityProfile`s ([#58](https://github.com/EcoJulia/Microbiome.jl/pull/58))
- Generic filtering of `CommunityProfile`s [#79](https://github.com/EcoJulia/Microbiome.jl/pull/79)

See [Release notes](https://github.com/EcoJulia/Microbiome.jl/releases/tag/v0.8.0) for full details.

## v0.7.0 Major overhaul

Lots of breaking changes in this one, including

- Whole new community profile structure using `AxisArrays` back-end
- New `AbstractSample` type, containing embedded metadata
- New `Taxon` and `GeneFunction` feature types
- Dropped dependencies on `SpatialEcology.jl` and `DataFrames.jl`,
  but added `Tables.jl` interface and `EcoBase.jl` interfaces
- Brought back convenience `braycurtis` and `pcoa` functions relying on
  `Distances.jl` and `MultivariateStats.jl` respectively

## v0.5.0

Major Changes

- Dropped Microbiome.jl-defined types and methods for DistanceMatrix. Use Distances.jl instead.
  [Requires this SpatialEcology PR](https://github.com/EcoJulia/SpatialEcology.jl/pull/36)
- Dropped Microbiome.jl-defined types and methods for PCoA. Use MDS from MultivariateStats.jl instead. [Requires this MultivariateStats PR](https://github.com/JuliaStats/MultivariateStats.jl/pull/85)
- Dropped Hclust leaf ordering. Added to Clustering.jl instead ([see Clustering.jl PR](https://github.com/JuliaStats/Clustering.jl/pull/170))

## v0.4.1

Major Changes:

- Dropped support for julia-0.7

## v0.3.0

Major changes

- Dropped support for julia-0.6, moved on to julia-0.7
- Moved Plotting and BioBakery utilities out to separate packages to reduce dependencies
  - These two new packages will be released shortly
- I'm keeping DataFrames dependency for now, though I'd like to move to IterableTables soonish (though constructors might be more sense in EcoBase)

## Older versions

I didn't have NEWS.md until v0.3. I doubt anyone cares, but I kinda tried to
keep track in release notes.
# Contributing to Microbiome.jl and BiobakeryUtils.jl

:+1::tada: First off, thanks for taking the time to contribute! :tada::+1:

## Etiquette and conduct

Contributors and maintainers for Microbiome.jl and BiobakeryUtils.jl
are expected to abide by the [Julia Community Standards](https://julialang.org/community/standards/)

## How can I contribute?

[(Adapted from the BioJulia's contributing guidelines)](https://github.com/BioJulia/Contributing)
### Reporting Bugs

If you follow the advice here, we will
better understand your report :pencil:, be able to reproduce the behaviour
:computer: :computer:, and identify related problems :mag_right:.

#### Before creating a bug report:

Please do the following:

1. Check the GitHub issue list for the package that is giving you problems.

2. If you find an issue already open for your problem, add a comment to let
  everyone know that you are experiencing the same issue.

3. If no **currently open** issue already exists for your problem that has already been
   then you should create a new issue.

   > **Note:** If you find a **Closed** issue that seems like it is the same thing
   > that you're experiencing, open a new issue and include a link to the original
   > issue in the body of your new one.

#### How to create a (good) new bug report:

Bugs are tracked as [GitHub issues](https://guides.github.com/features/issues/).

When you are creating a bug report, please do the following:

1. **Explain the problem**

   - *Use a clear and descriptive title* for the issue to identify the problem.
   - *Describe the exact steps which reproduce the problem* in as many details as possible.
     - Which function / method exactly you used?
     - What arguments or parameters were used?
     - *Provide a specific example*. (Includes links to pastebin, gists and so on.)
       If you're providing snippets in the issue, use
       [Markdown code blocks](https://help.github.com/articles/markdown-basics/#multiple-lines).

   - *Describe the behaviour you observed after following the steps*
     - Point out what exactly is the problem with that behaviour.
     - *Explain which behaviour you expected to see instead and why.*
     - *OPTIONALLY: Include screenshots and animated GIFs* which show you
       following the described steps and clearly demonstrate the problem.
       You can use [this tool](https://www.cockos.com/licecap/) to record GIFs on
       macOS and Windows, or [this tool](https://github.com/colinkeenan/silentcast)
       or [this tool](https://github.com/GNOME/byzanz) on Linux.

2. **Provide additional context for the problem (some of these may not always apply)**

   - *Did the problem start happening recently* (e.g. after updating to a new version)?
     - If the problem started recently, *can you reproduce the problem in older versions?*
     - Do you know the most recent package version in which the problem doesn't happen?

   - *Can you reliably reproduce the issue?* If not...
     - Provide details about how often the problem happens.
     - Provide details about under which conditions it normally happens.

   - Is the problem is related to *working with files*? If so....
     - Does the problem happen for all files and projects or only some?
     - Does the problem happen only when working with local or remote files?
     - Does the problem happen for files of a specific type, size, or encoding?
     - Is there anything else special about the files you are using?

3. **Include details about your configuration and environment**

- *Which version of the package are you using?*

- *What's the name and version of the OS you're using?*

- *Which julia packages do you have installed?*

- Are you using local configuration files to customize julia behaviour? If so...
  - Please provide the contents of those files, preferably in a
  [code block](https://help.github.com/articles/markdown-basics/#multiple-lines)
  or with a link to a [gist](https://gist.github.com/).


### Suggest an Enhancement

This section explains how to submit an enhancement proposal.
This includes completely new features, as well as minor improvements to existing functionality.
Following these suggestions will help maintainers and the community understand
your suggestion :pencil: and find related suggestions :mag_right:.

#### Before Submitting An Enhancement Proposal

1. **Perform a cursory issue search** to see if the enhancement has already been suggested.
2. If it has not, open a new issue as per the guidance below.
3. If it has...
  1. Add a comment to the existing issue instead of opening a new one.
  2. If it was closed, take the time to understand why this was so (it's ok to
     ask! :) ), and consider whether anything has changed that makes the reason
     outdated. If you can think of a convincing reason to reconsider the
     enhancement, feel free to open a new issue as per the guidance below.

#### How to submit a (good) new enhancement proposal

Enhancement proposals are tracked as
[GitHub issues](https://guides.github.com/features/issues/).

1. **Explain the enhancement**
   - *Use a clear and descriptive title* for the issue to identify the suggestion.
   - *Provide a step-by-step description of the suggested enhancement* in as many details as possible.
   - *Provide specific examples to demonstrate the steps*.
     Include copy/pasteable snippets which you use in those examples, as
     [Markdown code blocks](https://help.github.com/articles/markdown-basics/#multiple-lines).

   - If you want to change current behaviour...
     - Describe the *current* behaviour.
     - *Explain which behaviour you expected* to see instead and *why*.
     - *Will the proposed change alter APIs or existing exposed methods/types?*
       If so, this may cause dependency issues and breakages, so the maintainer
       will need to consider this when versioning the next release.

   - *OPTIONALLY: Include screenshots and animated GIFs*.
     You can use [this tool](https://www.cockos.com/licecap/) to record GIFs on
     macOS and Windows, and [this tool](https://github.com/colinkeenan/silentcast)
     or [this tool](https://github.com/GNOME/byzanz) on Linux.

2. **Provide additional context for the enhancement**

   - *Explain why this enhancement would be useful* to users and
     isn't something that can or should be implemented as a separate package.

   - *Do you know of other projects where this enhancement exists?*

3. **Include details about your configuration and environment**

   - Specify which *version of the package* you're using.

   - Specify the *name and version of the OS* you're using.

### Making Pull Requests

All julia packages can be developed locally.
For information on how to do this, see this section of the julia
[documentation](https://docs.julialang.org/en/stable/manual/packages/#Package-Development-1).

Before you start working on code, it is often a good idea to open an enhancement
[suggestion](#suggest-an-enhancement)

Once you decide to start working on code, the first thing you should do is make
yourself an account on [Github](https://github.com).
The chances are you already have one if you've done coding before and wanted to
make any scripts or software from a science project public.

The first step to contributing is to find the repository for the package you want to modify
(Microbiome.jl [is here](https://github.com/EcoJulia/Microbiome.jl),
BiobakeryUtils.jl [is here](https://github.com/EcoJulia/BiobakeryUtils.jl)).
Hit the 'Fork' button on the repositories page to create a forked copy of the
package for your own Github account. This is your blank slate to work on, and
will ensure your work and experiments won't hinder other users of the released
and stable package.

From there you can clone your fork of the package and work on it on your
machine using git.
Here's an example of cloning, assuming you already forked `Microbiome.jl`:

```sh
git clone https://github.com/<YOUR_GITHUB_USERNAME_HERE>/Microbiome.jl.git
```

Git will download or "clone" your fork and put it in a folder called
BioSequences.jl it creates in your current directory.

It is beyond the scope of this document to describe good git and github use in
more specific detail, as the folks at Git and GitHub have already done that wonderfully
on their own sites. If you have additional questions,
please feel free to start a [discussion on Microbiome.jl](https://github.com/EcoJulia/Microbiome.jl/discussions/new).

#### How to make (good) code contributions and new Pull-Requests

1. **In your code changes**

   - **Branch properly!**
     - If you are making a bug-fix, then you need to checkout your bug-fix branch
       from the last release tag.
     - If you are making a feature addition or other enhancement, checkout your
       branch from master.
     - See [here](#a-suggested-branching-model) for more information (or ask a package maintainer :smile:).

   - Follow the [julia style guide](https://docs.julialang.org/en/stable/manual/style-guide/).

   - Follow the [additional style suggestions](#additional-julia-code-style-suggestions).

   - Follow the [julia performance tips](https://docs.julialang.org/en/stable/manual/performance-tips/).

   - Update and add docstrings for new code, consistent with the [documentation styleguide](https://docs.julialang.org/en/stable/manual/documentation/).

   - Update information in the documentation located in the `docs/src/`
     folder of the package/repository if necessary.

   - Ensure that unit tests have been added which cover your code changes.

   - All changes should be compatible with the latest stable version of
     Julia.

   - Please comment liberally for complex pieces of internal code to facilitate comprehension.

2. **In your pull request**

   - *Describe* the changes in the pull request

   - Provide a *clear, simple, descriptive title*.

   - Do not include issue numbers in the PR title.

   - If you have implemented *new features* or behaviour
     - *Provide a description of the addition* in as many details as possible.
     - *Provide justification of the addition*.
     - *Provide a runnable example of use of your addition*. This lets reviewers
       and others try out the feature before it is merged or makes it's way to release.

   - If you have *changed current behaviour*...
     - *Describe the behaviour prior to you changes*
     - *Describe the behaviour after your changes* and justify why you have made the changes.
     - *Does your change alter APIs or existing exposed methods/types?*
       If so, this may cause dependency issues and breakages, so the maintainer
       will need to consider this when versioning the next release.
     - If you are implementing changes that are intended to increase performance, you
       should provide the results of a simple performance benchmark exercise
       demonstrating the improvement. Especially if the changes make code less legible.

#### Reviews and merging

You can open a pull request early on and push changes to it until it is ready,
or you can do all your editing locally and make a pull request only when it is
finished - it is up to you.

When your pull request is ready on Github,
mention one of the maintainers of the repo in a comment e.g. `@kescobo`
and ask them to review it.
You can also use Github's review feature.
They will review the code and documentation in the pull request,
and will assess it.

Your pull request will be accepted and merged if:

1. The dedicated package maintainers approve the pull request for merging.
2. The automated build system confirms that all unit tests pass without any issues.

It may also be that the reviewers or package maintainers will want to you to make
changes to your pull request before they will merge it.
Take the time to understand why any such request has been made,
and freely discuss it with the reviewers.
Feedback you receive should be constructive and considerate
(also see [here](#etiquette-and-conduct)).

## Styleguides

### Git Commit messages

* Use the present tense ("Add feature" not "Added feature").
* Use the imperative mood ("Move cursor to..." not "Moves cursor to...").
* Limit the first line to 72 characters or less.
* Reference issues and pull requests liberally after the first line.
* Consider starting the commit message with an applicable emoji:
    * :art: `:art:` when improving the format/structure of the code
    * :racehorse: `:racehorse:` when improving performance
    * :memo: `:memo:` when writing docs
    * :penguin: `:penguin:` when fixing something on Linux
    * :apple: `:apple:` when fixing something on macOS
    * :checkered_flag: `:checkered_flag:` when fixing something on Windows
    * :bug: `:bug:` when fixing a bug
    * :fire: `:fire:` when removing code or files
    * :green_heart: `:green_heart:` when fixing the CI build
    * :white_check_mark: `:white_check_mark:` when adding tests
    * :arrow_up: `:arrow_up:` when upgrading dependencies
    * :arrow_down: `:arrow_down:` when downgrading dependencies
    * :exclamation: `:exclamation:` when removing warnings or depreciations

### Additional julia style suggestions

- Indent with 4 spaces.

- For functions that are not a single expression, it is preferred to use an explicit `return`.
  Be aware that functions in julia implicitly return the the result of the last
  expression in the function, so plain `return` should be used to indicate that
  the function returns `nothing`.

- Type names are camel case, with the first letter capitalized. E.g.
  `SomeVeryUsefulType`.

- Module names should be camel case.

- Separate logical blocks of code with one blank line. Although it is common
  and acceptable for short single-line functions to be defined together on
  consecutive lines with no blank lines between them.

- Function names, apart from constructors, are all lowercase.
  Include underscores between words only if the name would be hard
  to read without.
  E.g.  `start`, `stop`, `find_letter` `find_last_digit`.
  It is good to separate concepts in a name with a `_`.

- Generally try to keep lines below 100-columns, unless splitting a long line
  onto multiple lines makes it harder to read.

- Files that declare modules should only declare the module, and import any
  modules that it requires. Any subsequent significant code should be included
  from separate files. E.g.

```julia
module AwesomeFeatures

using IntervalsTrees, JSON

include("feature1.jl")
include("feature2.jl")

end
```

- Files that declare modules should have the same name name of the module.
  E.g the module `SomeModule` is declared under the file `SomeModule.jl`.

- When extending method definitions, define the methods with a module name prefix. E.g.

```julia
function Base.start(iter::YourType)
  ...
end

Base.done(iter::YourType, state) = ...
```

- Functions that get or set variables in a struct should not be
  prefixed with 'get' or 'set'.
  The getter should be named for the variable it gets, and the setter
  should have the same name as the getter, with the suffix `!`.
  For example, for the variable `names`:

```julia
name(node) # get node name
name!(node, "somename") # set node name
```

- When using conditional branching, if code is statement-like, an
  if-else block should be used. However if the code is expression-like
  then julia's ternary operator should be used.
  ```julia
  matches == sketchlen ? 1.0 : matches / (2 * sketchlen - matches)
  ```
  Some simple checks and expressions are also expressed using the `&&` or `||`
  operators instead of if-else syntax. For example:
  ```julia
  isvalid(foo) || throw(ArgumentError("$foo is not valid"))
  ```
---
title: "Microbiome.jl and BiobakeryUtils.jl - Julia packages for working with microbial community data"
tags:
  - julia
  - biology
  - microbiome
  - microbial communities
authors:
  - name: Kevin S. Bonham^[Corresponding author]
    affiliation: 1
  - name: Annelle Abatoni Kayisire
    affiliation: 1
  - name: Anika S. Luo
    affiliation: 1
  - name: Vanja Klepac-Ceraj^[Corresponding author]
    affiliation: 1
affiliations:
  - name: Department of Biological Sciences, Wellesley College
    index: 1
date: 30 October 2021
bibliography: paper.bib
---

# Summary

`Microbiome.jl` is a julia package to facilitate analysis of microbial community data.
`BiobakeryUtils.jl` is built on top of `Microbiome.jl`,
and provides utilities for working with a suite of command line tools
(the bioBakery) that are widely used for converting raw metagenomic sequencing data
into tables of taxon and gene function counts.
Together, these packages provide an effective way to link microbial community data
with the power of julia's numerical, statistical, and plotting libraries.

# Statement of need

Complex microbial communities exist everywhere, including in and on the human body,
and have profound effects on the environment and human health [@LloydPrice2017].
Common methods for analyzing microbial communities (eg. 16S amplicon or metagenomic sequencing)
generate a large quantity of numerical data (eg. count or relative abundance data)
as well as metadata associated with biological samples (eg. locations, human subject data)
and microbial features (eg. taxa, gene functions) [@Mallick2017ExperimentalDA].

The julia programming language [@Bezanson2017-ud] is gaining increasing prominence in biological research
due to its speed and flexibility [@roesch2021julia],
and has a growing ecosystem of packages for working with biological and ecological data,
as well as libraries for Bayesian statistical analysis [@ge2018t],
scientific machine learning [@rackauckas2017differentialequations],
and plotting [@DanischKrumbiegel2021].
Julia's type system makes it incredibly easy for packages to interoperate,
making `Microbiome.jl` and `BiobakeryUtils.jl` an effective bridge between
microbial community data and julia's package ecosystem,
while remaining agnostic to downstream analysis.

# Functionality

At its most basic, microbial community data can be represented as a sparse matrix,
where one dimension is indexed by microbial `feature`s (eg. species),
and the other is indexed by biological `sample`s or observations (eg. a stool sample).
Together, the measured abundances of each `feature` in each `sample`
make up the taxonomic or function "profile."
Typically, additional information (`metadata`) about each `sample`
is also needed for downstream statistical analysis,
such as the location or human subject it was collected from,
data about that environment (salinity, temperature, etc. for environmental samples;
clinical covariates for human subjects),
and storage or processing details.
While the observed values for microbial `feature`s are uniformly numeric,
and can be efficiently stored in a sparse matrix of floating point numbers,
`metadata` can take many forms.
Further, `CommunityProfile`s may have hundreds to hundreds of thousands of features,
while typically only a few dozen metadata variables are necessary for a given analysis.

`Microbiome.jl` provides a convenient set of types and type constructors
to store and access this information (\autoref{fig1}).

- The `MicrobiomeSample` type contains `name` and `metadata` fields,
  and methods for efficiently adding and extracting stored metadata
- The `Taxon` type stores `name` and taxonomic `rank` (eg. `genus`, `phylum`) fields
- The `GeneFunction` type stores `name` and `taxon` fields,
  the later of which may be a `Taxon` (allowing taxonomically stratified gene functions).
- The `CommunityProfile` type is a wrapped `SparseMatrixCSC`,
  with `MicrobiomeSample`s as columns and features (`Taxon`s or `GeneFunction`s) as rows.
- `CommunityProfile`s can be indexed like normal julia arrays with integers,
  or with strings and regular expressions that will search on the `name`
  fields of the sample or feature dimensions.

Further, the `CommunityProfile` type implements the `Tables.jl` interface,
making it trivial to convert to other tabular representations,
in particular enabling round-tripping to and from column separated values (`.csv`) files
using `CSV.jl`.
Feature (`Taxon` and `GeneFunction`), `MicrobiomeSample`, and `CommunityProfile`
types are also implemented with the interface of `EcoBase.jl`,
potentially enabling integration with the wider EcoJulia family of packages.

![Concrete types provided by `Microbiome.jl` for storing information about features, samples, and whole communities.\label{fig1}](Microbiome-jl-fig1.png)

`BiobakeryUtils.jl` provides a julia interface for the command line utilities
from HUMAnN and MetaPhlAn, two widely-used tools
for using metagenomic sequencing reads to generate
functional and taxonomic profiles, respectively.
It also provides functionality to simplify installation of the tools
and I/O for the common file types used and produced by those tools.
Together, `Microbiome.jl` and `BiobakeryUtils.jl`
make it easy to load, manipulate, and analyze microbial community data (\autoref{fig2}).

![Microbial community analysis workflow using `Microbiome.jl` and `BiobakeryUtils.jl`\label{fig2}](Microbiome-jl-fig2.png)

# Limitations and future work

While `Microbiome.jl` and `BiobakeryUtils.jl` are already functional
and being used for research [@Tso2021-vv; @Lewis2021-be; @Peterson2021-mr],
there are several avenues for further development.

First, there are many additional tools in the bioBakery
whose interface and outputs could be incorporated into `BiobakeryUtils.jl`.
In particular, `StrainPhlAn` and `PanPhlAn` [@Beghini2021-xy],
which have tabular output distinct from but quite similar to that
of `HUMAnN` and `MetaPhlAn` could be supported.

Second, two of the largest plotting packages in the julia ecosystem,
`Plots.jl` and `Makie.jl` [@tom_breloff_2021_5566503; @DanischKrumbiegel2021]
share a common "recipes" system,
enabling package authors to provide instructions
for how to plot their types.
`Microbiome.jl` currently contains convenience functions to facilitate
the generation of easy-to-plot data structures,
but including plot recipes for things like ordinations (PCoA),
abundance bar plots, and other commonly used microbial community visualizations
would make it even easier to generate publication-quality figures.

Finally, better integration with EcoJulia would carry a host of benefits.
For example, `Diversity.jl` [@Diversity.jl-2016] provides a wide array
of alpha and beta diversity metrics that could be beneficial
for investigations of microbial diversity.
There are also several packages that provide functionality
around phylogenies and taxonomic information
that could enhance or replace `Taxon`,
making it easier to gain insight into the relationships
between different microbial taxa found in communities.

# Acknowledgements

The authors would like to thank the families participating in the RESONANCE cohort.
This work was funded in part by the NIH UG3 OD023313 (VK-C).

# References
```@meta
CurrentModule = Microbiome
DocTestSetup  = quote
    using Microbiome
    using Microbiome.SparseArrays

    open("taxprof.csv", "w") do io
        print(io, """
        features,s1,s2
        s__Escherichia_coli,1,0
        s__Bifidobacterium_longum,3,2
        s__Prevotella_copri,5,4
        """)
    end
end
```

# Working with microbial abundances

The primary type for working with microbial abundances is the [`CommunityProfile`](@ref),
which is a sparse matrix with [`MicrobiomeSample`](@ref)s as column indices
and features (eg [`Taxon`](@ref)s or [`GeneFunction`](@ref)s) as row indices.

For example, let's make a `CommunityProfile` with 3 samples,
5 species and 5 genera.

```jldoctest profiles
julia> samps = MicrobiomeSample.(["s1", "s2", "s3"])
3-element Vector{MicrobiomeSample}:
 MicrobiomeSample("s1", {})
 MicrobiomeSample("s2", {})
 MicrobiomeSample("s3", {})

julia> taxa = [[Taxon("s$i", :species) for i in 1:5]; [Taxon("g$i", :genus) for i in 1:5]]
10-element Vector{Taxon}:
 Taxon("s1", :species)
 Taxon("s2", :species)
 Taxon("s3", :species)
 Taxon("s4", :species)
 Taxon("s5", :species)
 Taxon("g1", :genus)
 Taxon("g2", :genus)
 Taxon("g3", :genus)
 Taxon("g4", :genus)
 Taxon("g5", :genus)

julia> using SparseArrays

julia> mat = spzeros(10, 3); # 10 x 3 matrix filled with zeros

julia> for i in 1:10, j in 1:3
           if all(isodd, (i,j))
               mat[i,j] = i+j
           elseif all(iseven, (i, j))
               mat[i,j] = i / j
           end
       end

julia> mat
10×3 SparseMatrixCSC{Float64, Int64} with 15 stored entries:
  2.0   ⋅    4.0
   ⋅   1.0    ⋅
  4.0   ⋅    6.0
   ⋅   2.0    ⋅
  6.0   ⋅    8.0
   ⋅   3.0    ⋅
  8.0   ⋅   10.0
   ⋅   4.0    ⋅
 10.0   ⋅   12.0
   ⋅   5.0    ⋅

julia> comm = CommunityProfile(mat, taxa, samps)
CommunityProfile{Float64, Taxon, MicrobiomeSample} with 10 features in 3 samples

Feature names:
s1, s2, s3...g4, g5

Sample names:
s1, s2, s3
```

## Accessing `CommunityProfile` contents

It is easy to get out the underlying abundance matrix,
features, and samples using
[`abundances`](@ref), [`features`](@ref), and [`samples`](@ref) respectively:

```jldoctest profiles
julia> abundances(comm)
10×3 SparseMatrixCSC{Float64, Int64} with 15 stored entries:
  2.0   ⋅    4.0
   ⋅   1.0    ⋅
  4.0   ⋅    6.0
   ⋅   2.0    ⋅
  6.0   ⋅    8.0
   ⋅   3.0    ⋅
  8.0   ⋅   10.0
   ⋅   4.0    ⋅
 10.0   ⋅   12.0
   ⋅   5.0    ⋅


julia> features(comm)
10-element Vector{Taxon}:
 Taxon("s1", :species)
 Taxon("s2", :species)
 Taxon("s3", :species)
 Taxon("s4", :species)
 Taxon("s5", :species)
 Taxon("g1", :genus)
 Taxon("g2", :genus)
 Taxon("g3", :genus)
 Taxon("g4", :genus)
 Taxon("g5", :genus)
 
julia> samples(comm)
3-element Vector{MicrobiomeSample}:
 MicrobiomeSample("s1", {})
 MicrobiomeSample("s2", {})
 MicrobiomeSample("s3", {})
```

You can also just get an array of feature and sample names
using [`featurenames`](@ref) and [`samplenames`](@ref) respectively.

```jldoctest profiles
julia> featurenames(comm)
10-element Vector{String}:
 "s1"
 "s2"
 "s3"
 "s4"
 "s5"
 "g1"
 "g2"
 "g3"
 "g4"
 "g5"

julia> samplenames(comm)
3-element Vector{String}:
 "s1"
 "s2"
 "s3"
```

Finally, you can pull out the metadata of all samples
or a subset using [`metadata`](@ref).
The returned value is a vector of `NamedTuple`s,
which is compliant with the [`Tables.jl`](https://github.com/JuliaData/Tables.jl) interface,
so it's easy to load into other formats (like [`DataFrames.jl`](https://github.com/JuliaData/DataFrames.jl) for example):

```jldoctest profiles
julia> metadata(comm)
3-element Vector{NamedTuple{(:sample,), Tuple{String}}}:
 (sample = "s1",)
 (sample = "s2",)
 (sample = "s3",)
```

Of course, for now, the metadata is pretty boring - 
jump ahead to [Working with metadata](@ref working-metadata)
to do some more fun things.

## Indexing and selecting

`CommunityProfile`s wrap a sparse matrix,
and you can access the values as you would a normal matrix.
In julia, you can pull out specific values using `[row, col]`.
So for example, to get the 3rd row, 2nd column, of matrix `mat`:

```jldoctest indexing
julia> mat = reshape(1:12, 4, 3) |> collect
4×3 Matrix{Int64}:
 1  5   9
 2  6  10
 3  7  11
 4  8  12

julia> mat[3,2]
7
```

You can also get "slices", eg to get rows 2-4, column 1:

```jldoctest indexing
julia> mat[2:4, 1]
3-element Vector{Int64}:
 2
 3
 4
```

To get all of one dimension, you can just use a bare `:`

```jldoctest indexing
julia> mat[:, 1:2]
4×2 Matrix{Int64}:
 1  5
 2  6
 3  7
 4  8
```

For `CommunityProfile`s, indexing with integer values
will return the value of the matrix at that position,
while indexing with slices will return a new `CommunityProfile`.
To get the values of a matrix slice, use `abundances` after indexing.

```jldoctest profiles
julia> comm[1,3]
4.0

julia> comm[6:8,3]
CommunityProfile{Float64, Taxon, MicrobiomeSample} with 3 features in 1 samples

Feature names:
g1, g2, g3

Sample names:
s3

julia> comm[1:3,3] |> abundances
3×1 SparseMatrixCSC{Float64, Int64} with 2 stored entries:
 4.0
  ⋅
 6.0
```

### Indexing with strings and regular expressions

It is often inconvenient to find the numerical index
of a particular feature or sample.
Instead, you can use strings or regular expressions
to get slices of a `CommunityProfile`,
which will match on the `name` field of the features or samples.
This kind of indexing always returns a `CommunityProfile`,
even if it only has 1 value.

```jldoctest profiles
julia> comm["g1", "s1"]
CommunityProfile{Float64, Taxon, MicrobiomeSample} with 1 features in 1 samples

Feature names:
g1

Sample names:
s1



julia> comm[r"[gs]1", "s1"]
CommunityProfile{Float64, Taxon, MicrobiomeSample} with 2 features in 1 samples

Feature names:
s1, g1

Sample names:
s1
```



## [Working with metadata](@id working-metadata)

The `metadata` dictionaries of `MicrobiomeSample`s can accessed
and updated inside a `CommunityProfile` using `insert!`, `delete!`, `set!`, and `unset!`.
either by directly acting on the underlying sample object,
or using special methods that take the sample name as the second argument:

```jldoctest profiles
julia> set!(samples(comm)[1], :subject, "kevin")
MicrobiomeSample("s1", {:subject = "kevin"})

julia> insert!(comm, "s2", :subject, "anika")
MicrobiomeSample("s2", {:subject = "anika"})
```

You can also retrieve all metadata associated with a table using [`metadata`](@ref),
which will return a `Tables.jl`-compatible vector or `NamedTuple`s,
where each "row" corresponds to one sample.
All metadata fields found in any sample will be returned in every row,
with the value `missing` in any samples that do not have that field set.


```jldoctest profiles
julia> metadata(comm)
3-element Vector{NamedTuple{(:sample, :subject)}}:
 (sample = "s1", subject = "kevin")
 (sample = "s2", subject = "anika")
 (sample = "s3", subject = missing)
```

And you can bulk-`insert!` or `set!` metadata by passing a similar Table-like object
with the a field (`:sample` by default) matching sample names found in the `CommunityProfile`.

```jldoctest profiles
julia> md = [(sample="s1", subject="kevin", foo="bar"), (sample="s3", subject="annelle", foo="baz")]
2-element Vector{NamedTuple{(:sample, :subject, :foo), Tuple{String, String, String}}}:
 (sample = "s1", subject = "kevin", foo = "bar")
 (sample = "s3", subject = "annelle", foo = "baz")

julia> set!(comm, md)

julia> md2 = [(name="s1", other="Hello, World!"), (name="s2", other="Goodbye!")]
2-element Vector{NamedTuple{(:name, :other), Tuple{String, String}}}:
 (name = "s1", other = "Hello, World!")
 (name = "s2", other = "Goodbye!")

julia> insert!(comm, md2; namecol=:name)

julia> metadata(comm)
3-element Vector{NamedTuple{(:sample, :subject, :foo, :other)}}:
 (sample = "s1", subject = "kevin", foo = "bar", other = "Hello, World!")
 (sample = "s2", subject = "anika", foo = missing, other = "Goodbye!")
 (sample = "s3", subject = "annelle", foo = "baz", other = missing)
```

## Importing from tabular data

The `CommunityProfile` type is compatible with `Tables.jl`,
so should be able to import data from tabular sources
with the help of some other packages.

For example, using the [`DataFrames.jl`](https://github.com/JuliaData/DataFrames.jl) package
(which you can install at the REPL with `] add DataFrames`):

```jldoctest
julia> using DataFrames, Microbiome

julia> df = DataFrame(features = ["gene$i|s__Some_species" for i in 1:5], s1 = 1:2:10, s2 = 0:2:9)
5×3 DataFrame
 Row │ features               s1     s2
     │ String                 Int64  Int64
─────┼─────────────────────────────────────
   1 │ gene1|s__Some_species      1      0
   2 │ gene2|s__Some_species      3      2
   3 │ gene3|s__Some_species      5      4
   4 │ gene4|s__Some_species      7      6
   5 │ gene5|s__Some_species      9      8

julia> gfs = genefunction.(df.features)
5-element Vector{GeneFunction}:
 GeneFunction("gene1", Taxon("Some_species", :species))
 GeneFunction("gene2", Taxon("Some_species", :species))
 GeneFunction("gene3", Taxon("Some_species", :species))
 GeneFunction("gene4", Taxon("Some_species", :species))
 GeneFunction("gene5", Taxon("Some_species", :species))

julia> mss = MicrobiomeSample.(names(df)[2:end])
2-element Vector{MicrobiomeSample}:
 MicrobiomeSample("s1", {})
 MicrobiomeSample("s2", {})

julia> CommunityProfile(Matrix(df[!, 2:end]), gfs, mss)
CommunityProfile{Int64, GeneFunction, MicrobiomeSample} with 5 features in 2 samples

Feature names:
gene1, gene2, gene3, gene4, gene5

Sample names:
s1, s2
```

Alternatively, with a CSV file, and the [`CSV.jl`](https://github.com/JuliaData/CSV.jl) package:

```jldoctest csvexample
julia> println.(eachline("taxprof.csv"));
features,s1,s2
s__Escherichia_coli,1,0
s__Bifidobacterium_longum,3,2
s__Prevotella_copri,5,4

julia> using CSV, CSV.Tables

julia> tbl = CSV.read("taxprof.csv", Tables.columntable);

julia> txs = taxon.(tbl[1])
3-element Vector{Taxon}:
 Taxon("Escherichia_coli", :species)
 Taxon("Bifidobacterium_longum", :species)
 Taxon("Prevotella_copri", :species)

julia> mss = [MicrobiomeSample(string(k)) for k in keys(tbl)[2:end]]
2-element Vector{MicrobiomeSample}:
 MicrobiomeSample("s1", {})
 MicrobiomeSample("s2", {})

julia> mat =  hcat([tbl[i] for i in 2:length(tbl)]...)
3×2 Matrix{Int64}:
 1  0
 3  2
 5  4

julia> CommunityProfile(mat, txs, mss)
CommunityProfile{Int64, Taxon, MicrobiomeSample} with 3 features in 2 samples

Feature names:
Escherichia_coli, Bifidobacterium_longum, Prevotella_copri

Sample names:
s1, s2
```

You may also be interested in [`BiobakeryUtils.jl`](https://github.com/EcoJulia/BiobakeryUtils.jl)
which has convenient functions for reading in file types generated by
the `bioBakery` suite of computational tools.

## Types and Methods

```@docs
CommunityProfile
samples
features
samplenames
**featurenames**
commjoin
relativeabundance
relativeabundance!
present
prevalence
prevalence_filter
```
# Microbiome.jl

## Description

Microbiome.jl is a package for manipulating and analyzing
microbiome and microbial community data.

## Installation

Install Microbiome from the Julia REPL:

```julia
julia> using Pkg

julia> pkg"add Microbiome"
```

If you are interested in the cutting edge of the development, please check out
the master branch to try new features before release.

```julia
julia> pkg"add Microbiome#master"
```

```@contents
Pages=["samples_features.md", "profiles.md", "diversity.md"]
Depth=2
``````@meta
CurrentModule = Microbiome
DocTestSetup  = quote
    using Microbiome
    using Microbiome.Dictionaries
end
```
# Samples and Features

Microbial profiles are made up of `AbstractSample`s and `AbstractFeature`s.
Typically, an `AbstractSample` is an individual biospecimen or other observation,
and contains some number of `AbstractFeature`s, such as taxa or gene functions.
`AbstractSample`s may also contain arbitrary metadata.

## Sample types

At its most basic, an `AbstractSample` simply encodes a `name`
(which should be a unique identifier), and a place to hold metadata.
The concrete type [`MicrobiomeSample`](@ref) is implemented with these two fields,
the latter of which is a `Dictionary` from [`Dictionaries.jl`](https://github.com/andyferris/Dictionaries.jl).


You can instantiate a `MicrobiomeSample` with just a name (in which case the metadata dictionary will be empty),
using keyword arguments for metadata entries,
or with existing metadata in the form of a dictionary (with keys of type `Symbol`) or a `NamedTuple`.

```jldoctest sampletypes
julia> s1 = MicrobiomeSample("sample1")
MicrobiomeSample("sample1", {})

julia> s2 = MicrobiomeSample("sample2"; age=37)
MicrobiomeSample("sample2", {:age = 37})

julia> s3 = MicrobiomeSample("sample3", Dict(:gender=>"female", :age=>23))
MicrobiomeSample("sample3", {:age = 23, :gender = "female"})
```

### Working with metadata

To change or add metadata, you can use the [same syntax](https://github.com/andyferris/Dictionaries.jl#accessing-dictionaries)
as working with a [`Dictionary`] directly,
though note that this is a bit different from the `Dict` type in base julia:

```jldoctest sampletypes
julia> insert!(s1, :age, 50)
MicrobiomeSample("sample1", {:age = 50})

julia> set!(s3, :gender, "nonbinary")
MicrobiomeSample("sample3", {:age = 23, :gender = "nonbinary"})

julia> delete!(s3, :gender)
MicrobiomeSample("sample3", {:age = 23})
```

You can access values of the dictionary using either `getindex`
or `getfield` syntax, that is:

```jldoctest sampletypes
julia> s3[:age]
23

julia> s3.age
23
```

Bulk addiction of metadata is also possible, by passing a `Dictionary` or `NamedTuple`
to `set!` or `insert!` (the latter will fail if any of the incoming keys are already found):

```jldoctest sampletypes
julia> insert!(s3, (gender = "nonbinary", genotype="XX"))
MicrobiomeSample("sample3", {:age = 23, :gender = "nonbinary", :genotype = "XX"})

julia> set!(s3, (genotype="XY", ses=7))
MicrobiomeSample("sample3", {:age = 23, :gender = "nonbinary", :genotype = "XY", :ses = 7})
```


## Feature Types

`AbstractFeature` types also have a `name`, but other fields are optional.
`Microbiome.jl` defines three concrete `AbstractFeature` types,
[`Taxon`](@ref), [`GeneFunction`](@ref), and [`Metabolite`](@ref).


### Taxon

The [`Taxon`](@ref) type contains a name and (optionally) a rank (eg `:phylum`).

```jldoctest taxon
julia> ecoli = Taxon("Escherichia_coli", :species)
Taxon("Escherichia_coli", :species)

julia> uncl = Taxon("Unknown_bug")
Taxon("Unknown_bug", missing)
```

You can access the name and rank fields using [`name`](@ref) and [`taxrank`](@ref) respectively, and also check whether the instance
has a rank with [`hasrank`](@ref), which returns `true` or `false`.

```jldoctest taxon
julia> hasrank(ecoli)
true

julia> hasrank(uncl)
false

julia> taxrank(ecoli)
:species

julia> taxrank(uncl)
missing

julia> name(ecoli)
"Escherichia_coli"

julia> name(uncl)
"Unknown_bug"
```

For compatibility with other tools, converting a `Taxon` to a `String`
will return the name prepended with the first letter of the taxonomic rank
and 2 underscores.
You can convert back using [`taxon`](@ref) (note the lowercase 't'):

```jldoctest taxon
julia> String(uncl)
"u__Unknown_bug"

julia> String(ecoli)
"s__Escherichia_coli"

julia> String(ecoli) |> Taxon
Taxon("s__Escherichia_coli", missing)

julia> String(ecoli) |> taxon
Taxon("Escherichia_coli", :species)
```

### GeneFunction

The [`GeneFunction`](@ref) type contains a name and (optionally) a [`Taxon`](@ref).
In addition to providing both a name and `Taxon`,
you can instantiate a `GeneFunction` with just a name (in which case the taxon will be `missing`),
or with the name of the taxon (in which case it will not have a `rank`).

```jldoctest genefunction
julia> gf1 = GeneFunction("gene1")
GeneFunction("gene1", missing)

julia> gf2 = GeneFunction("gene2", "Species_name")
GeneFunction("gene2", Taxon("Species_name", missing))

julia> gf3 = GeneFunction("gene2", Taxon("Species_name", :species))
GeneFunction("gene2", Taxon("Species_name", :species))
```

You can access or check for various fields using similar methods as for `Taxon`:

```jldoctest genefunction
julia> hastaxon(gf1)
false

julia> hastaxon(gf2)
true

julia> hasrank(gf2)
false

julia> hasrank(gf3)
true

julia> name(gf3)
"gene2"

julia> taxon(gf3)
Taxon("Species_name", :species)

julia> taxrank(gf3)
:species
```

For compatibility with other tools,
Converting a `GeneFunction` to a `String` if it has a `Taxon`
will include the taxon name separated by `|`.
Converting back can be done using [`genefunction`](@ref)
(note the lowercase g and f).

```jldoctest genefunction
julia> String(gf3)
"gene2|s__Species_name"

julia> genefunction(String(gf3))
GeneFunction("gene2", Taxon("Species_name", :species))
```


### Metabolites

The [`Metabolite`](@ref) type has a `name` and optionally
a `commonname`, a mass / charge ratio (`mz`), and retention time (`rt`).


```jldoctest
julia> m = Metabolite("name", "common", 1., 2.)
Metabolite("name", "common", 1.0, 2.0)

julia> name(m)
"name"

julia> commonname(m)
"common"

julia> masscharge(m)
1.0

julia> retentiontime(m)
2.0

julia> m2 = Metabolite("other name")
Metabolite("other name", missing, missing, missing)
```

## Types and Methods

```@docs
MicrobiomeSample
metadata
```

```@docs
Taxon
name
hasrank
taxrank
taxon
```

```@docs
GeneFunction
genefunction
```

```@docs
Metabolite
commonname
masscharge
retentiontime
``````@meta
CurrentModule = Microbiome
DocTestSetup  = quote
    using Microbiome
    using SparseArrays
end
```

# Diversity measures on communities and samples

## Alpha Diversity

```@docs
shannon
shannon!
ginisimpson
ginisimpson!
```

## Beta Diversity

Quite often, it's useful to boil stuff down to distances between samples.
`AbstractAbundanceTable`s take advantage of the `pairwise` function
from [`Distances.jl`](https://github.com/JuliaStats/Distances.jl)
to get a symetric distance matrix.

Right now, only Bray-Curtis, Jaccard, and Hellinger are implemented, 
but it would be straightforward to add any others.
Open an issue if you want them!

You can also get fit a principal coordinates analysis (PCoA) to your `AbstractAbundanceTable`
using the `fit(MDS, ...)` from `MultivariateStats.jl` under the hood.

```@docs
braycurtis
jaccard
hellinger
pcoa
```
