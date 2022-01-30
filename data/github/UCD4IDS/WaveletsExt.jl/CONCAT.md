# WaveletsExt.jl

| Docs | Build | Test |
|------|-------|------|
| [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://UCD4IDS.github.io/WaveletsExt.jl/stable) | [![CI](https://github.com/UCD4IDS/WaveletsExt.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/UCD4IDS/WaveletsExt.jl/actions) | [![codecov](https://codecov.io/gh/UCD4IDS/WaveletsExt.jl/branch/master/graph/badge.svg?token=U3EOscAvPE)](https://codecov.io/gh/UCD4IDS/WaveletsExt.jl) |
| [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://UCD4IDS.github.io/WaveletsExt.jl/dev) | | |

[![status](https://joss.theoj.org/papers/af5f6558736b9c3ec2bd3cf36b0cdf40/status.svg)](https://joss.theoj.org/papers/af5f6558736b9c3ec2bd3cf36b0cdf40)
[![deps](https://juliahub.com/docs/WaveletsExt/deps.svg)](https://juliahub.com/ui/Packages/WaveletsExt/iZ29j?t=2)
[![version](https://juliahub.com/docs/WaveletsExt/version.svg)](https://juliahub.com/ui/Packages/WaveletsExt/iZ29j)
[![pkgeval](https://juliahub.com/docs/WaveletsExt/pkgeval.svg)](https://juliahub.com/ui/Packages/WaveletsExt/iZ29j)

This package is a [Julia](https://github.com/JuliaLang/julia) extension package to
[Wavelets.jl](https://github.com/JuliaDSP/Wavelets.jl) (WaveletsExt is short for Wavelets
Extension). It contains additional functionalities that complement Wavelets.jl, namely
- Multi-dimensional wavelet transforms
- Redundant wavelet transforms
    - [Autocorrelation Wavelet Transforms (ACWT)](https://www.spiedigitallibrary.org/conference-proceedings-of-spie/1826/1/Wavelets-their-autocorrelation-functions-and-multiresolution-representations-of-signals/10.1117/12.131585.short)
    - [Stationary Wavelet Transforms (SWT)](https://doi.org/10.1007/978-1-4612-2544-7_17)
    - [Shift Invariant Wavelet Transforms (SIWT)](https://doi.org/10.1016/S0165-1684(97)00007-8)
- Best basis algorithms
    - [Joint best basis (JBB)](https://ieeexplore.ieee.org/document/119732)
    - [Least statistically dependent basis (LSDB)](https://ieeexplore.ieee.org/document/750958)
- Denoising methods
    - [Relative Error Shrink (RelErrorShrink)](https://ieeexplore.ieee.org/document/7752982)
    - [Stein Unbiased Risk Estimator Shrink (SUREShrink)](https://www.tandfonline.com/doi/abs/10.1080/01621459.1995.10476626)
- Wavelet transform based feature extraction techniques
    - [Local Discriminant Basis (LDB)](https://www.spiedigitallibrary.org/conference-proceedings-of-spie/2303/1/Local-discriminant-bases/10.1117/12.188763.short)

## Authors
This package is written and maintained by Zeng Fung Liew and Shozen Dan under the supervision of Professor Naoki Saito at the University of California, Davis.

## What's New (v0.1.13)
- **Changes in supported types in `denoise` and `denoiseall` functions.** For the `inputtype` positional argument, the initially supported arguments `:acwt` and `:acwpt` are now changed to `:acdwt` and `:acwpd` to match the function name change in `WaveletsExt.ACWT`.
- **2D Local Discriminant Basis now supported.** 2D version of LDB is now up and running without any changes in the syntax compared to the 1D version.

## What's New (v0.1.12)
- **Bug fixes on best basis algorithms** to allow compatibility when partial wavelet decomposition is run.
- **New function `plot_tfbdry2()` implemented.** Visual representation of leaf nodes for 2D best basis trees now available.

## Installation
The package is part of the official Julia Registry. It can be install via the Julia REPL.
```julia
(@1.7) pkg> add WaveletsExt
```
or
```julia
julia> using Pkg; Pkg.add("WaveletsExt")
```
## Usage
Load the WaveletsExt module along with Wavelets.jl.
```julia
using Wavelets, WaveletsExt
```

## References
[1] Coifman, R.R., Wickerhauser, M.V. (1992). *Entropy-based algorithms for best basis
selection*. DOI: [10.1109/18.119732](https://ieeexplore.ieee.org/document/119732) <br>
[2] Saito, N. (1998). *The least statistically-dependent basis and its applications*. DOI:
[10.1109/ACSSC.1998.750958](https://ieeexplore.ieee.org/document/750958) <br>
[3] Beylkin, G., Saito, N. (1992). *Wavelets, their autocorrelation functions, and
multiresolution representations of signals*. DOI:
[10.1117/12.131585](https://www.spiedigitallibrary.org/conference-proceedings-of-spie/1826/1/Wavelets-their-autocorrelation-functions-and-multiresolution-representations-of-signals/10.1117/12.131585.short)
<br>
[4] Nason, G.P., Silverman, B.W. (1995) *The Stationary Wavelet Transform and some
Statistical Applications*. DOI:
[10.1007/978-1-4612-2544-7_17](https://doi.org/10.1007/978-1-4612-2544-7_17) <br>
[5] Donoho, D.L., Johnstone, I.M. (1995). *Adapting to Unknown Smoothness via Wavelet
Shrinkage*. DOI:
[10.1080/01621459.1995.10476626](https://www.tandfonline.com/doi/abs/10.1080/01621459.1995.10476626)
<br>
[6] Saito, N., Coifman, R.R. (1994). *Local Discriminant Basis*. DOI:
[10.1117/12.188763](https://www.spiedigitallibrary.org/conference-proceedings-of-spie/2303/1/Local-discriminant-bases/10.1117/12.188763.short)
<br>
[7] Saito, N., Coifman, R.R. (1995). *Local discriminant basis and their applications*. DOI:
[10.1007/BF01250288](https://doi.org/10.1007/BF01250288) <br>
[8] Saito, N., Marchand, B. (2012). *Earth Mover's Distance-Based Local Discriminant Basis*.
DOI: [10.1007/978-1-4614-4145-8_12](https://doi.org/10.1007/978-1-4614-4145-8_12) <br>
[9] Cohen, I., Raz, S., Malah, D. (1997). *Orthonormal shift-invariant wavelet packet
decomposition and representation*. DOI:
[10.1016/S0165-1684(97)00007-8](https://doi.org/10.1016/S0165-1684(97)00007-8) <br>
[10] Irion, J., Saito, N. (2017). *Efficient Approximation and Denoising of Graph Signals
Using the Multiscale Basis Dictionaries*. DOI: [10.1109/TSIPN.2016.2632039](https://ieeexplore.ieee.org/document/7752982)


## TODO(long term):
* Inverse Transforms for Shift-Invariant WPT
* nD wavelet transforms for redundant and non-redundant versions# Contributing Code

## How to contribute
The preferred way to contribute to WaveletsExt.jl is to fork the [main repository](https://github.com/UCD4IDS/WaveletsExt.jl) on Github.

1. *Search or open an [issue](https://github.com/UCD4IDS/WaveletsExt.jl/issues) for this project*:   
    * Search for and click into a related issue that you would like to contribute to. If the issue has not been assigned to anyone, leave a comment and assign the issue to yourself by clicking on the 'Assignees' button on the right.

    * If no such issues have been raised, open a new issue by clicking on the 'New issue' on the top right corner. Write a brief description of the issue and wait for the feedbacks from core developers before moving on to the next steps.

1. *Fork the [project repository](https://github.com/UCD4IDS/WaveletsExt.jl)*: Click on the 'Fork' button on the top right corner of the page. This creates a copy of the code under your account on the GitHub server. 

1. *Clone this copy to your local disk*: Open up a terminal on your computer, navigate to your preferred directory, and copy the following.
```
$ git clone git@github.com:<YourGithubUsername>/WaveletsExt.jl.git
$ cd WaveletsExt
```

4. *Instantiate the project*: Open up your Julia REPL, and instantiate the current project in the package manager.
```
$ julia
```

```julia
julia> ]
(v1.7) pkg> activate .
(WaveletsExt) pkg> instantiate
```

5. *Create a branch to hold your changes*:
```
git checkout -b my-feature
```

6. Work on this copy on your computer using your preferred code editor such as VSCode. Make sure you add the [corresponding tests](https://docs.julialang.org/en/v1/stdlib/Test/) in the `test/` directory where appropriate. Use Git to do the version control. When you're done editing, do the following to record your changes in Git:
```
$ git add modified_files
$ git commit
```

7. *Push your changes* to Github with:
```
$ git push -u origin my-feature
```

8. Finally, go to the web page of your fork of the WaveletsExt.jl repository, and *click 'Pull request'* to send your changes to the maintainers for review. This will send an email to the committers.

(If any of the above seems like magic to you, then look up the [Git documentation](https://git-scm.com/doc) on the web.)

### Feature contribution notes
It is recommended to check that your contribution complies with the following rules before submitting a pull request:

- All public methods should have informative docstrings with sample usage presented.

You can also check for common programming errors by testing your code in the package manager mode:
```julia
(WaveletsExt) pkg> test
```

## Filing bugs and feature requests
We use Github issues to track all bugs and feature requests; feel free to [open an issue](https://github.com/UCD4IDS/WaveletsExt.jl/issues) if you have found a bug or wish to see a feature implemented.

It is recommended to check that your issue complies with the following rules before submitting:

* Verify that your issue is not being currently addressed by other [issues](https://github.com/UCD4IDS/WaveletsExt.jl/issues) or [pull requests](https://github.com/UCD4IDS/WaveletsExt.jl/pulls).

* Please ensure all code snippets and error messages are formatted in appropriate code blocks. See [Creating and highlighting code blocks](https://docs.github.com/en/github/writing-on-github/working-with-advanced-formatting/creating-and-highlighting-code-blocks).

* Please include your Julia and WaveletsExt.jl version. This information can be found by running the following code snippet:
```julia
import Pkg
println("Julia $VERSION")       # Julia version
Pkg.status("WaveletsExt")       # Package version
```

## Documentation
You can edit the documentation using any text editor and then generate the HTML output by doing:
```
$ julia --project=docs/ docs/make.jl
```
The resulting HTML files will be placed in `docs/build/` and are viewable in a web browser. See the documentation from [Documenter.jl](https://juliadocs.github.io/Documenter.jl/stable/man/guide/) for more information.

*Note: Do not commit any files from the `docs/build/` directory. These webpages will be generated automatically via Github Actions when the commit is merged with the master branch of the project repository.*

*Tip: Generally, only the `docs/make.jl` and the files in `docs/src/` need to be updated.*

## Note
This document was gleefully borrowed from [librosa](https://librosa.org/doc/latest/index.html) and [scikit-learn](https://scikit-learn.org/stable/).---
title: 'WaveletsExt.jl: Extending the boundaries of wavelets in Julia'
tags:
  - Julia
  - Wavelets
authors:
  - name: Zeng Fung Liew^[co-first author] 
    affiliation: "1" 
  - name: Shozen Dan^[co-first author, corresponding author]
    affiliation: "2"
  - name: Naoki Saito
    affiliation: "3"
affiliations:
 - name: Department of Statistics, University of California, Davis, United States
   index: 1
 - name: Department of Mathematics, Imperial College London, United Kingdom
   index: 2
 - name: Department of Mathematics, University of California, Davis, United States
   index: 3
date: 22 September 2021
bibliography: paper.bib
---

# Summary
Whether it is seismic surveys, ECG signals, stock market trends, or sensor data, the wavelet and wavelet packet transforms are powerful tools for signal analysis and classification with many advantages over the conventional Fourier methods. Primary among them is the ability to extract information localized in both time and frequency domains, enabling multiresolution analysis [@Daubechies:1992; @Mallat:2009]. As such, wavelets and wavelet packets have become popular tools for computational harmonic analysis. `WaveletsExt.jl` was developed to augment `Wavelets.jl` (the existing wavelet toolbox for Julia) by providing routines for wavelet analysis, wavelet packet analysis, and associated utilities.

# Statement of Need
Julia's principal package for wavelets is `Wavelets.jl` [@JuliaDSP:2021], which provides the essential building blocks for data analysis using wavelets. These include 1-D, 2-D, and 3-D wavelet transforms via filter banks or lifting, a range of thresholding functions, and other utilities. However, as a general-purpose package for wavelets, `Wavelets.jl` does not include many targeted and sophisticated methods present in the literature.

`WaveletsExt.jl` (Wavelets Extension) enlarges the wavelet toolbox for Julia by providing a host of useful wavelet-based functions such as Stationary Wavelet Transform [@Nason:1995], Autocorrelation Wavelet Transform [@Saito:1993], Local Discriminant Basis [@Saito:1995], and Shift-invariant Wavelet Packet Decomposition [@Cohen:1995]. The package also contains denoising utilities such as SureShrink [@Donoho:1995a] and Relative Error Shrink [@Irion:2017] as well as several data visualization utilities.

One of the most distinguishing features of `WaveletsExt.jl` is the presence of algorithms for handling an ensemble of input signals. Currently, `Wavelets.jl` implements best basis selection utilities for wavelet packets for a single input. However, it does not include methods for selecting a single best basis for a set of inputs with similar properties (e.g., signals or images belonging to the same class), which is valuable for feature extraction and data compression. To address this, `WaveletsExt.jl` implements the Joint Best Basis (JBB) [@Wickerhauser:1996] and the Least Statistically Dependent Basis (LSDB) [@Saito:2001], which provide approximations of the Principal Component Analysis (PCA) and Independent Component Analysis (ICA), respectively, in a computationally fast manner through a dictionary of orthonormal bases.

# Examples
## 1. Redundant Wavelet Transforms
`WaveletsExt.jl` implements several redundant wavelet transforms including Autocorrelation Wavelet Transform [@Saito:1993] and Stationary Wavelet Transform (SWT) [@Nason:1995]. These transformations can be performed using the `acdwt` and `sdwt` functions, and the resulting decomposition can be visualized with the `wiggle` function included in `WaveletsExt.jl`.

```julia
using Plots, Wavelets, WaveletsExt

x = zeros(1<<8) # Generate a unit impulse (dirac delta) signal
x[128] = 1
wt = wavelet(WT.db4)  # Construct Daubechies 4-tap wavelet filter

# ----- Autocorrelation Wavelet Transforms -----
y = acdwt(x, wt)
p1 = wiggle(y) |> p -> plot!(p, yticks=1:9, title="Autocorrelation WT")

# ----- Stationary Wavelet Transforms -----
y = sdwt(x, wt)
p2 = wiggle(y) |> p -> plot!(p, yticks=1:9, title="Stationary WT")

# Combine and save plot
p = plot(p1, p2, layout=(1,2), size=(600,300))
savefig(p, "transforms.png")
```
!["Wiggle" plots displaying the value of coefficients at each level of the autocorrelation and stationary wavelet transform for a unit impulse signal. \label{fig:transforms}](transforms.png)

## 2. Best Basis Algorithms
`WaveletsExt.jl` can select a best basis for a multiple signal input (i.e., an array of signals) through the Joint Best Basis (JBB) [@Wickerhauser:1996] or Least Statistically Dependent Basis (LSDB) [@Saito:2001] algorithms. The resulting best basis tree can be visualized using `plot_tfbdry` also included in `WaveletsExt.jl`.

```julia
using Plots, Wavelets, WaveletsExt

# Generate 100 noisy heavysine signals of length 2⁸
x = generatesignals(:heavysine, 8) |> 
    x -> duplicatesignals(x, 100, 2, true, 0.5)

# Wavelet packet decomposition of all signals
xw = wpdall(x, wt, 6)

# ----- Joint Best Basis (JBB)
tree = bestbasistree(xw, JBB())
p1 = plot_tfbdry(tree, 
                 node_color=:green, 
                 line_color=:black, 
                 background_color=:white) |> 
     p -> plot!(p, title="JBB")

# ----- Least Statistically Dependent Basis (LSDB)
tree = bestbasistree(xw, LSDB())
p2 = plot_tfbdry(tree, 
                 node_color=:green, 
                 line_color=:black, 
                 background_color=:white) |> 
     p -> plot!(p, title="LSDB")

# Combine and save plot
p = plot(p1, p2, layout=(1,2), size=(600,300))
savefig(p, "bestbasis.png")
```
![The best basis trees of 100 HeaviSine signals (A sinusoid + two Heaviside step functions) [@Donoho:1995a; @Donoho:1995b] selected by the JBB and LSDB algorithms. Each row represents a decomposition level, where level 0 is the original input signal, and each cell represents a frequency subband (low to high frequency from left to right). The colored cells indicate those subbands selected by the JBB (left) and the LSDB (right) algorithms.  \label{fig:bestbasis}](bestbasis.png)

## 3. Denoising Algorithms
`WaveletsExt.jl` contains two functions for denoising: `denoise` and `denoiseall`. The former denoises a single signal input whereas the latter denoises multiple signal input. For more examples of denoising algorithms in `WaveletsExt.jl`, see [@Liew:2021].

```julia
using Plots, Wavelets, WaveletsExt

# Generate 6 circularly shifted HeaviSine signals
x₀ = generatesignals(:heavisine, 8) |> 
     x -> duplicatesignals(x, 6, 2, false)
     
# Generate 6 noisy versions of the original signals
x = generatesignals(:heavisine, 8) |> 
    x -> duplicatesignals(x, 6, 2, true, 0.8)

# Decompose each noisy signal
xw = wpdall(x, wt)

# Get best basis tree from the decomposition of signals
bt = bestbasistree(xw, JBB())

# Get best basis coefficients based on best basis tree
y = bestbasiscoef(xw, bt)

# Denoise all signals based on computed best basis tree
x̂ = denoiseall(y, :wpt, wt, tree=bt)

# Plot results
p1 = plot(title="Noisy Signals")
wiggle!(x₀, sc=0.7, FaceColor=:white, ZDir=:reverse)
wiggle!(x, sc=0.7, FaceColor=:white, ZDir=:reverse)

p2 = plot(title="Denoised Signals")
wiggle!(x₀, sc=0.7, FaceColor=:white, ZDir=:reverse)
wiggle!(x̂, sc=0.7, FaceColor=:white, ZDir=:reverse)

# Combine and save plot
p = plot(p1, p2, layout=(1,2), size=(600,300))
savefig(p, "denoising.png")
```
![Left: HeaviSine signals with Gaussian noise. Black lines represent the original (non-noisy) signal. Right: Simultaneously denoised signals using the JBB algorithm with a universal thresholding constant determined by the VisuShrink method [@Donoho:1994]. \label{fig:denoising}](denoising.png)

## 4. Feature Extraction
For signal classification problems, users can extract distinguishing features localized in both the time and frequency domains using the Local Discriminant Basis (LDB) algorithm. Further details can be found in the original papers by Saito and his collaborators [@Saito:1995; @Saito:2002] as well as the interactive tutorial [@Dan:2021].

```julia
using Plots, Wavelets, WaveletsExt

# Generate 100 signals for each class of cylinder-bell-funnel
X, y = generateclassdata(ClassData(:cbf, 100, 100, 100))
# View sample signals and how each class differs from one another
cylinder = wiggle(X[:,1:5], sc=0.3, EdgeColor=:white, FaceColor=:white)
plot!(cylinder, yticks=1:5, ylabel="Cylinder")
bell = wiggle(X[:,101:105], sc=0.3, EdgeColor=:white, FaceColor=:white)
plot!(bell, yticks=1:5, ylabel="Bell")
funnel = wiggle(X[:,201:205], sc=0.3, EdgeColor=:white, FaceColor=:white)
plot!(funnel, yticks=1:5, ylabel="Funnel")
p1 = plot(cylinder, bell, funnel, layout=(3,1))

# Instantiate the LDB object
wt = wavelet(WT.coif4)
ldb = LocalDiscriminantBasis(
  wt=wt,
  max_dec_level=6,
  dm=SymmetricRelativeEntropy(),
  en=TimeFrequency(),
  dp=BasisDiscriminantMeasure(),
  top_k=10,
  n_features=10 # Number of features to extract
)
                            
# Extract features using LDB
X̂ = fit_transform(ldb, X, y)

# Plot the best basis for feature extraction
p2 = plot_tfbdry(ldb.tree, 
                 node_color=:green, 
                 line_color=:black, 
                 background_color=:white)
plot!(p2, title="Basis Selection using LDB")

# Combine and save plot
p = plot(p1, p2, size=(600,300))
savefig(p, "ldb.png")
```
![Left: Examples of Cylinder, Bell, and Funnel signals. Right: The best basis tree selected by the LDB algorithm for discriminating the three classes of signals. \label{fig:denoising}](ldb.png)

# Reproducible Research
`WaveletsExt.jl` was partially inspired by the `WaveLab` library in MATLAB, which was developed to enable reproducible wavelet research [@Donoho:1995b]. In this spirit, we wrote a series of tutorials, examples, and experiments using `Pluto.jl` [@Liew:2021; @Dan:2021], a platform with which Julia users can create and share reactive documents [@Fonsp:2021]. By downloading and running these so-called Pluto notebooks, researchers and students alike can reproduce the results of our research and interactively adjust parameters to see the changes in experiment outcomes.

# Acknowledgements
This project was partially supported by the following grants from the US National Science Foundation: DMS-1148643; DMS-1418779; DMS-1912747; and CCF-1934568.

# References
# JOSS paper

This directory contains the materials for the paper submitted to the Journal of Open Source Software (JOSS).
```@raw html
<div style="width:100%; height:150px;border-width:4px;border-style:solid;padding-top:25px;
        border-color:#000;border-radius:10px;text-align:center;background-color:#B3D8FF;
        color:#000">
    <h3 style="color: black;">Star us on GitHub!</h3>
    <a class="github-button" href="https://github.com/UCD4IDS/WaveletsExt.jl.git" data-icon="octicon-star" data-size="large" data-show-count="true" aria-label="Star UCD4IDS/WaveletsExt.jl on GitHub" style="margin:auto">Star</a>
    <script async defer src="https://buttons.github.io/buttons.js"></script>
</div>
```

# WaveletsExt.jl
This package is a [Julia](https://github.com/JuliaLang/julia) extension package to
[Wavelets.jl](https://github.com/JuliaDSP/Wavelets.jl) (WaveletsExt is short for Wavelets
Extension). It contains additional functionalities that complement Wavelets.jl, namely
- Multi-dimensional wavelet transforms
- Redundant wavelet transforms
    - [Autocorrelation Wavelet Transforms (ACWT)](https://www.spiedigitallibrary.org/conference-proceedings-of-spie/1826/1/Wavelets-their-autocorrelation-functions-and-multiresolution-representations-of-signals/10.1117/12.131585.short)
    - [Stationary Wavelet Transforms (SWT)](https://doi.org/10.1007/978-1-4612-2544-7_17)
    - [Shift Invariant Wavelet Transforms (SIWT)](https://doi.org/10.1016/S0165-1684(97)00007-8)
- Best basis algorithms
    - [Joint best basis (JBB)](https://ieeexplore.ieee.org/document/119732)
    - [Least statistically dependent basis (LSDB)](https://ieeexplore.ieee.org/document/750958)
- Denoising methods
    - [Relative Error Shrink (RelErrorShrink)](https://ieeexplore.ieee.org/document/7752982)
    - [Stein Unbiased Risk Estimator Shrink (SUREShrink)](https://www.tandfonline.com/doi/abs/10.1080/01621459.1995.10476626)
- Wavelet transform based feature extraction techniques
    - [Local Discriminant Basis (LDB)](https://www.spiedigitallibrary.org/conference-proceedings-of-spie/2303/1/Local-discriminant-bases/10.1117/12.188763.short).

## Authors
This package is written and maintained by Zeng Fung Liew and Shozen Dan under the supervision of Professor Naoki Saito at the University of California, Davis.

## Installation
The package is part of the official Julia Registry. It can be install via the Julia REPL.
```julia
(v1.6) pkg> add WaveletsExt
```
or
```julia
julia> using Pkg; Pkg.add("WaveletsExt")
```
## Usage
Load the WaveletsExt module along with Wavelets.jl.
```julia
using Wavelets, WaveletsExt
```

## References
[1] Coifman, R.R., Wickerhauser, M.V. (1992). *Entropy-based algorithms for best basis
selection*. DOI: [10.1109/18.119732](https://ieeexplore.ieee.org/document/119732)

[2] Saito, N. (1998). *The least statistically-dependent basis and its applications*. DOI:
[10.1109/ACSSC.1998.750958](https://ieeexplore.ieee.org/document/750958)

[3] Beylkin, G., Saito, N. (1992). *Wavelets, their autocorrelation functions, and
multiresolution representations of signals*. DOI:
[10.1117/12.131585](https://www.spiedigitallibrary.org/conference-proceedings-of-spie/1826/1/Wavelets-their-autocorrelation-functions-and-multiresolution-representations-of-signals/10.1117/12.131585.short)

[4] Nason, G.P., Silverman, B.W. (1995) *The Stationary Wavelet Transform and some
Statistical Applications*. DOI:
[10.1007/978-1-4612-2544-7_17](https://doi.org/10.1007/978-1-4612-2544-7_17) 

[5] Donoho, D.L., Johnstone, I.M. (1995). *Adapting to Unknown Smoothness via Wavelet
Shrinkage*. DOI:
[10.1080/01621459.1995.10476626](https://www.tandfonline.com/doi/abs/10.1080/01621459.1995.10476626) 

[6] Saito, N., Coifman, R.R. (1994). *Local Discriminant Basis*. DOI:
[10.1117/12.188763](https://www.spiedigitallibrary.org/conference-proceedings-of-spie/2303/1/Local-discriminant-bases/10.1117/12.188763.short)

[7] Saito, N., Coifman, R.R. (1995). *Local discriminant basis and their applications*. DOI:
[10.1007/BF01250288](https://doi.org/10.1007/BF01250288)

[8] Saito, N., Marchand, B. (2012). *Earth Mover's Distance-Based Local Discriminant Basis*.
DOI: [10.1007/978-1-4614-4145-8_12](https://doi.org/10.1007/978-1-4614-4145-8_12)

[9] Cohen, I., Raz, S., Malah, D. (1997). *Orthonormal shift-invariant wavelet packet
decomposition and representation*. DOI:
[10.1016/S0165-1684(97)00007-8](https://doi.org/10.1016/S0165-1684(97)00007-8)

[10] Irion, J., Saito, N. (2017). *Efficient Approximation and Denoising of Graph Signals
Using the Multiscale Basis Dictionaries*. DOI: [10.1109/TSIPN.2016.2632039](https://ieeexplore.ieee.org/document/7752982)
# [Extracting the best bases from signals](@id bestbasis_manual)
The Wavelets.jl's package contains a best basis algorithm (via [`bestbasistree`](@ref
WaveletsExt.BestBasis.bestbasistree)) that search for the basis tree within a single signal
$x$ such that the Shannon's Entropy or the Log Energy Entropy of the basis of the signal is
minimized, ie. ``\min_{b \in B} M_x(b)``, where  
- ``B`` is the collection of bases for the signal ``x``,
- ``b`` is a basis for the signal ``x``, and
- ``M_x(.)`` is the information cost (eg. Shannon's entropy or Log Energy entropy) for
the basis of signal ``x``. 

However, the challenge arises when there is a need to work on a group of signals $X$ and
find a single best basis $b$ that minimizes the information cost $M_x$ for all $x \in X$, ie. ``\min_{b \in B} \sum_{x \in X} M_x(b)``

A brute force search is not ideal as its computational time grows exponentially to the
number and size of the signals. Here, we have the two efficient
algorithms for the estimation of an overall best basis:
- Joint Best Basis (JBB),
- Least Statistically Dependent Basis (LSDB)

## Best basis representations
To represent a best basis tree in the most memory efficient way, Wavelets.jl uses Julia's
`BitVector` data structure to represent a binary tree that corresponds to the bases of 1D
signals. The indices of the `BitVector`s correspond to nodes of the trees as follows:

|Figure 1: Binary tree structure | Figure 2: Binary tree indexing |
|:---:|:---:|
|![](../fig/binary_tree.PNG) |![](../fig/binary_tree_indexing.PNG) |

where L corresponds to a low pass filter and H corresponds to a high pass filter.

Similarly, a quadtree, used to represent the basis of a 2D wavelet transform, uses a
`BitVector` with the following indexing:

|Figure 3: Quadtree structure | Figure 4: Quadtree indexing |
|:---:|:---:|
|![](../fig/quad_tree.PNG) |![](../fig/quad_tree_indexing.PNG) |

## Examples
### 1D Best Basis
Assume we are given a large amount of signals to transform to its best basis, one may use
the following approach, which uses the standard best basis algorithm to search for the best basis of each signal.

```@example wt
using Wavelets, WaveletsExt, Plots

# Generate 4 HeaviSine signals
x = generatesignals(:heavisine, 7)
X = duplicatesignals(x, 4, 2, true, 0.5)

# Construct wavelet
wt = wavelet(WT.haar)

# Decomposition of all signals
xw = wpdall(X, wt)

# Best basis trees, each column corresponds to 1 tree
trees = bestbasistreeall(xw, BB()); 
nothing # hide
```

Similarly, we can also find the best basis of the signals using Joint Best Basis (JBB) and Least Statistically Dependent Basis (LSDB). Note that these algorithms return 1 basis tree that generalizes the best basis over the entire set of signals.
```@example wt
# JBB
treeⱼ = bestbasistree(xw, JBB())

# LSDB
treeₗ = bestbasistree(xw, LSDB());
nothing # hide
```

One can then view the selected nodes from the best basis trees using the [`plot_tfbdry`](@ref WaveletsExt.Visualizations.plot_tfbdry) function as shown below.

```@example wt
# Plot each tree
p1 = plot_tfbdry(trees[:,1])
plot!(p1, title="Signal 1 best basis")
p2 = plot_tfbdry(trees[:,2])
plot!(p2, title="Signal 2 best basis")
p3 = plot_tfbdry(trees[:,3])
plot!(p3, title="Signal 3 best basis")
p4 = plot_tfbdry(trees[:,4])
plot!(p4, title="Signal 4 best basis")
p5 = plot_tfbdry(treeⱼ)
plot!(p5, title="JBB")
p6 = plot_tfbdry(treeₗ)
plot!(p6, title="LSDB")

# Draw all plots
plot(p1, p2, p3, p4, p5, p6, layout=(3,2))
```

### 2D Best Basis
Similar to the mechanics of the best basis algorithms for 1D signals, WaveletsExt also supports the best basis search for 2D signals.
```@example wt
using Images, TestImages

img = testimage("moonsurface");
y = convert(Array{Float64}, img);
Y = duplicatesignals(y, 4, 2, true, 0.5);

# Decomposition of all signals
yw = wpdall(Y, wt, 5)

# Best basis trees
trees = bestbasistreeall(yw, BB())
treeⱼ = bestbasistree(yw, JBB())
treeₗ = bestbasistree(yw, LSDB())

# Plot each tree
p1 = plot_tfbdry2(trees[:,1])
plot!(p1, title="Signal 1 best basis")
p2 = plot_tfbdry2(trees[:,2])
plot!(p2, title="Signal 2 best basis")
p3 = plot_tfbdry2(trees[:,3])
plot!(p3, title="Signal 3 best basis")
p4 = plot_tfbdry2(trees[:,4])
plot!(p4, title="Signal 4 best basis")
p5 = plot_tfbdry2(treeⱼ)
plot!(p5, title="JBB")
p6 = plot_tfbdry2(treeₗ)
plot!(p6, title="LSDB")

# Draw all plots
plot(p1, p2, p3, p4, p5, p6, layout=(3,2))
```

## [Best Basis of Shift-Invariant Wavelet Packet Decomposition](@id si_bestbasis)
One can think of searching for the best basis of the shift-invariant wavelet packet decomposition as a problem of finding ``\min_{b \in B} \sum_{x \in X} M_x(b)``, where ``X`` is all the possible shifted versions of an original signal ``y``. One can compute the best basis tree as follows:
```@example wt
xw = siwpd(x, wt)

# SIBB
tree = bestbasistree(xw, 7, SIBB());
nothing #hide
```

!!! warning 
    SIWPD is still undergoing large changes in terms of data structures and efficiency improvements. Syntax changes may occur in the next patch updates.
# [Signal Denoising](@id denoising_manual)
Wavelet denoising is an important step in signal analysis as it helps remove unnecessary high frequency noise while maintaining the most important features of the signal. Intuitively, signal denoising comes in the following simple steps:
1. Decompose a signal or a group of signals. One can choose to decompose signals into its best basis tree for more optimal results.
2. Find a suitable threshold value. There are many ways to do so, with VisuShrink (D. Donoho, I. Johnstone) being one of the most popular approaches. The VisuShrink implementation in Wavelets.jl, along with the RelErrorShrink and the SureShrink implementations in WaveletsExt.jl give users more threshold selection options.
3. Threshold the wavelet coefficients. There are various thresholding methods implemented in Wavelets.jl for this purpose, with Hard and Soft thresholding being the usual go-to method due to its simplistic approach.
4. Reconstruct the original signals using the thresholded coefficients.

For more information and examples on wavelet denoising using WaveletsExt.jl, visit [Wavelets Denoising Experiment repository under UCD4IDS](https://github.com/UCD4IDS/WaveletsDenoisingExperiment) for a step-by-step tutorial in a Pluto notebook. The following is a simple guide on denoisng using WaveletsExt.jl.

## Denoising a single signal
To denoise a single signal, one can use the `denoise` function from WaveletsExt.jl as shown below. Note the following key parameters:
- `x`: Input signal.
- `inputtype`: Type of input. One can input an original signal `:sig`, or first transform the signal and type in one of `:dwt`, `:wpt`, `:sdwt`, `:swpd`, `:acdwt`, and `:acwpd`.
- `wt`: Transform wavelet.
- `L`: Number of decomposition levels. Necessary for input types `:sig`, `:dwt`, and `:sdwt`.
- `tree`: Decomposition tree of the signals. Necessary for input types `:wpt` and `:swpd`.
- `dnt`: Denoise type. One should input either of VisuShrink, RelErrorShrink, or SureShrink.
- `estnoise`: Noise estimation. Can be a function or a value.
- `smooth`: Smoothing method used.
For more detailed information, visit the [denoising API](@ref denoising_api) page.

```@example
using Wavelets, WaveletsExt, Random, Plots

# define function and wavelet
x₀ = generatesignals(:heavisine, 8)
x = x₀ + 0.8*randn(256)
wt = wavelet(WT.db4)

# best basis tree
xw = wpd(x, wt)
bt = bestbasistree(xw, BB())
y = getbasiscoef(xw, bt)

# denoise
x̂ = denoise(y, :wpt, wt, tree=bt)

# plot results
nothing # hide
plot([x₀ x x̂], title="Denoising Example", label=["original" "noisy" "denoised"],
     lw=[3 1 2], lc=[:black :grey :red])
```

## Denoising a group of signals
Similar to the `denoise` function we saw previously, for denoising a group of signals, one can use the `denoiseall` function. The parameters used are the same, with the following addition:
- `bestTH`: Method to determine the best threshold value for a group of signals. One can choose each signal's individual best threshold value, or use a function such as `mean` or `median` to generalize an overall best threshold value.

For more detailed information, visit the [denoising API](@ref denoising_api) page.
```@example
using Wavelets, WaveletsExt, Random, Plots

# define function and wavelet
x = generatesignals(:heavisine, 8)
X₀ = duplicatesignals(x, 6, 2, false)
X = duplicatesignals(x, 6, 2, true, 0.8)
wt = wavelet(WT.db4)

# decomposition
coef = wpdall(X, wt)

# best basis tree
bt = bestbasistree(coef, JBB())
Y = getbasiscoefall(coef, bt)

# denoise
X̂ = denoiseall(Y, :wpt, wt, tree=bt)

# plot results
nothing # hide
wiggle(X₀, sc=0.7, FaceColor=:white, ZDir=:reverse)
wiggle!(X, sc=0.7, EdgeColor=:grey, FaceColor=:white, ZDir=:reverse)
wiggle!(X̂, sc=0.7, EdgeColor=:red, FaceColor=:white, ZDir=:reverse)
plot!(title="Group Denoising Example")
```# [Local Discriminant Basis](@id ldb_manual)

Local Discriminant Basis is a feature extraction technique developed by N. Saito and R. Coifman in 1995. This algorithm follows the following basic steps:

1. Decompose a set of multi-class signals using wavelet packet decomposition. A wavelet packet decomposition decomposes a signal into multiple nodes which resembles a binary tree.
2. Based on the decomposed wavelet coefficients, build an energy map based on time-frequency or probability density.
3. Using the energy map, compute the discriminant measure and select a basis tree that best discriminates the different classes of signals.
4. Based on the selected basis tree, extract the corresponding wavelet coefficients for each signal.
5. Compute the discriminant power of each coefficient index. Select the top k set of coefficients to be used as features to be passed onto a classifier such as Linear Discriminant Analysis (LDA) and Classification and Regression Trees (CART).

A more in-depth tutorial can be found in the Pluto notebook [here](https://github.com/ShozenD/LDBExperiments). For more information on LDB, please refer to the original paper "Local Discriminant Basis and their Applications" by Saito and Coifman [here](https://www.math.ucdavis.edu/~saito/publications/saito_ldb_jmiv.pdf).

## Example
We first generate a multi-class dataset. WaveletsExt.jl has 2 built-in multi-class signals dataset, namely the triangular signals (`:tri`) and the cylinder-bell-funnel signals (`:cbf`).
```@example ldb_tutorial
using Wavelets, WaveletsExt, Plots

# generates 100 signals for each class of cylinder-bell-funnel
X, y = generateclassdata(ClassData(:cbf, 100, 100, 100));

# view sample signals and how each class differs from one another
cylinder = wiggle(X[:,1:5], sc=0.3)
plot!(cylinder, title="Cylinder signals")
bell = wiggle(X[:,101:105], sc=0.3)
plot!(bell, title="Bell signals")
funnel = wiggle(X[:,201:205], sc=0.3)
plot!(funnel, title="Funnel signals")
plot(cylinder, bell, funnel, layout=(3,1))
```

Next, we define the parameters for our Local Discriminant Basis object. Here are a few key parameters to note:
* `wt`: Type of wavelet used. Default is `wavelet(WT.haar)`.
* `max_dec_level`: Maximum decomposition level. Default is to decompose each signal all the way to its maximum possible depth.
* `dm`: Type of discriminant measure. Available choices are:
    - `AsymmetricRelativeEntropy()` (default)
    - `SymmetricRelativeEntropy()`
    - `LpEntropy()`
    - `HellingerDistance()`
* `en`: Type of energy map. Available choices are:
    - `TimeFrequency()` (default)
    - `ProbabilityDensity()`
* `dp`: Type of discriminant power. Available choices are:
    - `BasisDiscriminantMeasure()` (default)
    - `FishersClassSeparability()`
    - `RobustFishersClassSeparability()`
* `top_k`: Max number of coefficients used in each node for the computation of discriminant power. The default setting uses all available coefficients for the computation.
* `n_features`: Number of features to be returned. All features/coefficients will be returned by default.
```@example ldb_tutorial
wt = wavelet(WT.coif4);
ldb = LocalDiscriminantBasis(
    wt=wt, 
    max_dec_level=7,
    dm=SymmetricRelativeEntropy(), 
    en=TimeFrequency(),
    dp=BasisDiscriminantMeasure(),
    top_k=10,
    n_features=10
);

# transform and extract the features using LDB
X̂ = fit_transform(ldb, X, y);
nothing # hide
```

After fitting our data, we will then also be able to conduct our own analysis. We can observe where the best basis is selected from using the `plot_tfbdry` function.
```@example ldb_tutorial
plot_tfbdry(ldb.tree)
```

Another thing we can do is observe the heatmap produced by the discriminant measure (`ldb.DM`).
```@example ldb_tutorial
heatmap(1:ldb.sz[1], 0:ldb.max_dec_level, ldb.DM);
plot!(title="Discriminant Measure Heatmap")
```

To decide how many features we should select, we can use the elbow rule on the discriminant powers (`ldb.DP`). From the plot below, we can see that approximately 6 features should be chosen for the classification step.
```@example ldb_tutorial
plot(ldb.DP[ldb.order], labels="discriminant power");
plot!(title="Plot of LDB Discriminant Power")
```

Knowing the 6 features we want to select, we can go one step further and examine the basis vectors generated by the coefficients of these 6 indices by defining the function below. In the illustration purpose of this tutorial, the basis vectors generated by the coefficients of the top 10 features are plotted below.
```@example ldb_tutorial
function get_basisvectors(n::Integer, wt::DiscreteWavelet, tree::BitVector,
        idx::Vector{<:Integer})

    k = length(idx)
    y = Array{Float64,2}(undef, (n,k))
    for (i,j) in enumerate(idx)
        x = zeros(n)
        x[j] = 1
        y[:,i] = iwpt(x, wt, tree)
    end
    return y
end

bases = get_basisvectors(128, ldb.wt, ldb.tree, ldb.order[1:10]);
nothing # hide
wiggle(bases, sc=0.3, ZDir=:reverse);
plot!(title="Top 10 LDB vectors")
```

Since we have decided that 6 features are optimum for classification purposes, we can use the `change_nfeatures` function as below.
```@example ldb_tutorial
X̂ = change_nfeatures(ldb, X̂, 6);
nothing # hide
```

If we are curious, we can use the `inverse_transform` function to observe how the signals look like if they're generated from these 6 features.
```@example ldb_tutorial
X̃  = inverse_transform(ldb, X̂);

# view sample signals and how each class differs from one another
nothing # hide
cylinder = wiggle(X̃[:,1:5], sc=0.3)
plot!(cylinder, title="Cylinder signals")
bell = wiggle(X̃[:,101:105], sc=0.3)
plot!(bell, title="Bell signals")
funnel = wiggle(X̃[:,201:205], sc=0.3)
plot!(funnel, title="Funnel signals")
plot(cylinder, bell, funnel, layout=(3,1))
```

With that said, we are essentially done with the LDB step, and we can move on to the model fitting step using packages such as [MLJ.jl](https://alan-turing-institute.github.io/MLJ.jl/stable/) and [MultivariateStats.jl](https://multivariatestatsjl.readthedocs.io/en/latest/).# [Wavelet Transforms](@id transforms_manual)
Wavelet transform is a feature extraction process for decomposing signals into high and low
frequency segments. Using a pair of orthonormal wavelet ``\psi \in L^2(\mathbb{R})``, where
``L^2(\mathbb{R})`` is a Hilbert space of square integrable functions, one can compute
- ``y_{low} = g(x)`` where ``x`` is the signal of interest, ``g`` is the low pass filter
  corresponding to ``\psi``, and ``y_{low}`` is the output when ``x`` passes through ``g``.
- ``y_{high} = h(x)`` where ``x`` is the signal of interest, ``h`` is the high pass filter
  corresponding to ``\psi``, and ``y_{high}`` is the output when ``x`` passes through ``h``.

The wavelet transform can be thought of as an improvement over the Fourier transform due to
its ability to preserve information in both the time and frequency domains. It has vast
applications in fields such as signal analysis and image compression.

As an extension to [Wavelets.jl](https://github.com/JuliaDSP/Wavelets.jl), WaveletsExt.jl
offers additional (redundant) wavelet transform techniques via [autocorrelation wavelet
transforms](@ref ac_transforms) (Beylkin, Saito), [stationary wavelet transforms](@ref
s_transforms) (Nason, Silverman), and [shift invariant wavelet transform](@ref
si_transforms) (Cohen et. al.).

## Wavelet Transform Methods
There are essentially 3 methods of wavelet transforms: discrete wavelet transforms, wavelet
packet transforms, and wavelet packet decomposition. The overall idea of signals being
decomposed into high and low frequency segments remain the same, but the number of levels of
decomposition for each segment may vary.
### Discrete Wavelet Transforms (DWT)
The discrete wavelet transfrom only iteratively decomposes the approximation coefficients at each level, ie. iteratively transforms the output from the low pass filter. The coefficients of the leaf nodes are returned. See Figure 1 for a visualization of the DWT transform process and output.
### Wavelet Packet Transforms (WPT)
The wavelet packet transform takes the decomposition process one step further and itereatively decomposes on both the approximation and detail coefficients, ie. both the outputs from the low pass filter and high pass filter are being iteratively decomposed. The coefficients of the leaf nodes are returned.

An extension to WPT is that one can decompose a signal based on a given tree. See Figure 1 for better visualization of the transform process.

### Wavelet Packet Decomposition (WPD)
The wavelet packet decomposition functions similarly to WPT, except that all the coefficients (regardless of whether they're at the lead node) are retained in the output. The WPD is useful for selecting the [wavelet best basis](@ref bestbasis_manual) and feature extraction algorithms such as [Local Discriminant Basis](@ref ldb_manual).

| Discrete Wavelet Transform | Wavelet Packet Transform | Wavelet Packet Decomposition|
|:---:|:---:|:---:|
| ![](../fig/dwt.PNG) | ![](../fig/wpt1.PNG) | ![](../fig/wpd.PNG) |
|| ![](../fig/wpt2.PNG) ||

Figure 1: Decomposition method for DWT, WPT, and WPD respectively. Coefficient outputs from DWT, WPT, and WPD are highlighted in red.
## Types of Wavelet Transforms and Their Examples in WaveletsExt.jl
### Regular Wavelet Transform
The standard wavelet transform (DWT and WPT) from Wavelets.jl and the WPD can be performed
as follows:
```@example wt
using Wavelets, WaveletsExt

# Define function and wavelet
x = generatesignals(:heavisine, 8)
wt = wavelet(WT.db4)

# ----- Discrete Wavelet Transform (DWT) -----
y = dwt(x, wt)      # Forward transform
z = idwt(y, wt)     # Inverse transform

# ----- Wavelet Packet Transform (WPT) -----
y = wpt(x, wt)      # Forward transform
z = iwpt(y, wt)     # Inverse transform

# ----- Wavelet Packet Decomposition (WPD) -----
y = wpd(x, wt)      # Decompose into L levels
nothing # hide
```

In the case where there are multiple signals to transform, one may opt for [`dwtall`](@ref WaveletsExt.DWT.dwtall), [`wptall`](@ref WaveletsExt.DWT.wptall), and [`wpdall`](@ref WaveletsExt.DWT.wpdall).

### [Stationary Wavelet Transforms] (@id s_transforms)
The stationary wavelet transform is a redundant type of wavelet transform. This means that there are no downsampling involved unlike the standard transforms, resulting in an exponentially larger number of coefficients compared to that of the standard transforms. A strength of the stationary wavelet transform is its ability to retain more information, thereby being more useful in certain signal analysis applications such as denoising. However, it also takes an exponentially larger amount of time and space to decompose a signal compared to the standard transforms.

```@example wt
# ----- Discrete Wavelet Transform (DWT) -----
y = sdwt(x, wt)     # Forward transform
z = isdwt(y, wt)    # Inverse transform

# ----- Wavelet Packet Decomposition (WPD) -----
y = swpd(x, wt)     # Decompose into L levels
nothing # hide
```

### [Autocorrelation Wavelet Transforms] (@id ac_transforms)
The autocorrelation wavelet transforms, similar to the stationary transforms, is a redundant type of wavelet transform. This also means that there is no downsampling involved unlike the standard transforms, and that more information is retained. 

While the decomposition process is still slower than that of the standard transform, its reconstruction process is extremely quick as it only requires the iterative summation of approximation and detail coefficients.
```@example wt
# ----- Discrete Wavelet Transform (DWT) -----
y = acdwt(x, wt)    # Forward transform
z = iacdwt(y)       # Inverse transform

# ----- Wavelet Packet Decomposition (WPD) -----
y = acwpd(x, wt)    # Decompose into L levels
nothing # hide
```

### Comparisons between standard, autocorrelation, and stationary wavelet transforms
- 1D Example:
```@example dwt
using Wavelets, WaveletsExt, Plots

# Define signal and wavelet
x = zeros(256); x[128] = 1;
wt = wavelet(WT.db4);

# Wavelet transforms
xw0 = dwt(x, wt, 4);
xw0 = [repeat(xw0[1:16],inner=16) repeat(xw0[17:32], inner=16) repeat(xw0[33:64], inner=8) repeat(xw0[65:128], inner=4) repeat(xw0[129:256], inner=2)]; nothing # hide
xw1 = sdwt(x, wt, 4);
xw2 = acdwt(x, wt, 4);

# Wiggle plots
p0 = wiggle(xw0, Overlap=false) 
plot!(p0, yticks=1:5, title="Standard WT")
p1 = wiggle(xw1, Overlap=false) 
plot!(p1, yticks=1:5, title="Stationary WT")
p2 = wiggle(xw2, Overlap=false)
plot!(p2, yticks=1:5, title="Autocorrelation WT")
plot(p0, p1, p2, layout=(1,3))
```

- 2D Example
```@example dwt
using TestImages

img = testimage("cameraman");
x = convert(Array{Float64}, img);

# Wavelet Transforms
xw0 = dwt(x, wt, 1);
xw1 = sdwt(x, wt, 1);
xw2 = acdwt(x, wt, 1);

# Outputs
p0 = heatmap(xw0, yflip=true, color=:greys, legend=false, xaxis=false, yaxis=false, xticks=false, yticks=false);
plot!(p0, title="Standard WT")
p1 = heatmap([xw1[:,:,1] xw1[:,:,2]; xw1[:,:,3] xw1[:,:,4]], yflip=true, color=:greys, legend=false, xaxis=false, yaxis=false, xticks=false, yticks=false)
plot!(p1, title="Stationary WT")
p2 = heatmap([xw2[:,:,1] xw2[:,:,2]; xw2[:,:,3] xw2[:,:,4]], yflip=true, color=:greys, legend=false, xaxis=false, yaxis=false, xticks=false, yticks=false)
plot!(p2, title="Autocorrelation WT")
plot(plot(img, title="Original"), p0, p1, p2, layout=(2,2))
```

### [Shift Invariant Wavelet Packet Decomposition] (@id si_transforms)
The [Shift-Invariant Wavelet Decomposition (SIWPD)](https://israelcohen.com/wp-content/uploads/2018/05/ICASSP95.pdf) is developed by Cohen et. al.. While it is also a type of redundant transform, it does not follow the same methodology as the SWT and the ACWT. Cohen's main goal for developing this algorithm was to obtain a global minimum entropy from a signal and all its shifted versions. See [its best basis implementation](@ref si_bestbasis) for more information.

One can compute the SIWPD of a single signal as follows.
```@example wt
# decomposition
xw = siwpd(x, wt);
nothing # hide
```

!!! note 
    As of right now, there is not too many functions written based on the SIWPD, as it does not follow the conventional style of wavelet transforms. There is a lot of ongoing work to develop more functions catered for the SIWPD such as it's inverse transforms and group-implementations.




# Stationary Wavelet Transform

```@index
Modules = [SWT]
```

## Transforms on 1 Signal
```@docs
SWT.sdwt
SWT.sdwt!
SWT.isdwt
SWT.isdwt!
SWT.swpt
SWT.swpt!
SWT.iswpt
SWT.iswpt!
SWT.swpd
SWT.swpd!
SWT.iswpd
SWT.iswpd!
```

## Transforms on Multiple Signals
```@docs
SWT.sdwtall
SWT.isdwtall
SWT.swptall
SWT.iswptall
SWT.swpdall
SWT.iswpdall
```

## Single Step Transforms
```@docs
SWT.sdwt_step
SWT.sdwt_step!
SWT.isdwt_step
SWT.isdwt_step!
```# Utils

```@index
Modules = [Utils]
```

## Useful wavelet/signal utilities
```@docs
Wavelets.Util.maxtransformlevels
Utils.getbasiscoef
Utils.getbasiscoefall
Utils.getrowrange
Utils.getcolrange
Utils.nodelength
Utils.coarsestscalingrange
Utils.finestdetailrange
```

## Tree traversing functions
```@docs
Wavelets.Util.isvalidtree
Wavelets.Util.maketree
Utils.getchildindex
Utils.getparentindex
Utils.getleaf
Utils.getdepth
Utils.gettreelength
```

## Metrics
```@docs
Utils.relativenorm
Utils.psnr
Utils.snr
Utils.ssim
```

## Dataset generation
```@docs
Utils.ClassData
Utils.duplicatesignals
Utils.generatesignals
Utils.generateclassdata
```

## Miscellaneous
```@docs
Utils.main2depthshift
```# Best Basis

```@index
Modules = [BestBasis]
```

## Cost functions and computations
```@docs
BestBasis.CostFunction
BestBasis.LSDBCost
BestBasis.JBBCost
BestBasis.BBCost
BestBasis.LoglpCost
BestBasis.NormCost
BestBasis.DifferentialEntropyCost
BestBasis.ShannonEntropyCost
BestBasis.LogEnergyEntropyCost
BestBasis.coefcost
BestBasis.tree_costs
BestBasis.tree_costs(::AbstractMatrix{T}, ::AbstractVector{BitVector}, ::SIBB) where T<:Number
```

# Best Basis Tree Selection
```@docs
BestBasis.bestbasis_treeselection
BestBasis.bestbasis_treeselection(::AbstractVector{Tc}, ::AbstractVector{Tt}) where {Tc<:AbstractVector{<:Union{Number,Nothing}}, Tt<:BitVector}
BestBasis.delete_subtree!
```

## Best basis computation
```@docs
BestBasis.BestBasisType
BestBasis.LSDB
BestBasis.JBB
BestBasis.BB
BestBasis.SIBB
Wavelets.Threshold.bestbasistree
Wavelets.Threshold.bestbasistree(::AbstractMatrix{T}, ::Integer, ::SIBB) where T<:Number
BestBasis.bestbasistreeall
```
# [Denoising](@id denoising_api)

```@index
Modules = [Denoising]
```

## Shrinking Types and Constructors
```@docs
Denoising.RelErrorShrink
Denoising.SureShrink
Denoising.SureShrink(::AbstractArray{T}, ::Bool, ::Union{BitVector, Nothing}, ::Wavelets.Threshold.THType) where T<:Number
Wavelets.Threshold.VisuShrink
```

## Threshold Determination and Noise Estimation
```@docs
Wavelets.Threshold.noisest
Denoising.relerrorthreshold
```

## Denoising Functions
```@docs
Wavelets.Threshold.denoise
Denoising.denoiseall
```

## Helper Functions for Threshold Determination and Noise Estimation
```@docs
Denoising.surethreshold
Denoising.orth2relerror
Denoising.findelbow
Denoising.relerrorplot
```# Visualizations

```@index
Modules = [Visualizations]
```

## Public API
```@autodocs
Modules = [Visualizations]
Private = false
```

## Private API
```@autodocs
Modules = [Visualizations]
Public = false
```# Standard Wavelet Transforms

```@index
Modules = [DWT]
```

## Transforms on 1 Signal
```@docs
Wavelets.Transforms.wpt
Wavelets.Transforms.wpt(::AbstractArray{T,2}, ::OrthoFilter, ::Integer; ::Bool) where T<:Number
Wavelets.Transforms.wpt!
Wavelets.Transforms.wpt!(::AbstractArray{T,2}, ::AbstractArray{T,2}, ::OrthoFilter, ::Integer; ::Bool) where T<:Number
Wavelets.Transforms.iwpt
Wavelets.Transforms.iwpt(::AbstractArray{T,2}, ::OrthoFilter, ::Integer; ::Bool) where T<:Number
Wavelets.Transforms.iwpt!
Wavelets.Transforms.iwpt!(::AbstractArray{T,2}, ::AbstractArray{T,2}, ::OrthoFilter, ::Integer; ::Bool) where T<:Number
DWT.wpd
DWT.wpd!
DWT.iwpd
DWT.iwpd!
```

## Transforms on Multiple Signals
```@docs
DWT.dwtall
DWT.idwtall
DWT.wptall
DWT.iwptall
DWT.wpdall
DWT.iwpdall
```

## Single Step Transforms
```@docs
DWT.dwt_step
DWT.dwt_step!
DWT.idwt_step
DWT.idwt_step!
```# Autocorrelation Wavelet Transform

```@index
Modules = [ACWT]
```

## Transforms on 1 Signal
```@docs
ACWT.acdwt
ACWT.acdwt!
ACWT.iacdwt
ACWT.iacdwt!
ACWT.acwpt
ACWT.acwpt!
ACWT.iacwpt
ACWT.iacwpt!
ACWT.acwpd
ACWT.acwpd!
ACWT.iacwpd
ACWT.iacwpd!
```

## Transforms on Multiple Signals
```@docs
ACWT.acdwtall
ACWT.iacdwtall
ACWT.acwptall
ACWT.iacwptall
ACWT.acwpdall
ACWT.iacwpdall
```

## Utilities
```@docs
ACWT.autocorr
ACWT.pfilter
ACWT.qfilter
ACWT.make_acqmfpair
ACWT.make_acreverseqmfpair
```

## Single Step Transforms
```@docs
ACWT.acdwt_step
ACWT.acdwt_step!
ACWT.iacdwt_step
ACWT.iacdwt_step!
```# Shift Invariant Wavelet Packet Decomposition

```@index
Modules = [SIWPD]
```

## Public API
```@autodocs
Modules = [SIWPD]
Private = false
```
# Local Discriminant Basis

```@index
Modules = [LDB]
```

## Energy Maps
```@docs
LDB.EnergyMap
LDB.TimeFrequency
LDB.ProbabilityDensity
LDB.Signatures
LDB.energy_map
```

## Discriminant Measures
```@docs
LDB.DiscriminantMeasure
LDB.ProbabilityDensityDM
LDB.SignaturesDM
LDB.AsymmetricRelativeEntropy
LDB.SymmetricRelativeEntropy
LDB.HellingerDistance
LDB.LpDistance
LDB.EarthMoverDistance
LDB.discriminant_measure
LDB.pairwise_discriminant_measure
```

## Computation of Discriminant Powers
```@docs
LDB.DiscriminantPower
LDB.BasisDiscriminantMeasure
LDB.FishersClassSeparability
LDB.RobustFishersClassSeparability
LDB.discriminant_power
```

## Feature Extraction and Transformation
```@docs
LDB.LocalDiscriminantBasis
LDB.fit!
LDB.fitdec!
LDB.fit_transform
LDB.transform
LDB.inverse_transform
LDB.change_nfeatures
```
