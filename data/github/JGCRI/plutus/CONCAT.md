[![docs](https://github.com/JGCRI/plutus/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/JGCRI/plutus/actions/workflows/pkgdown.yaml)
[![build](https://github.com/JGCRI/plutus/actions/workflows/rcmd.yml/badge.svg)](https://github.com/JGCRI/plutus/actions/workflows/rcmd.yml)
[![codecov](https://codecov.io/gh/JGCRI/plutus/branch/main/graph/badge.svg?token=1PK34KIHKE)](https://codecov.io/gh/JGCRI/plutus)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03212/status.svg)](https://doi.org/10.21105/joss.03212)

# plutus
Plutus is designed for GCAM v5.3 (excluding GCAM-USA).
<br />

<!-------------------------->
<!-------------------------->
# <a name="Contents"></a>Contents
<!-------------------------->
<!-------------------------->

- [Key Links](#KeyLinks)
- [Introduction](#Introduction)
- [Citation](#Citation)
- [Installation Guide](#InstallGuides)
- [How-to Guides](#How-toGuides) 

<br />

<!-------------------------->
<!-------------------------->
# <a name="KeyLinks"></a>Key Links
<!-------------------------->
<!-------------------------->

- Github: https://github.com/JGCRI/plutus
- Webpage: https://jgcri.github.io/plutus/

[Back to Contents](#Contents)

<br />

<!-------------------------->
<!-------------------------->
# <a name="Introduction"></a>Introduction
<!-------------------------->
<!-------------------------->

`plutus` post-processes outputs from the Global Change Analysis Model (GCAM) to calculate the electricity investment costs and stranded asset costs associated with GCAM projections of future power sector energy generation by technology.


[Back to Contents](#Contents)

<br />

<!-------------------------->
<!-------------------------->
# <a name="Citation"></a>Citation
<!-------------------------->
<!-------------------------->

Zhao, M., Binsted, M., Wild, T.B., Khan, Z., Yarlagadda, B., Iyer, G., Vernon, C., Patel, P., Santos da Silva, S.R., Calvin, K.V., (2021). plutus - An R package to calculate electricity investments and stranded assets from the Global Change Analysis Model (GCAM). Journal of Open Source Software, 6(65), 3212, https://doi.org/10.21105/joss.03212


[Back to Contents](#Contents)

<br />


<!-------------------------->
<!-------------------------->
# <a name="InstallationGuides"></a>Installation Guides
<!-------------------------->
<!-------------------------->

1. Download and install:

    - R (https://www.r-project.org/)
    - R studio (https://www.rstudio.com/)

2. For Linux users, install following libraries:

```
sudo apt install build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev
sudo apt-get install libxml2-dev
```
    
3. Open R studio:

```
install.packages('devtools')
devtools::install_github('JGCRI/rgcam')
devtools::install_github('JGCRI/plutus')
```

4. `Metis` installation

`Metis` provides functions to visualize the outputs from `plutus`. The installation guide for `Metis` can be accessed at [Metis Github Page](https://github.com/JGCRI/metis).

[Back to Contents](#Contents)

<br />


<!-------------------------->
<!-------------------------->
# <a name="How-toGuides"></a>How-to Guides
<!-------------------------->
<!-------------------------->
`plutus::gcamInvest` provides all-in-one workflow that reads gcamdata, processes queries, and estimates stranded assets and capital investments. Please visit the followings for detailed instructions.

- [Instruction on `plutus::gcamInvest`](https://jgcri.github.io/plutus/articles/gcamInvest.html)
- [Case tutorial](https://jgcri.github.io/plutus/articles/CaseTutorial.html)

[Back to Contents](#Contents)

<br />
<!-- ------------------------>
<!-- ------------------------>
# plutus 0.1.0 (Under development)
<p align="center"> <img src="READMEfigs/plutus.PNG"></p>
<!-- ------------------------>
<!-- ------------------------>

## New features

## Bug Fixes
# How to Contribute to `plutus`

Thank you for taking the time to contribute and help us advance the science and architecture of `plutus`. We provide few guidelines that we ask contributors to follow. The guidelines aim to ease the maintainers' organizational and logistical duties, while encouraging development by others.

Before you start:

* Make sure you have a [GitHub account](https://github.com/signup/free).
* Trivial changes to comments or documentation do not require creating a new issue.

## Did you find a bug?

* Make sure the bug was not already reported in the Github [Issues](https://github.com/JGCRI/plutus/issues).
* [Open an issue](https://github.com/JGCRI/plutus/issues/new) and clearly describe the issue with as much information as possible. A code sample or an executable test case are recommended.
  
## Did you plan to write a patch that fixes a bug?

  * [Open an issue](https://github.com/JGCRI/plutus/issues/new) and clearly describes the problem and discuss how your solution will affect `plutus`.
  * Fork the repository on GitHub to work on the patch.
  * Interact with the project maintainers to refine/change/prioritize your issue.

## Making changes

* Start your work on your fork of the repository.
* Check for unnecessary whitespace with `git diff --check` and format code.
* Make sure your commit messages are descriptive but succinct, describing what was changed and why, and **reference the relevant issue number**. Make commits of logical units.
* Make sure you have added the necessary tests for your changes.
* Run **all** the tests to assure nothing else was accidentally broken.

## Submitting changes

* Submit a pull request with clear documentation of the methodology to the main `plutus` repository.
* **Your pull request should include one of the following two statements**:
   * You own the copyright on the code being contributed, and you hereby grant PNNL unlimited license to use this code in this version or any future version of `plutus`. You reserve all other rights to the code.
   * Somebody else owns the copyright on the code being contributed (e.g., your employer because you did it as part of your work for them); you are authorized by that owner to grant PNNL an unlimited license to use this code in this version or any future version of `plutus`, and you hereby do so. All other rights to the code are reserved by the copyright owner.
* The core team looks at Pull Requests, and will respond as soon as possible.

# Additional Resources

* [General GitHub documentation](http://help.github.com/)
* [GitHub pull request documentation](http://help.github.com/send-pull-requests/)
