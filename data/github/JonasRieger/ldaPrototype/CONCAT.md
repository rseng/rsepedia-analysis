# ldaPrototype
[![status](https://joss.theoj.org/papers/ecce89f69453dc2ee0c697fbcb776de8/status.svg)](https://joss.theoj.org/papers/10.21105/joss.02181)
[![CRAN](https://www.r-pkg.org/badges/version/ldaPrototype)](https://cran.r-project.org/package=ldaPrototype)
[![R build status](https://github.com/JonasRieger/ldaPrototype/workflows/R-CMD-check/badge.svg)](https://github.com/JonasRieger/ldaPrototype/actions)
[![Build status](https://ci.appveyor.com/api/projects/status/6nm48mamn79n6xf4?svg=true)](https://ci.appveyor.com/project/JonasRieger/ldaprototype)
[![codecov](https://codecov.io/gh/JonasRieger/ldaPrototype/branch/master/graph/badge.svg?token=IUWLTOW2HV)](https://codecov.io/gh/JonasRieger/ldaPrototype)
[![DOI](https://zenodo.org/badge/187803702.svg)](https://zenodo.org/badge/latestdoi/187803702)

## Prototype of Multiple Latent Dirichlet Allocation Runs
Determine a Prototype from a number of runs of Latent Dirichlet Allocation (LDA) measuring its similarities with S-CLOP: A procedure to select the LDA run with highest mean pairwise similarity, which is measured by S-CLOP (Similarity of multiple sets by Clustering with Local Pruning), to all other runs. LDA runs are specified by its assignments leading to estimators for distribution parameters. Repeated runs lead to different results, which we encounter by choosing the most representative LDA run as prototype.

## Citation
Please cite the [JOSS](https://doi.org/10.21105/joss.02181) paper using the BibTeX entry
```
@article{<placeholder>,
    title = {{ldaPrototype}: A method in {R} to get a Prototype of multiple Latent Dirichlet Allocations},
    author = {Jonas Rieger},
    journal = {Journal of Open Source Software},
    year = {2020},
    volume = {5},
    number = {51},
    pages = {2181},
    doi = {10.21105/joss.02181},
    url = {https://doi.org/10.21105/joss.02181}
  }
```
which is also obtained by the call ``citation("ldaPrototype")``.

## References
* Rieger, J., Jentsch, C. & Rahnenführer, J.: LDAPrototype: A Model Selection Algorithm to Improve Reliability of Latent Dirichlet Allocation. [working paper](https://www.statistik.tu-dortmund.de/fileadmin/user_upload/rieger_ldaproto.pdf)
* Rieger, J. (2020). ldaPrototype: A method in R to get a Prototype of multiple Latent Dirichlet Allocations. [Journal of Open Source Software](https://doi.org/10.21105/joss.02181), 5(51), 2181.
* Rieger, J., Rahnenführer, J. & Jentsch, C. (2020). Improving Latent Dirichlet Allocation: On Reliability of the Novel Method LDAPrototype. [Natural Language Processing and Information Systems, NLDB 2020.](https://doi.org/10.1007/978-3-030-51310-8_11) LNCS 12089, pp. 118-125.
* Rieger, J., Koppers, L., Jentsch, C. & Rahnenführer, J. (2020). Improving Reliability of Latent Dirichlet Allocation by Assessing Its Stability using Clustering Techniques on Replicated Runs. [arXiv](https://arxiv.org/abs/2003.04980)

## Related Software
* [tm](https://CRAN.R-project.org/package=tm) is useful for preprocessing text data.
* [lda](https://CRAN.R-project.org/package=lda) offers a fast implementation of the Latent Dirichlet Allocation and is used by ``ldaPrototype``.
* [quanteda](https://quanteda.io/) is a framework for "Quantitative Analysis of Textual Data".
* [stm](https://www.structuraltopicmodel.com/) is a framework for Structural Topic Models.
* [tosca](https://github.com/Docma-TU/tosca) is a framework for statistical methods in content analysis including visualizations and validation techniques. It is also useful for managing and manipulating text data to a structure requested by ``ldaPrototype``.
* [topicmodels](https://CRAN.R-project.org/package=topicmodels) is another framework for various topic models based on the Latent Dirichlet Allocation and Correlated Topics Models.
* [ldatuning](https://github.com/nikita-moor/ldatuning) is a framework for finding the optimal number of topics using various metrics.

## Contribution
This R package is licensed under the [GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html).
For bug reports (lack of documentation, misleading or wrong documentation, unexpected behaviour, ...) and feature requests please use the [issue tracker](https://github.com/JonasRieger/ldaPrototype/issues).
Pull requests are welcome and will be included at the discretion of the author.

## Installation
```{R}
install.packages("ldaPrototype")
```
For the development version use [devtools](https://cran.r-project.org/package=devtools):
```{R}
devtools::install_github("JonasRieger/ldaPrototype")
```

## (Quick Start) Example
Load the package and the example dataset from Reuters consisting of 91 articles - [tosca::LDAprep](https://github.com/Docma-TU/tosca/blob/master/R/LDAprep.R) can be used to manipulate text data to the format requested by ``ldaPrototype``.
```{R}
library("ldaPrototype")
data(reuters_docs)
data(reuters_vocab)
```
Run the shortcut function to create a LDAPrototype object. It consists of the LDAPrototype of 4 LDA runs (with specified seeds) with 10 topics each. The LDA selected by the algorithm can be retrieved using ``getPrototype`` or ``getLDA``.
```{R}
res = LDAPrototype(docs = reuters_docs, vocabLDA = reuters_vocab, n = 4, K = 10, seeds = 1:4)
proto = getPrototype(res) #= getLDA(res)
```
The same result can also be achieved by executing the following lines of code in several steps, which can be useful for interim evaluations.
```{R}
reps = LDARep(docs = reuters_docs, vocab = reuters_vocab,
  n = 4, K = 10, seeds = 1:4)
topics = mergeTopics(reps, vocab = reuters_vocab)
jacc = jaccardTopics(topics)
sclop = SCLOP.pairwise(jacc)
res2 = getPrototype(reps, sclop = sclop)

proto2 = getPrototype(res2) #= getLDA(res2)

identical(res, res2)
```
There is also the option to use similarity measures other than the Jaccard coefficient. Currently, the measures cosine similarity (``cosineTopics``), Jensen-Shannon divergence (``jsTopics``) and rank-biased overlap (``rboTopics``) are implemented in addition to the standard Jaccard coefficient (``jaccardTopics``).

To get an overview of the workflow, the associated functions and getters for each type of object, the following call is helpful:
```{R}
?`ldaPrototype-package`
```

## (Slightly more detailed) Example
Similar to the quick start example, the shortcut of one single call is again compared with the step-by-step procedure. 
We model 5 LDAs with ``K = 12`` topics, hyperparameters ``alpha = eta = 0.1`` and seeds ``1:5``. We want to calculate the log likelihoods for the 20 iterations after 5 burn-in iterations and topic similarities should be based on ``atLeast = 3`` words (see Step 3 below). In addition, we want to keep all interim calculations, which would be discarded by default to save memory space.
```{R}
res = LDAPrototype(docs = reuters_docs, vocabLDA = reuters_vocab,
  n = 5, K = 12, alpha = 0.1, eta = 0.1, compute.log.likelihood = TRUE,
  burnin = 5, num.iterations = 20, atLeast = 3, seeds = 1:5,
  keepLDAs = TRUE, keepSims = TRUE, keepTopics = TRUE)
```
Based on ``res`` we can have a look at several getter functions:
```{R}
getID(res)
getPrototypeID(res)

getParam(res)
getParam(getLDA(res))

getLDA(res, all = TRUE)
getLDA(res)

est = getEstimators(getLDA(res))
est$phi[,1:3]
est$theta[,1:3]
getLog.likelihoods(getLDA(res))

getSCLOP(res)
getSimilarity(res)[1:5, 1:5]
tosca::topWords(getTopics(getLDA(res)), 5)
```
#### Step 1: LDA Replications
In the first step we simply run the LDA procedure five times with the given parameters. This can also be done with support of [batchtools](https://github.com/mllg/batchtools) using ``LDABatch`` instead of ``LDARep`` or [parallelMap](https://github.com/mlr-org/parallelMap) setting the ``pm.backend`` and (optionally) ``ncpus`` argument(s).
```{R}
reps = LDARep(docs = reuters_docs, vocab = reuters_vocab,
  n = 5, K = 12, alpha = 0.1, eta = 0.1, compute.log.likelihood = TRUE,
  burnin = 5, num.iterations = 20, seeds = 1:5)
```
#### Step 2: Merging Topic Matrices of Replications
The topic matrices of all replications are merged and reduced to the vocabulary given in ``vocab``. By default the vocabulary of the first topic matrix is used as a simplification of the case that all LDAs contain the same vocabulary set.
```{R}
topics = mergeTopics(reps, vocab = reuters_vocab)
```
#### Step 3: Topic Similarities
We use the merged topic matrix to calculate pairwise topic similarites using the Jaccard coefficient with parameters adjusting the consideration of words. A word is taken as relevant for a topic if its count passes thresholds given by ``limit.rel`` and ``limit.abs``. A word is considered for calculation of similarities if it's relevant for the topic or if it belongs to the (``atLeast =``) 3 most common words in the corresponding topic. Alternatively, the similarities can also be calculated considering the cosine similarity (``cosineTopics``), Jensen-Shannon divergence (``jsTopics`` - parameter ``epsilon`` to ensure computability) or rank-biased overlap (``rboTopics`` - parameter ``k`` for maximum depth of evaluation and ``p`` as weighting parameter).
```{R}
jacc = jaccardTopics(topics, limit.rel = 1/500, limit.abs = 10, atLeast = 3)
getSimilarity(jacc)[1:3, 1:3]
```
We can check the number of relevant and considered words using the ad-hoc getter. The difference between ``n1`` and ``n2`` can become larger than (``atLeast =``) 3 if there are ties in the count of words, which is negligible for large sample sizes.
```{R}
n1 = getRelevantWords(jacc)
n2 = getConsideredWords(jacc)
(n2-n1)[n2-n1 != 0]
```
#### Step 3.1: Representation of Topic Similarities as Dendrogram
It is possible to represent the calulcated pairwise topic similarities as dendrogram using ``dendTopics`` and related ``plot`` options.
```{R}
dend = dendTopics(jacc)
plot(dend)
```
The S-CLOP algorithm results in a pruning state of the dendrogram, which can be retrieved calling ``pruneSCLOP``. By default each of the topics is colorized by its LDA run belonging; but the cluster belongings can also be visualized by the colors or by vertical lines with freely chosen parameters.
```{R}
pruned = pruneSCLOP(dend)
plot(dend, pruned)
plot(dend, pruning = pruned, pruning.par = list(type = "both", lty = 1, lwd = 2, col = "red"))
```
#### Step 4: Pairwise LDA Model Similarities (S-CLOP)
For determination of the LDAPrototype the pairwise S-CLOP similarities of the 5 LDA runs are needed.
```{R}
sclop = SCLOP.pairwise(jacc)
```
#### Step 5: Determine LDAPrototype
In the last step the LDAPrototype itself is determined by maximizing the mean pairwise S-CLOP per LDA.
```{R}
res2 = getPrototype(reps, sclop = sclop)
```
There are several possibilites for using shortcut functions to summarize steps of the procedure. For example, we can determine the LDAPrototype after Step 1:
```{R}
res3 = getPrototype(reps, atLeast = 3)
```
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
reported by contacting the project team at jonas.rieger@tu-dortmund.de. All
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
title: 'ldaPrototype: A method in R to get a Prototype of multiple Latent Dirichlet Allocations'
authors:
- affiliation: 1
  name: Jonas Rieger
  orcid: 0000-0002-0007-4478
date: "10 March 2020"
output: pdf_document
bibliography: paper.bib
tags:
- R
- topic modeling
- natural language processing
- text data
- tuning
- model selection
- high performance computing
- parallelization
affiliations:
- index: 1
  name: TU Dortmund University
---

# Summary

Topic Modeling [@blei2012] is one of the biggest subjects in the field of text data analysis. Here, the Latent Dirichlet Allocation [@blei2003] takes a special position. A large part of scientific text data analyses are based on this model (LDA). The LDA method has a far-reaching disadvantage. Random initialization and conditional reassignments within the iterative process of the Gibbs sampler [@griffiths2004] can result in fundamentally different models when executed several times on the same data and with identical parameter sets. This fact greatly limits the scientific reproducibility.

Up to now, the so-called eye-balling method has been used in practice to select suitable results. From a set of models, subjective decisions are made to select the model that seems to fit the data best or, in the worst case, the result that best supports one's hypothesis is chosen. The latter contradicts basically good scientific practice. A different method of objective and automated selection has also become established. A model from a set of LDAs can be determined optimizing the log-likelihood using the perplexity on held-out data. The [`R`](https://www.r-project.org/) [@R] package [`topicmodels`](https://CRAN.R-project.org/package=topicmodels) [@topicmodels] provides a workflow for this procedure. As an extension, @nguyen2014 proposed to average different iterations of the Gibbs sampling procedure to achieve an increase of perplexity. The averaging technique has the weakness, that the user does not get token specific assignments to topics, but only averaged topic counts or proportions per text. In addition, @chang2009 were able to show that selection mechanisms aiming for optimizing likelihood-based measures do not correspond to the human perception of a well-adapted model of text data. Instead, the authors propose a so-called intruder procedure based on human codings. The corresponding methodology is implemented in the package [`tosca`](https://github.com/Docma-TU/tosca) [@tosca].

The [`R`](https://www.r-project.org/) package [`ldaPrototype`](https://github.com/JonasRieger/ldaPrototype) on the other hand determines a prototypical LDA by automated selection from a set of LDAs. The method improves reliability of findings drawn from LDA results [@rieger2020], which is achieved following a typical statistical approach. For a given combination of parameters, a number of models is calculated (usually about 100), from which that LDA is determined that is most similar to all other LDAs from a set of models. For this purpose pairwise model similarities are calculated using the S-CLOP measure (**S**imilarity of Multiple Sets by **C**lustering with **Lo**cal **P**runing), which can be determined by a clustering procedure of the individual topic units based on topic similarities of the two LDA results considered. The package offers visualization possibilities for comparisons of LDA models based on the clustering of the associated topics. Furthermore, the package supports the repetition of the modeling procedure of the LDA by a simple calculation of the repeated LDA runs.

In addition to the possibility of local parallel computation by connecting to the package [`parallelMap`](https://github.com/berndbischl/parallelMap) [@parallelMap], there is the possibility to calculate using batch systems on high performance computing (HPC) clusters by integrating helpful functions from the package [`batchtools`](https://github.com/mllg/batchtools) [@batchtools]. This is especially helpful if the text corpora contains several hundred of thousands articles and the sequential calculation of 100 or more LDA runs would extend over several days. The modeling of single LDA runs is done with the help of the computation time optimized [`R`](https://www.r-project.org/) package [`lda`](https://CRAN.R-project.org/package=lda) [@lda], which implements the calculation in `C++` code. In general, the package [`ldaPrototype`](https://github.com/JonasRieger/ldaPrototype) is based on `S3` objects and thus extends the packages [`lda`](https://CRAN.R-project.org/package=lda) and [`tosca`](https://github.com/Docma-TU/tosca) by user-friendly display and processing options. Other [`R`](https://www.r-project.org/) packages for estimating LDA are [`topicmodels`](https://CRAN.R-project.org/package=topicmodels) and  [`mallet`](https://github.com/mimno/RMallet) [@mallet], whereas [`stm`](https://www.structuraltopicmodel.com/) [@stm] offers a powerful framework for Structural Topic Models and [`quanteda`](https://quanteda.io/) [@quanteda] is a popular framework for preprocessing and quantitative analysis of text data.

# References
