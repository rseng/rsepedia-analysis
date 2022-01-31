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
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mergeRepTopics.R
\name{mergeRepTopics}
\alias{mergeRepTopics}
\alias{mergeRepTopics.LDARep}
\alias{mergeRepTopics.default}
\title{Merge LDA Topic Matrices}
\usage{
mergeRepTopics(...)

\method{mergeRepTopics}{LDARep}(x, vocab, progress = TRUE, ...)

\method{mergeRepTopics}{default}(lda, vocab, id, progress = TRUE, ...)
}
\arguments{
\item{...}{additional arguments}

\item{x}{[\code{named list}]\cr
\code{\link{LDARep}} object. Alternatively \code{lda} and \code{id} can be passed.}

\item{vocab}{[\code{character}]\cr
Vocabularies taken into consideration for merging topic matrices. Default is
the vocabulary of the first LDA.}

\item{progress}{[\code{logical(1)}]\cr
Should a nice progress bar be shown? Turning it off, could lead to significantly
faster calculation. Default ist \code{TRUE}.}

\item{lda}{[\code{named list}]\cr
List of \code{\link{LDA}} objects, named by the corresponding "job.id".}

\item{id}{[\code{character(1)}]\cr
Name for the computation. Default is "LDARep".}
}
\value{
[\code{named matrix}] with the count of vocabularies (row wise) in topics (column wise).
}
\description{
Collects LDA results from a list of replicated runs and merges their topic
matrices for a given set of vocabularies.
}
\details{
For details and examples see \code{\link{mergeTopics}}.
}
\seealso{
Other merge functions: 
\code{\link{mergeBatchTopics}()},
\code{\link{mergeTopics}()}

Other replication functions: 
\code{\link{LDAPrototype}()},
\code{\link{LDARep}()},
\code{\link{as.LDARep}()},
\code{\link{getJob}()}
}
\concept{merge functions}
\concept{replication functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LDABatch.R
\name{LDABatch}
\alias{LDABatch}
\title{LDA Replications on a Batch System}
\usage{
LDABatch(
  docs,
  vocab,
  n = 100,
  seeds,
  id = "LDABatch",
  load = FALSE,
  chunk.size = 1,
  resources,
  ...
)
}
\arguments{
\item{docs}{[\code{list}]\cr
Documents as received from \code{\link[tosca]{LDAprep}}.}

\item{vocab}{[\code{character}]\cr
Vocabularies passed to \code{\link[lda]{lda.collapsed.gibbs.sampler}}.
For additional (and necessary) arguments passed, see ellipsis (three-dot argument).}

\item{n}{[\code{integer(1)}]\cr
Number of Replications.}

\item{seeds}{[\code{integer(n)}]\cr
Random Seeds for each Replication.}

\item{id}{[\code{character(1)}]\cr
Name for the registry's folder.}

\item{load}{[\code{logical(1)}]\cr
If a folder with name \code{id} exists: should the existing registry be loaded?}

\item{chunk.size}{[\code{integer(1)}]\cr
Requested chunk size for each single chunk. See \code{\link[batchtools]{chunk}}.}

\item{resources}{[\code{named list}]\cr
Computational resources for the jobs to submit. See \code{\link[batchtools]{submitJobs}}.}

\item{...}{additional arguments passed to \code{\link[lda]{lda.collapsed.gibbs.sampler}}.
Arguments will be coerced to a vector of length \code{n}.
Default parameters are \code{alpha = eta = 1/K} and \code{num.iterations = 200}.
There is no default for \code{K}.}
}
\value{
[\code{named list}] with entries \code{id} for the registry's folder name,
\code{jobs} for the submitted jobs' ids and its parameter settings and
\code{reg} for the registry itself.
}
\description{
Performs multiple runs of Latent Dirichlet Allocation on a batch system using
the \code{\link[batchtools]{batchtools-package}}.
}
\details{
The function generates multiple LDA runs with the possibility of
using a batch system. The integration is done by the
\code{\link[batchtools]{batchtools-package}}. After all jobs of the
corresponding registry are terminated, the whole registry can be ported to
your local computer for further analysis.

The function returns a \code{LDABatch} object. You can receive results and
all other elements of this object with getter functions (see \code{\link{getJob}}).
}
\examples{
\dontrun{
batch = LDABatch(docs = reuters_docs, vocab = reuters_vocab, n = 4, K = 15)
batch
getRegistry(batch)
getJob(batch)
getLDA(batch, 2)

batch2 = LDABatch(docs = reuters_docs, vocab = reuters_vocab, K = 15, chunk.size = 20)
batch2
head(getJob(batch2))
}

}
\seealso{
Other batch functions: 
\code{\link{as.LDABatch}()},
\code{\link{getJob}()},
\code{\link{mergeBatchTopics}()}

Other LDA functions: 
\code{\link{LDARep}()},
\code{\link{LDA}()},
\code{\link{getTopics}()}
}
\concept{LDA functions}
\concept{batch functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getTopics.R
\name{getTopics}
\alias{getTopics}
\alias{getAssignments}
\alias{getDocument_sums}
\alias{getDocument_expects}
\alias{getLog.likelihoods}
\alias{getParam}
\alias{getK}
\alias{getAlpha}
\alias{getEta}
\alias{getNum.iterations}
\alias{getEstimators}
\title{Getter for LDA}
\usage{
getTopics(x)

getAssignments(x)

getDocument_sums(x)

getDocument_expects(x)

getLog.likelihoods(x)

getParam(x)

getK(x)

getAlpha(x)

getEta(x)

getNum.iterations(x)

getEstimators(x)
}
\arguments{
\item{x}{[\code{named list}]\cr
\code{\link{LDA}} object.}
}
\description{
Returns the corresponding element of a \code{\link{LDA}} object.
\code{getEstimators} computes the estimators for \code{phi} and \code{theta}.
}
\details{
The estimators for \code{phi} and \code{theta} in
\deqn{w_n^{(m)} \mid T_n^{(m)}, \bm\phi_k  \sim \textsf{Discrete}(\bm\phi_k),}
\deqn{\bm\phi_k  \sim \textsf{Dirichlet}(\eta),}
\deqn{T_n^{(m)} \mid \bm\theta_m  \sim \textsf{Discrete}(\bm\theta_m),}
\deqn{\bm\theta_m  \sim \textsf{Dirichlet}(\alpha)}
are calculated referring to Griffiths and Steyvers (2004) by
\deqn{\hat{\phi}_{k, v} = \frac{n_k^{(v)} + \eta}{n_k + V \eta},}
\deqn{\hat{\theta}_{m, k} = \frac{n_k^{(m)} + \alpha}{N^{(m)} + K \alpha}}
with \eqn{V} is the vocabulary size, \eqn{K} is the number of modeled topics;
\eqn{n_k^{(v)}} is the count of assignments of the \eqn{v}-th word to
the \eqn{k}-th topic. Analogously, \eqn{n_k^{(m)}} is the count of assignments
of the \eqn{m}-th text to the \eqn{k}-th topic. \eqn{N^{(m)}} is the total
number of assigned tokens in text \eqn{m} and \eqn{n_k} the total number of
assigned tokens to topic \eqn{k}.
}
\references{
Griffiths, Thomas L. and Mark Steyvers (2004). "Finding scientific topics".
In: \emph{Proceedings of the National Academy of Sciences} \bold{101} (suppl 1), pp.5228--5235,
\doi{10.1073/pnas.0307752101}.
}
\seealso{
Other getter functions: 
\code{\link{getJob}()},
\code{\link{getSCLOP}()},
\code{\link{getSimilarity}()}

Other LDA functions: 
\code{\link{LDABatch}()},
\code{\link{LDARep}()},
\code{\link{LDA}()}
}
\concept{LDA functions}
\concept{getter functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pruneSCLOP.R
\name{pruneSCLOP}
\alias{pruneSCLOP}
\alias{plot.PruningSCLOP}
\alias{pruning.par}
\title{Local Pruning State of Topic Dendrograms}
\usage{
pruneSCLOP(dend)

\method{plot}{PruningSCLOP}(x, dend, pruning.par, ...)

pruning.par(pruning)
}
\arguments{
\item{dend}{[\code{\link[stats]{dendrogram}}]\cr
\code{\link[=dendTopics]{TopicDendrogram}}
(and \code{\link[stats]{dendrogram}}) object of all considered topics as the
output from \code{\link{dendTopics}}.}

\item{x}{an R object.}

\item{pruning.par}{[\code{list}]\cr
List of parameters to mark the pruning. See section "Details" at \code{\link{dendTopics}}
for default parameters. Types for marking the pruning state are \code{"abline"},
\code{"color"} and \code{"both"}.}

\item{...}{additional arguments.}

\item{pruning}{[\code{list of \link[stats]{dendrogram}s}]\cr
\code{\link[=pruneSCLOP]{PruningSCLOP}} object specifying the best possible
local pruning state.}
}
\value{
[\code{list of \link[stats]{dendrogram}s}]
\code{\link[=pruneSCLOP]{PruningSCLOP}} object specifying the best possible
local pruning state.
}
\description{
The function \code{\link{SCLOP}} calculates the S-CLOP value for the best possible
local pruning state of a dendrogram from \code{\link{dendTopics}}.
The function \code{pruneSCLOP} supplies the corresponding pruning state itself.
}
\details{
For details of computing the S-CLOP values see \code{\link{SCLOP}}.

For details and examples of plotting the pruning state see \code{\link{dendTopics}}.
}
\seealso{
Other plot functions: 
\code{\link{dendTopics}()}

Other SCLOP functions: 
\code{\link{SCLOP}()}
}
\concept{SCLOP functions}
\concept{plot functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{reuters}
\alias{reuters}
\alias{reuters_docs}
\alias{reuters_vocab}
\alias{docs}
\alias{vocab}
\title{A Snippet of the Reuters Dataset}
\format{
\code{reuters_docs} is a list of documents of length 91 prepared by \code{\link[tosca]{LDAprep}}.

\code{reuters_vocab} is

An object of class \code{character} of length 2141.
}
\source{
temporarily unavailable: http://ronaldo.cs.tcd.ie/esslli07/data/reuters21578-xml/
}
\usage{
data(reuters_docs)

data(reuters_vocab)
}
\description{
Example Dataset from Reuters consisting of 91 articles. It can be used to
familiarize with the bunch of functions offered by this package.
}
\references{
Lewis, David (1997). \emph{Reuters-21578 Text Categorization Collection Distribution 1.0}.
\url{http://kdd.ics.uci.edu/databases/reuters21578/reuters21578.html}

Luz, Saturnino. \emph{XML-encoded version of Reuters-21578}.
http://ronaldo.cs.tcd.ie/esslli07/data/reuters21578-xml/ (temporarily unavailable)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jaccardTopics.R
\name{jaccardTopics}
\alias{jaccardTopics}
\title{Pairwise Jaccard Coefficients}
\usage{
jaccardTopics(
  topics,
  limit.rel,
  limit.abs,
  atLeast,
  progress = TRUE,
  pm.backend,
  ncpus
)
}
\arguments{
\item{topics}{[\code{named matrix}]\cr
The counts of vocabularies/words (row wise) in topics (column wise).}

\item{limit.rel}{[0,1]\cr
A relative lower bound limit for which words are taken into account. Those words
are taken as relevant for a topic that have a count higher than \code{limit.rel}
multiplied by the total count of the given topic. Default is \code{1/500}.}

\item{limit.abs}{[\code{integer(1)}]\cr
An absolute lower bound limit for which words are taken into account. All words
are taken as relevant for a topic that have a count higher than \code{limit.abs}.
Default is \code{10}.}

\item{atLeast}{[\code{integer(1)}]\cr
An absolute count of how many words are at least considered as relevant for a topic.
Default is \code{0}.}

\item{progress}{[\code{logical(1)}]\cr
Should a nice progress bar be shown? Turning it off, could lead to significantly
faster calculation. Default is \code{TRUE}.
If \code{pm.backend} is set, parallelization is done and no progress bar will be shown.}

\item{pm.backend}{[\code{character(1)}]\cr
One of "multicore", "socket" or "mpi".
If \code{pm.backend} is set, \code{\link[parallelMap]{parallelStart}} is
called before computation is started and \code{\link[parallelMap]{parallelStop}}
is called after.}

\item{ncpus}{[\code{integer(1)}]\cr
Number of (physical) CPUs to use. If \code{pm.backend} is passed,
default is determined by \code{\link[future]{availableCores}}.}
}
\value{
[\code{named list}] with entries
\describe{
  \item{\code{sims}}{[\code{lower triangular named matrix}] with all pairwise
  jaccard similarities of the given topics.}
  \item{\code{wordslimit}}{[\code{integer}] with counts of words determined as
  relevant based on \code{limit.rel} and \code{limit.abs}.}
  \item{\code{wordsconsidered}}{[\code{integer}] with counts of considered
  words for similarity calculation. Could differ from \code{wordslimit}, if
  \code{atLeast} is greater than zero.}
  \item{\code{param}}{[\code{named list}] with parameter specifications for
  \code{type} [\code{character(1)}] \code{= "Jaccard Coefficient"},
  \code{limit.rel} [0,1], \code{limit.abs} [\code{integer(1)}] and
  \code{atLeast} [\code{integer(1)}]. See above for explanation.}
}
}
\description{
Calculates the similarity of all pairwise topic combinations using a modified
Jaccard Coefficient.
}
\details{
The modified Jaccard Coefficient for two topics \eqn{\bm z_{i}} and
\eqn{\bm z_{j}} is calculated by
\deqn{J_m(\bm z_{i}, \bm z_{j} \mid \bm c) = \frac{\sum_{v = 1}^{V} 1_{\left\{n_{i}^{(v)} > c_i ~\wedge~ n_{j}^{(v)} > c_j\right\}}\left(n_{i}^{(v)}, n_{j}^{(v)}\right)}{\sum_{v = 1}^{V} 1_{\left\{n_{i}^{(v)} > c_i ~\vee~ n_{j}^{(v)} > c_j\right\}}\left(n_{i}^{(v)}, n_{j}^{(v)}\right)}}
with \eqn{V} is the vocabulary size and \eqn{n_k^{(v)}} is the count of
assignments of the \eqn{v}-th word to the \eqn{k}-th topic. The threshold vector \eqn{\bm c}
is determined by the maximum threshold of the user given lower bounds \code{limit.rel}
and \code{limit.abs}. In addition, at least \code{atLeast} words per topic are
considered for calculation. According to this, if there are less than
\code{atLeast} words considered as relevant after applying \code{limit.rel}
and \code{limit.abs} the \code{atLeast} most common words per topic are taken
to determine topic similarities.

The procedure of determining relevant words is executed for each topic individually.
The values \code{wordslimit} and \code{wordsconsidered} describes the number
of relevant words per topic.
}
\examples{
res = LDARep(docs = reuters_docs, vocab = reuters_vocab, n = 4, K = 10, num.iterations = 30)
topics = mergeTopics(res, vocab = reuters_vocab)
jacc = jaccardTopics(topics, atLeast = 2)
jacc

n1 = getConsideredWords(jacc)
n2 = getRelevantWords(jacc)
(n1 - n2)[n1 - n2 != 0]

sim = getSimilarity(jacc)
dim(sim)

# Comparison to Cosine and Jensen-Shannon (more interesting on large datasets)
cosine = cosineTopics(topics)
js = jsTopics(topics)

sims = list(jaccard = sim, cosine = getSimilarity(cosine), js = getSimilarity(js))
pairs(do.call(cbind, lapply(sims, as.vector)))

}
\seealso{
Other TopicSimilarity functions: 
\code{\link{cosineTopics}()},
\code{\link{dendTopics}()},
\code{\link{getSimilarity}()},
\code{\link{jsTopics}()},
\code{\link{rboTopics}()}

Other workflow functions: 
\code{\link{LDARep}()},
\code{\link{SCLOP}()},
\code{\link{dendTopics}()},
\code{\link{getPrototype}()},
\code{\link{mergeTopics}()}
}
\concept{TopicSimilarity functions}
\concept{workflow functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getSimilarity.R
\name{getSimilarity}
\alias{getSimilarity}
\alias{getRelevantWords}
\alias{getConsideredWords}
\alias{getParam.TopicSimilarity}
\title{Getter for TopicSimilarity}
\usage{
getSimilarity(x)

getRelevantWords(x)

getConsideredWords(x)

\method{getParam}{TopicSimilarity}(x)
}
\arguments{
\item{x}{[\code{named list}]\cr
\code{\link[=jaccardTopics]{TopicSimilarity}} object.}
}
\description{
Returns the corresponding element of a \code{\link[=jaccardTopics]{TopicSimilarity}} object.
}
\seealso{
Other getter functions: 
\code{\link{getJob}()},
\code{\link{getSCLOP}()},
\code{\link{getTopics}()}

Other TopicSimilarity functions: 
\code{\link{cosineTopics}()},
\code{\link{dendTopics}()},
\code{\link{jaccardTopics}()},
\code{\link{jsTopics}()},
\code{\link{rboTopics}()}
}
\concept{TopicSimilarity functions}
\concept{getter functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LDARep.R
\name{LDARep}
\alias{LDARep}
\title{LDA Replications}
\usage{
LDARep(docs, vocab, n = 100, seeds, id = "LDARep", pm.backend, ncpus, ...)
}
\arguments{
\item{docs}{[\code{list}]\cr
Documents as received from \code{\link[tosca]{LDAprep}}.}

\item{vocab}{[\code{character}]\cr
Vocabularies passed to \code{\link[lda]{lda.collapsed.gibbs.sampler}}.
For additional (and necessary) arguments passed, see ellipsis (three-dot argument).}

\item{n}{[\code{integer(1)}]\cr
Number of Replications.}

\item{seeds}{[\code{integer(n)}]\cr
Random Seeds for each Replication.}

\item{id}{[\code{character(1)}]\cr
Name for the computation.}

\item{pm.backend}{[\code{character(1)}]\cr
One of "multicore", "socket" or "mpi".
If \code{pm.backend} is set, \code{\link[parallelMap]{parallelStart}} is
called before computation is started and \code{\link[parallelMap]{parallelStop}}
is called after.}

\item{ncpus}{[\code{integer(1)}]\cr
Number of (physical) CPUs to use. If \code{pm.backend} is passed,
default is determined by \code{\link[future]{availableCores}}.}

\item{...}{additional arguments passed to \code{\link[lda]{lda.collapsed.gibbs.sampler}}.
Arguments will be coerced to a vector of length \code{n}.
Default parameters are \code{alpha = eta = 1/K} and \code{num.iterations = 200}.
There is no default for \code{K}.}
}
\value{
[\code{named list}] with entries \code{id} for computation's name,
\code{jobs} for the parameter settings and \code{lda} for the results itself.
}
\description{
Performs multiple runs of Latent Dirichlet Allocation.
}
\details{
The function generates multiple LDA runs with the possibility of
using parallelization. The integration is done by the
\code{\link[parallelMap:parallelMap]{parallelMap-package}}.

The function returns a \code{LDARep} object. You can receive results and
all other elements of this object with getter functions (see \code{\link{getJob}}).
}
\examples{
res = LDARep(docs = reuters_docs, vocab = reuters_vocab, n = 4, seeds = 1:4,
   id = "myComputation", K = 7:10, alpha = 1, eta = 0.01, num.iterations = 20)
res
getJob(res)
getID(res)
getLDA(res, 4)

\donttest{
LDARep(docs = reuters_docs, vocab = reuters_vocab,
   K = 10, num.iterations = 100, pm.backend = "socket")
}

}
\seealso{
Other replication functions: 
\code{\link{LDAPrototype}()},
\code{\link{as.LDARep}()},
\code{\link{getJob}()},
\code{\link{mergeRepTopics}()}

Other LDA functions: 
\code{\link{LDABatch}()},
\code{\link{LDA}()},
\code{\link{getTopics}()}

Other workflow functions: 
\code{\link{SCLOP}()},
\code{\link{dendTopics}()},
\code{\link{getPrototype}()},
\code{\link{jaccardTopics}()},
\code{\link{mergeTopics}()}
}
\concept{LDA functions}
\concept{replication functions}
\concept{workflow functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as.LDARep.R
\name{as.LDARep}
\alias{as.LDARep}
\alias{as.LDARep.default}
\alias{as.LDARep.LDARep}
\alias{is.LDARep}
\title{LDARep Constructor}
\usage{
as.LDARep(...)

\method{as.LDARep}{default}(lda, job, id, ...)

\method{as.LDARep}{LDARep}(x, ...)

is.LDARep(obj, verbose = FALSE)
}
\arguments{
\item{...}{additional arguments}

\item{lda}{[\code{named list}]\cr
List of \code{\link{LDA}} objects, named by the corresponding "job.id" (\code{integerish}).
If list is unnamed, names are set.}

\item{job}{[\code{\link{data.frame}} or \code{named vector}]\cr
A data.frame or data.table with named columns (at least)
"job.id" (\code{integerish}), "K", "alpha", "eta" and "num.iterations"
or a named vector with entries (at least) "K", "alpha", "eta" and "num.iterations".
If not passed, it is interpreted from \code{param} of each LDA.}

\item{id}{[\code{character(1)}]\cr
A name for the computation. If not passed, it is set to "LDARep".}

\item{x}{[\code{named list}]\cr
\code{\link{LDABatch}} or \code{\link{LDARep}} object.}

\item{obj}{[\code{R} object]\cr
Object to test.}

\item{verbose}{[\code{logical(1)}]\cr
Should test information be given in the console?}
}
\value{
[\code{named list}] with entries \code{id} for computation's name,
\code{jobs} for the parameter settings and \code{lda} for the results themselves.
}
\description{
Constructs a \code{\link{LDARep}} object for given elements \code{lda},
\code{job} and \code{id}.
}
\details{
Given a list of \code{\link{LDA}} objects the function returns
a \code{\link{LDARep}} object, which can be handled using the getter functions
at \code{\link{getJob}}.
}
\examples{
res = LDARep(docs = reuters_docs, vocab = reuters_vocab, n = 4, K = 7, num.iterations = 20)
lda = getLDA(res)

res2 = as.LDARep(lda, id = "newName")
res2
getJob(res2)
getJob(res)

\dontrun{
batch = LDABatch(docs = reuters_docs, vocab = reuters_vocab, n = 4, id = "TEMP", K = 30)
res3 = as.LDARep(batch)
res3
getJob(res3)
}

}
\seealso{
Other constructor functions: 
\code{\link{LDA}()},
\code{\link{as.LDABatch}()}

Other replication functions: 
\code{\link{LDAPrototype}()},
\code{\link{LDARep}()},
\code{\link{getJob}()},
\code{\link{mergeRepTopics}()}
}
\concept{constructor functions}
\concept{replication functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rboTopics.R
\name{rboTopics}
\alias{rboTopics}
\title{Pairwise RBO Similarities}
\usage{
rboTopics(topics, k, p, progress = TRUE, pm.backend, ncpus)
}
\arguments{
\item{topics}{[\code{named matrix}]\cr
The counts of vocabularies/words (row wise) in topics (column wise).}

\item{k}{[\code{integer(1)}]\cr
Maximum depth for evaluation. Words down to this rank are considered for the calculation of similarities.}

\item{p}{[0,1]\cr
Weighting parameter. Lower values emphasizes top ranked words while values
that go towards 1 correspond to equal weights for each evaluation depth.}

\item{progress}{[\code{logical(1)}]\cr
Should a nice progress bar be shown? Turning it off, could lead to significantly
faster calculation. Default is \code{TRUE}.
If \code{pm.backend} is set, parallelization is done and no progress bar will be shown.}

\item{pm.backend}{[\code{character(1)}]\cr
One of "multicore", "socket" or "mpi".
If \code{pm.backend} is set, \code{\link[parallelMap]{parallelStart}} is
called before computation is started and \code{\link[parallelMap]{parallelStop}}
is called after.}

\item{ncpus}{[\code{integer(1)}]\cr
Number of (physical) CPUs to use. If \code{pm.backend} is passed,
default is determined by \code{\link[future]{availableCores}}.}
}
\value{
[\code{named list}] with entries
\describe{
  \item{\code{sims}}{[\code{lower triangular named matrix}] with all pairwise
  similarities of the given topics.}
  \item{\code{wordslimit}}{[\code{integer}] = vocabulary size. See
  \code{\link{jaccardTopics}} for original purpose.}
  \item{\code{wordsconsidered}}{[\code{integer}] = vocabulary size. See
  \code{\link{jaccardTopics}} for original purpose.}
  \item{\code{param}}{[\code{named list}] with parameter
  \code{type} [\code{character(1)}] \code{= "RBO Similarity"},
  \code{k} [\code{integer(1)}] and \code{p} [0,1]. See above for explanation.}
}
}
\description{
Calculates the similarity of all pairwise topic combinations using the
rank-biased overlap (RBO) Similarity.
}
\details{
The RBO Similarity for two topics \eqn{\bm z_{i}} and \eqn{\bm z_{j}}
is calculated by
\deqn{RBO(\bm z_{i}, \bm z_{j} \mid k, p) = 2p^k\frac{\left|Z_{i}^{(k)} \cap Z_{j}^{(k)}\right|}{\left|Z_{i}^{(k)}\right| + \left|Z_{j}^{(k)}\right|} + \frac{1-p}{p} \sum_{d=1}^k 2 p^d\frac{\left|Z_{i}^{(d)} \cap Z_{j}^{(d)}\right|}{\left|Z_{i}^{(d)}\right| + \left|Z_{j}^{(d)}\right|}}
with \eqn{Z_{i}^{(d)}} is the vocabulary set of topic \eqn{\bm z_{i}} down to
rank \eqn{d}. Ties in ranks are resolved by taking the minimum.

The value \code{wordsconsidered} describes the number of words per topic
ranked at rank \eqn{k} or above.
}
\examples{
res = LDARep(docs = reuters_docs, vocab = reuters_vocab, n = 4, K = 10, num.iterations = 30)
topics = mergeTopics(res, vocab = reuters_vocab)
rbo = rboTopics(topics, k = 12, p = 0.9)
rbo

sim = getSimilarity(rbo)
dim(sim)

}
\references{
Webber, William, Alistair Moffat and Justin Zobel (2010).
"A similarity measure for indefinite rankings".
In: \emph{ACM Transations on Information Systems} 28(4), p.20:1–-20:38,
DOI 10.1145/1852102.1852106,
URL \url{https://doi.acm.org/10.1145/1852102.1852106}
}
\seealso{
Other TopicSimilarity functions: 
\code{\link{cosineTopics}()},
\code{\link{dendTopics}()},
\code{\link{getSimilarity}()},
\code{\link{jaccardTopics}()},
\code{\link{jsTopics}()}
}
\concept{TopicSimilarity functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ldaPrototype-package.R
\docType{package}
\name{ldaPrototype-package}
\alias{ldaPrototype}
\alias{ldaPrototype-package}
\title{ldaPrototype: Prototype of Multiple Latent Dirichlet Allocation Runs}
\description{
Determine a Prototype from a number of runs of Latent Dirichlet
Allocation (LDA) measuring its similarities with S-CLOP: A procedure to select
the LDA run with highest mean pairwise similarity, which is measured by S-CLOP
(Similarity of multiple sets by Clustering with Local Pruning), to all other
runs. LDA runs are specified by its assignments leading to estimators for
distribution parameters. Repeated runs lead to different results, which we
encounter by choosing the most representative LDA run as prototype.\cr
For bug reports and feature requests please use the issue tracker:
\url{https://github.com/JonasRieger/ldaPrototype/issues}. Also have a look at
the (detailed) example at \url{https://github.com/JonasRieger/ldaPrototype}.
}
\section{Data}{

\code{\link{reuters}} Example Dataset (91 articles from Reuters) for testing.
}

\section{Constructor}{

\code{\link{LDA}} LDA objects used in this package.\cr
\code{\link{as.LDARep}} LDARep objects.\cr
\code{\link{as.LDABatch}} LDABatch objects.
}

\section{Getter}{

\code{\link{getTopics}} Getter for \code{\link{LDA}} objects.\cr
\code{\link{getJob}} Getter for \code{\link{LDARep}} and \code{\link{LDABatch}} objects.\cr
\code{\link{getSimilarity}} Getter for \code{\link[=jaccardTopics]{TopicSimilarity}} objects.\cr
\code{\link{getSCLOP}} Getter for \code{\link[=getPrototype]{PrototypeLDA}} objects.\cr
\code{\link{getPrototype}} Determine the Prototype LDA.
}

\section{Performing multiple LDAs}{

\code{\link{LDARep}} Performing multiple LDAs locally (using parallelization).\cr
\code{\link{LDABatch}} Performing multiple LDAs on Batch Systems.
}

\section{Calculation Steps (Workflow) to determine the Prototype LDA}{

\code{\link{mergeTopics}} Merge topic matrices from multiple LDAs.\cr
\code{\link{jaccardTopics}} Calculate topic similarities using the Jaccard coefficient (see Similarity Measures for other possible measures).\cr
\code{\link{dendTopics}} Create a dendrogram from topic similarities.\cr
\code{\link{SCLOP}} Determine various S-CLOP values.\cr
\code{\link{pruneSCLOP}} Prune \code{\link[=dendTopics]{TopicDendrogram}} objects.
}

\section{Similarity Measures}{

\code{\link{cosineTopics}} Cosine Similarity.\cr
\code{\link{jaccardTopics}} Jaccard Coefficient.\cr
\code{\link{jsTopics}} Jensen-Shannon Divergence.\cr
\code{\link{rboTopics}} rank-biased overlap.
}

\section{Shortcuts}{

\code{\link{getPrototype}} Shortcut which includes all calculation steps.\cr
\code{\link{LDAPrototype}} Shortcut which performs multiple LDAs and
determines their Prototype.
}

\references{
Rieger, Jonas (2020). "ldaPrototype: A method in R to get a Prototype of multiple Latent
Dirichlet Allocations". Journal of Open Source Software, \bold{5}(51), 2181,
\doi{10.21105/joss.02181}.

Rieger, Jonas, Jörg Rahnenführer and Carsten Jentsch (2020).
"Improving Latent Dirichlet Allocation: On Reliability of the Novel Method LDAPrototype".
In: \emph{Natural Language Processing and Information Systems, NLDB 2020.} LNCS 12089, pp. 118--125,
\doi{10.1007/978-3-030-51310-8_11}.

Rieger, Jonas, Lars Koppers, Carsten Jentsch and Jörg Rahnenführer (2020).
"Improving Reliability of Latent Dirichlet Allocation by Assessing Its Stability using Clustering Techniques on Replicated Runs".
arXiv 2003.04980, URL \url{https://arxiv.org/abs/2003.04980}.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/JonasRieger/ldaPrototype}
  \item Report bugs at \url{https://github.com/JonasRieger/ldaPrototype/issues}
}

}
\author{
\strong{Maintainer}: Jonas Rieger \email{jonas.rieger@tu-dortmund.de} (\href{https://orcid.org/0000-0002-0007-4478}{ORCID})

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SCLOP.R
\name{SCLOP}
\alias{SCLOP}
\alias{disparitySum}
\alias{SCLOP.pairwise}
\title{Similarity/Stability of multiple sets of Objects using Clustering with Local Pruning}
\usage{
SCLOP(dend)

disparitySum(dend)

SCLOP.pairwise(sims)
}
\arguments{
\item{dend}{[\code{\link[stats]{dendrogram}}]\cr
Output from \code{\link{dendTopics}}.}

\item{sims}{[\code{\link[=jaccardTopics]{TopicSimilarity}} object
or \code{lower triangular named matrix}]\cr
\code{\link[=jaccardTopics]{TopicSimilarity}} object or
pairwise jaccard similarities of underlying topics as the \code{sims} element
from \code{\link[=jaccardTopics]{TopicSimilarity}} objects. The topic names should be
formatted as <\emph{Run X}>.<\emph{Topic Y}>, so that the name before the
first dot identifies the LDA run.}
}
\value{
\describe{
  \item{\code{SCLOP}}{[0,1] value specifying the S-CLOP for the best possible
  local pruning state of the given dendrogram.}
  \item{\code{disparitySum}}{[\code{numeric(1)}] value specifying the least
  possible sum of disparities on the given dendrogram.}
  \item{\code{SCLOP.pairwise}}{[\code{symmetrical named matrix}] with all
  pairwise S-CLOP scores of the given LDA runs.}
}
}
\description{
The function \code{SCLOP} calculates the S-CLOP value for the best possible
local pruning state of a dendrogram from \code{\link{dendTopics}}.
The function \code{\link{pruneSCLOP}} supplies the corresponding pruning state itself.\cr
To get all pairwise S-CLOP scores of two LDA runs, the function \code{SCLOP.pairwise}
can be used. It returns a matrix of the pairwise S-CLOP scores.\cr
All three functions use the function \code{disparitySum} to calculate the
least possible sum of disparities (on the best possible local pruning state)
on a given dendrogram.
}
\details{
For one specific cluster \eqn{g} and \eqn{R} LDA Runs the disparity is calculated by
\deqn{U(g) := \frac{1}{R} \sum_{r=1}^R \vert t_r^{(g)} - 1 \vert \cdot \sum_{r=1}^R t_r^{(g)},}
while \eqn{\bm t^{(g)} = (t_1^{(g)}, ..., t_R^{(g)})^T}
contains the number of topics that belong to the different LDA runs and that
occur in cluster \eqn{g}.

The function \code{disparitySum} returns the least possible sum of disparities
\eqn{U_{\Sigma}(G^*)} for the best possible pruning state \eqn{G^*}
with \eqn{U_{\Sigma}(G) = \sum_{g \in G} U(g) \to \min}.
The highest possible value for \eqn{U_{\Sigma}(G^*)} is limited by
\deqn{U_{\Sigma,\textsf{max}} := \sum_{g \in \tilde{G}} U(g) = N \cdot \frac{R-1}{R},}
with \eqn{\tilde{G}} denotes the corresponding worst case pruning state. This worst
case scenario is useful for normalizing the SCLOP scores.

The function \code{SCLOP} then calculates the value
\deqn{\textsf{S-CLOP}(G^*) := 1 - \frac{1}{U_{\Sigma,\textsf{max}}} \cdot \sum_{g \in G^*} U(g) ~\in [0,1],}
where \eqn{\sum\limits_{g \in G^*} U(g) = U_{\Sigma}(G^*)}.
}
\examples{
res = LDARep(docs = reuters_docs, vocab = reuters_vocab, n = 4, K = 10, num.iterations = 30)
topics = mergeTopics(res, vocab = reuters_vocab)
jacc = jaccardTopics(topics, atLeast = 2)
dend = dendTopics(jacc)

SCLOP(dend)
disparitySum(dend)

SCLOP.pairwise(jacc)
SCLOP.pairwise(getSimilarity(jacc))

}
\seealso{
Other SCLOP functions: 
\code{\link{pruneSCLOP}()}

Other workflow functions: 
\code{\link{LDARep}()},
\code{\link{dendTopics}()},
\code{\link{getPrototype}()},
\code{\link{jaccardTopics}()},
\code{\link{mergeTopics}()}
}
\concept{SCLOP functions}
\concept{workflow functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mergeTopics.R
\name{mergeTopics}
\alias{mergeTopics}
\title{Merge LDA Topic Matrices}
\usage{
mergeTopics(x, vocab, progress = TRUE)
}
\arguments{
\item{x}{[\code{named list}]\cr
\code{\link{LDARep}} or \code{\link{LDABatch}} object.}

\item{vocab}{[\code{character}]\cr
Vocabularies taken into consideration for merging topic matrices.}

\item{progress}{[\code{logical(1)}]\cr
Should a nice progress bar be shown? Turning it off, could lead to significantly
faster calculation. Default ist \code{TRUE}.}
}
\value{
[\code{named matrix}] with the count of vocabularies (row wise) in topics (column wise).
}
\description{
Generic function, which collects LDA results and merges their topic matrices
for a given set of vocabularies.
}
\details{
This function uses the function \code{\link{mergeRepTopics}} or
\code{\link{mergeBatchTopics}}. The topic matrices are transponed and cbinded,
so that the resulting matrix contains the counts of vocabularies/words (row wise)
in topics (column wise).
}
\examples{
res = LDARep(docs = reuters_docs, vocab = reuters_vocab, n = 4, K = 10, num.iterations = 30)
topics = mergeTopics(res, vocab = reuters_vocab)
dim(topics)
length(reuters_vocab)

\dontrun{
res = LDABatch(docs = reuters_docs, vocab = reuters_vocab, n = 4, K = 10, num.iterations = 30)
topics = mergeTopics(res, vocab = reuters_vocab)
dim(topics)
length(reuters_vocab)
}
}
\seealso{
Other merge functions: 
\code{\link{mergeBatchTopics}()},
\code{\link{mergeRepTopics}()}

Other workflow functions: 
\code{\link{LDARep}()},
\code{\link{SCLOP}()},
\code{\link{dendTopics}()},
\code{\link{getPrototype}()},
\code{\link{jaccardTopics}()}
}
\concept{merge functions}
\concept{workflow functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mergeBatchTopics.R
\name{mergeBatchTopics}
\alias{mergeBatchTopics}
\alias{mergeBatchTopics.LDABatch}
\alias{mergeBatchTopics.default}
\title{Merge LDA Topic Matrices}
\usage{
mergeBatchTopics(...)

\method{mergeBatchTopics}{LDABatch}(x, vocab, progress = TRUE, ...)

\method{mergeBatchTopics}{default}(vocab, reg, job, id, progress = TRUE, ...)
}
\arguments{
\item{...}{additional arguments}

\item{x}{[\code{named list}]\cr
\code{\link{LDABatch}} object. Alternatively \code{job}, \code{reg} and
\code{id} can be passed or their defaults are taken.}

\item{vocab}{[\code{character}]\cr
Vocabularies taken into consideration for merging topic matrices. Default is
the vocabulary of the first LDA.}

\item{progress}{[\code{logical(1)}]\cr
Should a nice progress bar be shown? Turning it off, could lead to significantly
faster calculation. Default ist \code{TRUE}.}

\item{reg}{[\code{\link[batchtools:makeRegistry]{Registry}}]\cr
Registry. See \code{\link[batchtools]{reduceResultsList}}.}

\item{job}{[\code{\link{data.frame}} or \code{integer}]\cr
A data.frame or data.table with a column named "job.id" or a vector of integerish job ids.
See \code{\link[batchtools]{reduceResultsList}}.}

\item{id}{[\code{character(1)}]\cr
A name for the registry. If not passed, the folder's name is extracted from \code{reg}.}
}
\value{
[\code{named matrix}] with the count of vocabularies (row wise) in topics (column wise).
}
\description{
Collects LDA results from a given registry and merges their topic matrices for
a given set of vocabularies.
}
\details{
For details and examples see \code{\link{mergeTopics}}.
}
\seealso{
Other merge functions: 
\code{\link{mergeRepTopics}()},
\code{\link{mergeTopics}()}

Other batch functions: 
\code{\link{LDABatch}()},
\code{\link{as.LDABatch}()},
\code{\link{getJob}()}
}
\concept{batch functions}
\concept{merge functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getSCLOP.R
\name{getSCLOP}
\alias{getSCLOP}
\alias{getSimilarity.PrototypeLDA}
\alias{getRelevantWords.PrototypeLDA}
\alias{getConsideredWords.PrototypeLDA}
\alias{getMergedTopics}
\alias{getPrototypeID}
\alias{getLDA.PrototypeLDA}
\alias{getID.PrototypeLDA}
\alias{getParam.PrototypeLDA}
\alias{getJob.PrototypeLDA}
\title{Getter for PrototypeLDA}
\usage{
getSCLOP(x)

\method{getSimilarity}{PrototypeLDA}(x)

\method{getRelevantWords}{PrototypeLDA}(x)

\method{getConsideredWords}{PrototypeLDA}(x)

getMergedTopics(x)

getPrototypeID(x)

\method{getLDA}{PrototypeLDA}(x, job, reduce = TRUE, all = FALSE)

\method{getID}{PrototypeLDA}(x)

\method{getParam}{PrototypeLDA}(x)

\method{getJob}{PrototypeLDA}(x)
}
\arguments{
\item{x}{[\code{named list}]\cr
\code{\link[=getPrototype]{PrototypeLDA}} object.}

\item{job}{[\code{\link{data.frame}} or \code{integer}]\cr
A data.frame or data.table with a column named "job.id" or a vector of
integerish job ids. Default is the (integerish) ID of the Prototype LDA.}

\item{reduce}{[\code{logical(1)}]\cr
If the list of LDAs contains only one element, should the list be reduced and
the single (unnamed) element be returned? Default is \code{TRUE}.
Not considered, if \code{all} is \code{TRUE}.}

\item{all}{[\code{logical(1)}]\cr
Shortcut for \code{job}: Should all stored LDAs be returned?}
}
\description{
Returns the corresponding element of a \code{\link[=getPrototype]{PrototypeLDA}} object.
}
\seealso{
Other getter functions: 
\code{\link{getJob}()},
\code{\link{getSimilarity}()},
\code{\link{getTopics}()}

Other PrototypeLDA functions: 
\code{\link{LDAPrototype}()},
\code{\link{getPrototype}()}
}
\concept{PrototypeLDA functions}
\concept{getter functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LDAPrototype.R
\name{LDAPrototype}
\alias{LDAPrototype}
\title{Determine the Prototype LDA}
\usage{
LDAPrototype(
  docs,
  vocabLDA,
  vocabMerge = vocabLDA,
  n = 100,
  seeds,
  id = "LDARep",
  pm.backend,
  ncpus,
  limit.rel,
  limit.abs,
  atLeast,
  progress = TRUE,
  keepTopics = FALSE,
  keepSims = FALSE,
  keepLDAs = FALSE,
  ...
)
}
\arguments{
\item{docs}{[\code{list}]\cr
Documents as received from \code{\link[tosca]{LDAprep}}.}

\item{vocabLDA}{[\code{character}]\cr
Vocabularies passed to \code{\link[lda]{lda.collapsed.gibbs.sampler}}.
For additional (and necessary) arguments passed, see ellipsis (three-dot argument).}

\item{vocabMerge}{[\code{character}]\cr
Vocabularies taken into consideration for merging topic matrices.}

\item{n}{[\code{integer(1)}]\cr
Number of Replications.}

\item{seeds}{[\code{integer(n)}]\cr
Random Seeds for each Replication.}

\item{id}{[\code{character(1)}]\cr
Name for the computation.}

\item{pm.backend}{[\code{character(1)}]\cr
One of "multicore", "socket" or "mpi".
If \code{pm.backend} is set, \code{\link[parallelMap]{parallelStart}} is
called before computation is started and \code{\link[parallelMap]{parallelStop}}
is called after.}

\item{ncpus}{[\code{integer(1)}]\cr
Number of (physical) CPUs to use. If \code{pm.backend} is passed,
default is determined by \code{\link[future]{availableCores}}.}

\item{limit.rel}{[0,1]\cr
See \code{\link{jaccardTopics}}. Default is \code{1/500}.}

\item{limit.abs}{[\code{integer(1)}]\cr
See \code{\link{jaccardTopics}}. Default is \code{10}.}

\item{atLeast}{[\code{integer(1)}]\cr
See \code{\link{jaccardTopics}}. Default is \code{0}.}

\item{progress}{[\code{logical(1)}]\cr
Should a nice progress bar be shown for the steps of \code{\link{mergeTopics}}
and \code{\link{jaccardTopics}}? Turning it off, could lead to significantly
faster calculation. Default ist \code{TRUE}.}

\item{keepTopics}{[\code{logical(1)}]\cr
Should the merged topic matrix from \code{\link{mergeTopics}} be kept?}

\item{keepSims}{[\code{logical(1)}]\cr
Should the calculated topic similarities matrix from \code{\link{jaccardTopics}} be kept?}

\item{keepLDAs}{[\code{logical(1)}]\cr
Should the considered LDAs be kept?}

\item{...}{additional arguments passed to \code{\link[lda]{lda.collapsed.gibbs.sampler}}.
Arguments will be coerced to a vector of length \code{n}.
Default parameters are \code{alpha = eta = 1/K} and \code{num.iterations = 200}.
There is no default for \code{K}.}
}
\value{
[\code{named list}] with entries
 \describe{
  \item{\code{id}}{[\code{character(1)}] See above.}
  \item{\code{protoid}}{[\code{character(1)}] Name (ID) of the determined Prototype LDA.}
  \item{\code{lda}}{List of \code{\link{LDA}} objects of the determined Prototype LDA
  and - if \code{keepLDAs} is \code{TRUE} - all considered LDAs.}
  \item{\code{jobs}}{[\code{data.table}] with parameter specifications for the LDAs.}
  \item{\code{param}}{[\code{named list}] with parameter specifications for
  \code{limit.rel} [0,1], \code{limit.abs} [\code{integer(1)}] and
  \code{atLeast} [\code{integer(1)}]. See above for explanation.}
  \item{\code{topics}}{[\code{named matrix}] with the count of vocabularies
  (row wise) in topics (column wise).}
  \item{\code{sims}}{[\code{lower triangular named matrix}] with all pairwise
  jaccard similarities of the given topics.}
  \item{\code{wordslimit}}{[\code{integer}] with counts of words determined as
  relevant based on \code{limit.rel} and \code{limit.abs}.}
  \item{\code{wordsconsidered}}{[\code{integer}] with counts of considered
  words for similarity calculation. Could differ from \code{wordslimit}, if
  \code{atLeast} is greater than zero.}
  \item{\code{sclop}}{[\code{symmetrical named matrix}] with all pairwise
  S-CLOP scores of the given LDA runs.}
}
}
\description{
Performs multiple runs of LDA and computes the Prototype LDA of
this set of LDAs.
}
\details{
While \code{LDAPrototype} marks the overall shortcut for performing
multiple LDA runs and choosing the Prototype of them, \code{\link{getPrototype}}
just hooks up at determining the Prototype. The generation of multiple LDAs
has to be done before use of \code{getPrototype}.

To save memory a lot of interim calculations are discarded by default.

If you use parallel computation, no progress bar is shown.

For details see the details sections of the workflow functions at \code{\link{getPrototype}}.
}
\examples{
res = LDAPrototype(docs = reuters_docs, vocabLDA = reuters_vocab,
   n = 4, K = 10, num.iterations = 30)
res
getPrototype(res) # = getLDA(res)
getSCLOP(res)

res = LDAPrototype(docs = reuters_docs, vocabLDA = reuters_vocab,
   n = 4, K = 10, num.iterations = 30, keepLDAs = TRUE)
res
getLDA(res, all = TRUE)
getPrototypeID(res)
getParam(res)

}
\seealso{
Other shortcut functions: 
\code{\link{getPrototype}()}

Other PrototypeLDA functions: 
\code{\link{getPrototype}()},
\code{\link{getSCLOP}()}

Other replication functions: 
\code{\link{LDARep}()},
\code{\link{as.LDARep}()},
\code{\link{getJob}()},
\code{\link{mergeRepTopics}()}
}
\concept{PrototypeLDA functions}
\concept{replication functions}
\concept{shortcut functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LDA.R
\name{LDA}
\alias{LDA}
\alias{as.LDA}
\alias{is.LDA}
\title{LDA Object}
\usage{
LDA(
  x,
  param,
  assignments,
  topics,
  document_sums,
  document_expects,
  log.likelihoods
)

as.LDA(
  x,
  param,
  assignments,
  topics,
  document_sums,
  document_expects,
  log.likelihoods
)

is.LDA(obj, verbose = FALSE)
}
\arguments{
\item{x}{[\code{named list}]\cr
Output from \code{\link[lda]{lda.collapsed.gibbs.sampler}}. Alternatively each
element can be passed for individual results. Individually set elements
overwrite elements from \code{x}.}

\item{param}{[\code{named list}]\cr
Parameters of the function call \code{\link[lda]{lda.collapsed.gibbs.sampler}}.
List always should contain names "K", "alpha", "eta" and "num.iterations".}

\item{assignments}{Individual element for LDA object.}

\item{topics}{Individual element for LDA object.}

\item{document_sums}{Individual element for LDA object.}

\item{document_expects}{Individual element for LDA object.}

\item{log.likelihoods}{Individual element for LDA object.}

\item{obj}{[\code{R} object]\cr
Object to test.}

\item{verbose}{[\code{logical(1)}]\cr
Should test information be given in the console?}
}
\value{
[\code{named list}] LDA object.
}
\description{
Constructor for LDA objects used in this package.
}
\details{
The functions \code{LDA} and \code{as.LDA} do exactly the same. If you call
\code{LDA} on an object \code{x} which already is of the structure of an
\code{LDA} object (in particular a \code{LDA} object itself),
the additional arguments \code{param, assignments, ...}
may be used to override the specific elements.
}
\examples{
res = LDARep(docs = reuters_docs, vocab = reuters_vocab, n = 1, K = 10)
lda = getLDA(res)

LDA(lda)
# does not change anything

LDA(lda, assignments = NULL)
# creates a new LDA object without the assignments element

LDA(param = getParam(lda), topics = getTopics(lda))
# creates a new LDA object with elements param and topics

}
\seealso{
Other constructor functions: 
\code{\link{as.LDABatch}()},
\code{\link{as.LDARep}()}

Other LDA functions: 
\code{\link{LDABatch}()},
\code{\link{LDARep}()},
\code{\link{getTopics}()}
}
\concept{LDA functions}
\concept{constructor functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jsTopics.R
\name{jsTopics}
\alias{jsTopics}
\title{Pairwise Jensen-Shannon Similarities (Divergences)}
\usage{
jsTopics(topics, epsilon = 1e-06, progress = TRUE, pm.backend, ncpus)
}
\arguments{
\item{topics}{[\code{named matrix}]\cr
The counts of vocabularies/words (row wise) in topics (column wise).}

\item{epsilon}{[\code{numeric(1)}]\cr
Numerical value added to \code{topics} to ensure computability. See details.
Default is \code{1e-06}.}

\item{progress}{[\code{logical(1)}]\cr
Should a nice progress bar be shown? Turning it off, could lead to significantly
faster calculation. Default is \code{TRUE}.
If \code{pm.backend} is set, parallelization is done and no progress bar will be shown.}

\item{pm.backend}{[\code{character(1)}]\cr
One of "multicore", "socket" or "mpi".
If \code{pm.backend} is set, \code{\link[parallelMap]{parallelStart}} is
called before computation is started and \code{\link[parallelMap]{parallelStop}}
is called after.}

\item{ncpus}{[\code{integer(1)}]\cr
Number of (physical) CPUs to use. If \code{pm.backend} is passed,
default is determined by \code{\link[future]{availableCores}}.}
}
\value{
[\code{named list}] with entries
\describe{
  \item{\code{sims}}{[\code{lower triangular named matrix}] with all pairwise
  similarities of the given topics.}
  \item{\code{wordslimit}}{[\code{integer}] = vocabulary size. See
  \code{\link{jaccardTopics}} for original purpose.}
  \item{\code{wordsconsidered}}{[\code{integer}] = vocabulary size. See
  \code{\link{jaccardTopics}} for original purpose.}
  \item{\code{param}}{[\code{named list}] with parameter specifications for
  \code{type} [\code{character(1)}] \code{= "Cosine Similarity"} and
  \code{epsilon} [\code{numeric(1)}]. See above for explanation.}
}
}
\description{
Calculates the similarity of all pairwise topic combinations using the
Jensen-Shannon Divergence.
}
\details{
The Jensen-Shannon Similarity for two topics \eqn{\bm z_{i}} and
\eqn{\bm z_{j}} is calculated by
\deqn{JS(\bm z_{i}, \bm z_{j}) = 1 - \left( KLD\left(\bm p_i, \frac{\bm p_i + \bm p_j}{2}\right) + KLD\left(\bm p_j, \frac{\bm p_i + \bm p_j}{2}\right) \right)/2}
\deqn{= 1 - KLD(\bm p_i, \bm p_i + \bm p_j)/2 - KLD(\bm p_j, \bm p_i + \bm p_j)/2 - \log(2)}
with \eqn{V} is the vocabulary size, \eqn{\bm p_k = \left(p_k^{(1)}, ..., p_k^{(V)}\right)},
and \eqn{p_k^{(v)}} is the proportion of assignments of the
\eqn{v}-th word to the \eqn{k}-th topic. KLD defines the Kullback-Leibler
Divergence calculated by
\deqn{KLD(\bm p_{k}, \bm p_{\Sigma}) = \sum_{v=1}^{V} p_k^{(v)} \log{\frac{p_k^{(v)}}{p_{\Sigma}^{(v)}}}.}

There is an \code{epsilon} added to every \eqn{n_k^{(v)}}, the count
(not proportion) of assignments to ensure computability with respect to zeros.
}
\examples{
res = LDARep(docs = reuters_docs, vocab = reuters_vocab, n = 4, K = 10, num.iterations = 30)
topics = mergeTopics(res, vocab = reuters_vocab)
js = jsTopics(topics)
js

sim = getSimilarity(js)
dim(sim)

js1 = jsTopics(topics, epsilon = 1)
sim1 = getSimilarity(js1)
summary((sim1-sim)[lower.tri(sim)])
plot(sim, sim1, xlab = "epsilon = 1e-6", ylab = "epsilon = 1")

}
\seealso{
Other TopicSimilarity functions: 
\code{\link{cosineTopics}()},
\code{\link{dendTopics}()},
\code{\link{getSimilarity}()},
\code{\link{jaccardTopics}()},
\code{\link{rboTopics}()}
}
\concept{TopicSimilarity functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as.LDABatch.R
\name{as.LDABatch}
\alias{as.LDABatch}
\alias{is.LDABatch}
\title{LDABatch Constructor}
\usage{
as.LDABatch(reg, job, id)

is.LDABatch(obj, verbose = FALSE)
}
\arguments{
\item{reg}{[\code{\link[batchtools:makeRegistry]{Registry}}]\cr
Registry. See \code{\link[batchtools:findJobs]{findDone}}.}

\item{job}{[\code{\link{data.frame}} or \code{integer}]\cr
A data.frame or data.table with a column named "job.id" or a vector of integerish job ids.
See \code{\link[batchtools]{reduceResultsList}}.}

\item{id}{[\code{character(1)}]\cr
A name for the registry. If not passed, the folder's name is extracted from \code{reg}.}

\item{obj}{[\code{R} object]\cr
Object to test.}

\item{verbose}{[\code{logical(1)}]\cr
Should test information be given in the console?}
}
\value{
[\code{named list}] with entries \code{id} for the registry's folder name,
\code{jobs} for the submitted jobs' ids and its parameter settings and
\code{reg} for the registry itself.
}
\description{
Constructs a \code{\link{LDABatch}} object for given elements \code{reg},
\code{job} and \code{id}.
}
\details{
Given a \code{\link[batchtools:makeRegistry]{Registry}} the function returns
a \code{\link{LDABatch}} object, which can be handled using the getter functions
at \code{\link{getJob}}.
}
\examples{
\dontrun{
batch = LDABatch(docs = reuters_docs, vocab = reuters_vocab, K = 15, chunk.size = 20)
batch

batch2 = as.LDABatch(reg = getRegistry(batch))
batch2
head(getJob(batch2))

batch3 = as.LDABatch()
batch3

### one way of loading an existing registry ###
batchtools::loadRegistry("LDABatch")
batch = as.LDABatch()
}

}
\seealso{
Other constructor functions: 
\code{\link{LDA}()},
\code{\link{as.LDARep}()}

Other batch functions: 
\code{\link{LDABatch}()},
\code{\link{getJob}()},
\code{\link{mergeBatchTopics}()}
}
\concept{batch functions}
\concept{constructor functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dendTopics.R
\name{dendTopics}
\alias{dendTopics}
\alias{plot.TopicDendrogram}
\title{Topic Dendrogram}
\usage{
dendTopics(sims, ind, method = "complete")

\method{plot}{TopicDendrogram}(x, pruning, pruning.par, ...)
}
\arguments{
\item{sims}{[\code{\link[=jaccardTopics]{TopicSimilarity}} object
or \code{lower triangular named matrix}]\cr
\code{\link[=jaccardTopics]{TopicSimilarity}} object or
pairwise jaccard similarities of underlying topics as the \code{sims} element
from \code{\link[=jaccardTopics]{TopicSimilarity}} objects. The topic names should be
formatted as <\emph{Run X}>.<\emph{Topic Y}>, so that the name before the
first dot identifies the LDA run.}

\item{ind}{[\code{integer}, \code{logical} or \code{character}]\cr
An integerish vector (or logical of the same length as the number of rows and columns)
for specifying the topics taken into account. Alternatively
a character vector can be passed. Then, all topics are taken for which the name
contain at least one of the phrases in \code{ind} (see \code{\link[=grep]{grepl}}).
By default all topics are considered.}

\item{method}{[\code{character(1)}]\cr
The agglomeration method. See \code{\link[stats]{hclust}}.}

\item{x}{an R object.}

\item{pruning}{[\code{list of \link[stats]{dendrogram}s}]\cr
\code{\link[=pruneSCLOP]{PruningSCLOP}} object specifying the best possible
local pruning state.}

\item{pruning.par}{[\code{list}]\cr
List of parameters to mark the pruning. See section "Details" at \code{\link{dendTopics}}
for default parameters. Types for marking the pruning state are \code{"abline"},
\code{"color"} and \code{"both"}.}

\item{...}{additional arguments.}
}
\value{
[\code{\link[stats]{dendrogram}}] \code{\link[=dendTopics]{TopicDendrogram}} object
(and \code{\link[stats]{dendrogram}} object) of all considered topics.
}
\description{
Builds a dendrogram for topics based on their pairwise similarities using the
cluster algorithm \code{\link[stats]{hclust}}.
}
\details{
The label´s colors are determined based on their Run belonging using
\code{\link[colorspace]{rainbow_hcl}} by default. Colors can be manipulated
using \code{\link[dendextend]{labels_colors}}. Analogously, the labels
themself can be  manipulated using \code{\link[dendextend:labels.hclust]{labels}}.
For both the function \code{\link[stats]{order.dendrogram}} is useful.

The resulting \code{\link[stats]{dendrogram}} can be plotted. In addition,
it is possible to mark a pruning state in the plot, either by color or by
separator lines (or both) setting \code{pruning.par}. For the default values
of \code{pruning.par} call the corresponding function on any
\code{\link[=pruneSCLOP]{PruningSCLOP}} object.
}
\examples{
res = LDARep(docs = reuters_docs, vocab = reuters_vocab, n = 4, K = 10, num.iterations = 30)
topics = mergeTopics(res, vocab = reuters_vocab)
jacc = jaccardTopics(topics, atLeast = 2)
sim = getSimilarity(jacc)

dend = dendTopics(jacc)
dend2 = dendTopics(sim)

\donttest{
plot(dend)
plot(dendTopics(jacc, ind = c("Rep2", "Rep3")))
}

pruned = pruneSCLOP(dend)
\donttest{
plot(dend, pruning = pruned)
plot(dend, pruning = pruned, pruning.par = list(type = "color"))
plot(dend, pruning = pruned, pruning.par = list(type = "both", lty = 1, lwd = 2, col = "red"))

dend2 = dendTopics(jacc, ind = c("Rep2", "Rep3"))
plot(dend2, pruning = pruneSCLOP(dend2), pruning.par = list(lwd = 2, col = "darkgrey"))
}

}
\seealso{
Other plot functions: 
\code{\link{pruneSCLOP}()}

Other TopicSimilarity functions: 
\code{\link{cosineTopics}()},
\code{\link{getSimilarity}()},
\code{\link{jaccardTopics}()},
\code{\link{jsTopics}()},
\code{\link{rboTopics}()}

Other workflow functions: 
\code{\link{LDARep}()},
\code{\link{SCLOP}()},
\code{\link{getPrototype}()},
\code{\link{jaccardTopics}()},
\code{\link{mergeTopics}()}
}
\concept{TopicSimilarity functions}
\concept{plot functions}
\concept{workflow functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cosineTopics.R
\name{cosineTopics}
\alias{cosineTopics}
\title{Pairwise Cosine Similarities}
\usage{
cosineTopics(topics, progress = TRUE, pm.backend, ncpus)
}
\arguments{
\item{topics}{[\code{named matrix}]\cr
The counts of vocabularies/words (row wise) in topics (column wise).}

\item{progress}{[\code{logical(1)}]\cr
Should a nice progress bar be shown? Turning it off, could lead to significantly
faster calculation. Default is \code{TRUE}.
If \code{pm.backend} is set, parallelization is done and no progress bar will be shown.}

\item{pm.backend}{[\code{character(1)}]\cr
One of "multicore", "socket" or "mpi".
If \code{pm.backend} is set, \code{\link[parallelMap]{parallelStart}} is
called before computation is started and \code{\link[parallelMap]{parallelStop}}
is called after.}

\item{ncpus}{[\code{integer(1)}]\cr
Number of (physical) CPUs to use. If \code{pm.backend} is passed,
default is determined by \code{\link[future]{availableCores}}.}
}
\value{
[\code{named list}] with entries
\describe{
  \item{\code{sims}}{[\code{lower triangular named matrix}] with all pairwise
  similarities of the given topics.}
  \item{\code{wordslimit}}{[\code{integer}] = vocabulary size. See
  \code{\link{jaccardTopics}} for original purpose.}
  \item{\code{wordsconsidered}}{[\code{integer}] = vocabulary size. See
  \code{\link{jaccardTopics}} for original purpose.}
  \item{\code{param}}{[\code{named list}] with parameter
  \code{type} [\code{character(1)}] \code{= "Cosine Similarity"}.}
}
}
\description{
Calculates the similarity of all pairwise topic combinations using the
Cosine Similarity.
}
\details{
The Cosine Similarity for two topics \eqn{\bm z_{i}} and \eqn{\bm z_{j}}
is calculated by
\deqn{ \cos(\theta | \bm z_{i}, \bm z_{j}) = \frac{ \sum_{v=1}^{V}{n_{i}^{(v)} n_{j}^{(v)}} }{ \sqrt{\sum_{v=1}^{V}{\left(n_{i}^{(v)}\right)^2}} \sqrt{\sum_{v=1}^{V}{\left(n_{j}^{(v)}\right)^2}} }}
with \eqn{\theta} determining the angle between the corresponding
count vectors \eqn{\bm z_{i}} and \eqn{\bm z_{j}},
\eqn{V} is the vocabulary size and \eqn{n_k^{(v)}} is the count of
assignments of the \eqn{v}-th word to the \eqn{k}-th topic.
}
\examples{
res = LDARep(docs = reuters_docs, vocab = reuters_vocab, n = 4, K = 10, num.iterations = 30)
topics = mergeTopics(res, vocab = reuters_vocab)
cosine = cosineTopics(topics)
cosine

sim = getSimilarity(cosine)
dim(sim)

}
\seealso{
Other TopicSimilarity functions: 
\code{\link{dendTopics}()},
\code{\link{getSimilarity}()},
\code{\link{jaccardTopics}()},
\code{\link{jsTopics}()},
\code{\link{rboTopics}()}
}
\concept{TopicSimilarity functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getPrototype.R
\name{getPrototype}
\alias{getPrototype}
\alias{getPrototype.LDARep}
\alias{getPrototype.LDABatch}
\alias{getPrototype.default}
\title{Determine the Prototype LDA}
\usage{
getPrototype(...)

\method{getPrototype}{LDARep}(
  x,
  vocab,
  limit.rel,
  limit.abs,
  atLeast,
  progress = TRUE,
  pm.backend,
  ncpus,
  keepTopics = FALSE,
  keepSims = FALSE,
  keepLDAs = FALSE,
  sclop,
  ...
)

\method{getPrototype}{LDABatch}(
  x,
  vocab,
  limit.rel,
  limit.abs,
  atLeast,
  progress = TRUE,
  pm.backend,
  ncpus,
  keepTopics = FALSE,
  keepSims = FALSE,
  keepLDAs = FALSE,
  sclop,
  ...
)

\method{getPrototype}{default}(
  lda,
  vocab,
  id,
  job,
  limit.rel,
  limit.abs,
  atLeast,
  progress = TRUE,
  pm.backend,
  ncpus,
  keepTopics = FALSE,
  keepSims = FALSE,
  keepLDAs = FALSE,
  sclop,
  ...
)
}
\arguments{
\item{...}{additional arguments}

\item{x}{[\code{named list}]\cr
\code{\link{LDABatch}} or \code{\link{LDARep}} object.}

\item{vocab}{[\code{character}]\cr
Vocabularies taken into consideration for merging topic matrices.
Not considered, if \code{sclop} is passed. Default is the vocabulary of the first LDA.}

\item{limit.rel}{[0,1]\cr
See \code{\link{jaccardTopics}}. Default is \code{1/500}.
Not considered for calculation, if \code{sclop} is passed. But should be
passed determining the correct value for the resulting object.}

\item{limit.abs}{[\code{integer(1)}]\cr
See \code{\link{jaccardTopics}}. Default is \code{10}.
Not considered for calculation, if \code{sclop} is passed. But should be
passed determining the correct value for the resulting object.}

\item{atLeast}{[\code{integer(1)}]\cr
See \code{\link{jaccardTopics}}. Default is \code{0}.
Not considered for calculation, if \code{sclop} is passed. But should be
passed determining the correct value for the resulting object.}

\item{progress}{[\code{logical(1)}]\cr
Should a nice progress bar be shown for the steps of \code{\link{mergeTopics}}
and \code{\link{jaccardTopics}}? Turning it off, could lead to significantly
faster calculation. Default ist \code{TRUE}.
Not considered, if \code{sclop} is passed.}

\item{pm.backend}{[\code{character(1)}]\cr
One of "multicore", "socket" or "mpi".
If \code{pm.backend} is set, \code{\link[parallelMap]{parallelStart}} is
called before computation is started and \code{\link[parallelMap]{parallelStop}}
is called after.
Not considered, if \code{sclop} is passed.}

\item{ncpus}{[\code{integer(1)}]\cr
Number of (physical) CPUs to use. If \code{pm.backend} is passed,
default is determined by \code{\link[future]{availableCores}}.
Not considered, if \code{sclop} is passed.}

\item{keepTopics}{[\code{logical(1)}]\cr
Should the merged topic matrix from \code{\link{mergeTopics}} be kept?
Not considered, if \code{sclop} is passed.}

\item{keepSims}{[\code{logical(1)}]\cr
Should the calculated topic similarities matrix from \code{\link{jaccardTopics}}
be kept? Not considered, if \code{sclop} is passed.}

\item{keepLDAs}{[\code{logical(1)}]\cr
Should the considered LDAs be kept?}

\item{sclop}{[\code{symmetrical named matrix}]\cr
(optional) All pairwise S-CLOP scores of the given LDA runs determined by
\code{\link{SCLOP.pairwise}}. Matching of names is not implemented yet, so order matters.}

\item{lda}{[\code{named list}]\cr
List of \code{\link{LDA}} objects, named by the corresponding "job.id".}

\item{id}{[\code{character(1)}]\cr
A name for the computation. If not passed, it is set to "LDARep".
Not considered for \code{\link{LDABatch}} or \code{\link{LDARep}} objects.}

\item{job}{[\code{\link{data.frame}} or \code{named vector}]\cr
A data.frame or data.table with named columns (at least)
"job.id" (\code{integerish}), "K", "alpha", "eta" and "num.iterations"
or a named vector with entries (at least) "K", "alpha", "eta" and "num.iterations".
If not passed, it is interpreted from \code{param} of each LDA.
Not considered for \code{\link{LDABatch}} or \code{\link{LDARep}} objects.}
}
\value{
[\code{named list}] with entries
 \describe{
  \item{\code{id}}{[\code{character(1)}] See above.}
  \item{\code{protoid}}{[\code{character(1)}] Name (ID) of the determined Prototype LDA.}
  \item{\code{lda}}{List of \code{\link{LDA}} objects of the determined Prototype LDA
  and - if \code{keepLDAs} is \code{TRUE} - all considered LDAs.}
  \item{\code{jobs}}{[\code{data.table}] with parameter specifications for the LDAs.}
  \item{\code{param}}{[\code{named list}] with parameter specifications for
  \code{limit.rel} [0,1], \code{limit.abs} [\code{integer(1)}] and
  \code{atLeast} [\code{integer(1)}]. See above for explanation.}
  \item{\code{topics}}{[\code{named matrix}] with the count of vocabularies
  (row wise) in topics (column wise).}
  \item{\code{sims}}{[\code{lower triangular named matrix}] with all pairwise
  jaccard similarities of the given topics.}
  \item{\code{wordslimit}}{[\code{integer}] with counts of words determined as
  relevant based on \code{limit.rel} and \code{limit.abs}.}
  \item{\code{wordsconsidered}}{[\code{integer}] with counts of considered
  words for similarity calculation. Could differ from \code{wordslimit}, if
  \code{atLeast} is greater than zero.}
  \item{\code{sclop}}{[\code{symmetrical named matrix}] with all pairwise
  S-CLOP scores of the given LDA runs.}
}
}
\description{
Returns the Prototype LDA of a set of LDAs. This set is given as
\code{\link{LDABatch}} object, \code{\link{LDARep}} object, or as list of LDAs.
If the matrix of S-CLOP scores \code{sclop} is passed, no calculation is needed/done.
}
\details{
While \code{\link{LDAPrototype}} marks the overall shortcut for performing
multiple LDA runs and choosing the Prototype of them, \code{getPrototype}
just hooks up at determining the Prototype. The generation of multiple LDAs
has to be done before use of this function. The function is flexible enough
to use it at at least two steps/parts of the analysis: After generating the
LDAs (no matter whether as LDABatch or LDARep object) or after determing
the pairwise SCLOP values.

To save memory a lot of interim calculations are discarded by default.

If you use parallel computation, no progress bar is shown.

For details see the details sections of the workflow functions.
}
\examples{
res = LDARep(docs = reuters_docs, vocab = reuters_vocab,
   n = 4, K = 10, num.iterations = 30)
topics = mergeTopics(res, vocab = reuters_vocab)
jacc = jaccardTopics(topics, atLeast = 2)
dend = dendTopics(jacc)
sclop = SCLOP.pairwise(jacc)

getPrototype(lda = getLDA(res), sclop = sclop)

proto = getPrototype(res, vocab = reuters_vocab, keepSims = TRUE,
   limit.abs = 20, atLeast = 10)
proto
getPrototype(proto) # = getLDA(proto)
getConsideredWords(proto)
# > 10 if there is more than one word which is the 10-th often word (ties)
getRelevantWords(proto)
getSCLOP(proto)
}
\seealso{
Other shortcut functions: 
\code{\link{LDAPrototype}()}

Other PrototypeLDA functions: 
\code{\link{LDAPrototype}()},
\code{\link{getSCLOP}()}

Other workflow functions: 
\code{\link{LDARep}()},
\code{\link{SCLOP}()},
\code{\link{dendTopics}()},
\code{\link{jaccardTopics}()},
\code{\link{mergeTopics}()}
}
\concept{PrototypeLDA functions}
\concept{shortcut functions}
\concept{workflow functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getJob.R
\name{getJob}
\alias{getJob}
\alias{getID}
\alias{getRegistry}
\alias{getLDA}
\alias{setFileDir}
\title{Getter and Setter for LDARep and LDABatch}
\usage{
getJob(x)

getID(x)

getRegistry(x)

getLDA(x, job, reduce, all)

setFileDir(x, file.dir)
}
\arguments{
\item{x}{[\code{named list}]\cr
\code{\link{LDABatch}} or \code{\link{LDARep}} object.}

\item{job}{[\code{\link{data.frame}} or \code{integer}]\cr
A data.frame or data.table with a column named "job.id" or a vector of integerish job ids.}

\item{reduce}{[\code{logical(1)}]\cr
If the list of LDAs contains only one element, should the list be reduced and
the single (unnamed) element be returned? Default is \code{TRUE}.}

\item{all}{not implemented for \code{\link{LDABatch}} and \code{\link{LDARep}}
object. See \code{\link[=getSCLOP]{getLDA}}}

\item{file.dir}{[Vector to be coerced to a \code{\link[fs]{fs_path}} object.]\cr
New file directory to overwrite the registry's old one. This can be useful
if the registry is transferred from a batch system.}
}
\description{
Returns the job ids and its parameter set (\code{getJob}) or the (registry's)
id (\code{getID}) for a \code{\link{LDABatch}} or \code{\link{LDARep}} object.
\code{getRegistry} returns the registry itself for a \code{\link{LDABatch}}
object. \code{getLDA} returns the list of \code{\link{LDA}} objects for a
\code{\link{LDABatch}} or \code{\link{LDARep}} object. In addition, you can
specify one or more LDAs by their id(s).\cr
\code{setFilDir} sets the registry's file directory for a
\code{\link{LDABatch}} object. This is useful if you move the registry´s folder,
e.g. if you do your calculations on a batch system, but want to do your
evaluation on your desktop computer.
}
\seealso{
Other getter functions: 
\code{\link{getSCLOP}()},
\code{\link{getSimilarity}()},
\code{\link{getTopics}()}

Other replication functions: 
\code{\link{LDAPrototype}()},
\code{\link{LDARep}()},
\code{\link{as.LDARep}()},
\code{\link{mergeRepTopics}()}

Other batch functions: 
\code{\link{LDABatch}()},
\code{\link{as.LDABatch}()},
\code{\link{mergeBatchTopics}()}
}
\concept{batch functions}
\concept{getter functions}
\concept{replication functions}
