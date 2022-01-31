# omicwas

Cell-Type-Specific Disease Association Testing in Bulk Omics Experiments

## Installation in R

In order to install the most recent version:

    install.packages("devtools")
    devtools::install_github("fumi-github/omicwas", build_vignettes = TRUE)

To install from CRAN archive (possibly a version older than github):

    install.packages("omicwas")

To uninstall package:

    remove.packages("omicwas")

## Usage

    library(omicwas)
    vignette("intro", package = "omicwas")

## Information

Please cite

[Takeuchi, F., Kato, N. Nonlinear ridge regression improves cell-type-specific differential expression analysis. BMC Bioinformatics 22, 141 (2021)](https://doi.org/10.1186/s12859-021-03982-3)
---
title: "Introduction to omicwas"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Introduction to omicwas}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## Installation

```{r, eval = FALSE}
install.packages("devtools")
devtools::install_github("fumi-github/omicwas")
```

If you encounter dependency error for `sva` package, install it from Bioconductor:

```{r, eval = FALSE}
if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
BiocManager::install("sva")
```

## Usage

`omicwas` is a package for cell-type-specific disease association testing,
using bulk tissue data as input.
The package accepts DNA methylation data for epigenome-wide association studies (EWAS),
as well as gene expression data for differential gene expression analyses.
The main function is `ctassoc`

```{r setup}
library(omicwas)
```

See description.

```{r}
?ctassoc
```

## Analyzing DNA methylation

Let's load a sample data.

```{r}
data(GSE42861small)
X = GSE42861small$X
W = GSE42861small$W
Y = GSE42861small$Y
Y = Y[seq(1, 20), ] # for brevity
C = GSE42861small$C
```

See description.

```{r}
?GSE42861small
```

The conventional way is to use ordinary linear regression.

```{r}
result = ctassoc(X, W, Y, C = C,
                 test = "full")
result$coefficients
```

We recommend nonlinear regression with ridge regularization.
For DNA methylation, we use the **logit** function for normalization,
and the test option is `nls.logit`

```{r}
result = ctassoc(X, W, Y, C = C,
                 test = "nls.logit",
                 regularize = TRUE)
print(result$coefficients, n = 20)
```

The first 19 lines show the result for CpG site cg10543797.
Line 1 shows that the basal methylation level in CD4.
(actually CD4+ T cells) is
plogis(1.498) = 0.817, so this cell type is 81% methylated.
Line 8 shows that the CD4.-specific effect of the disease
is 7.10e-04	(in logit scale).
Since the p.value is 0.64, this is not significant.
Line 15 shows that the effect of sexF (female compared to male) is
-4.14e-02 with a small p.value 9.15e-05.
Since sexF is a covariate that has uniform effect across cell types,
the celltype column is NA.

## Analyzing gene expression

Let's load a sample data.

```{r}
data(GTExsmall)
X = GTExsmall$X
W = GTExsmall$W
Y = GTExsmall$Y + 1
Y = Y[seq(1, 20), ] # for brevity
C = GTExsmall$C
```

See description.

```{r}
?GTExsmall
```

The conventional way is to use ordinary linear regression.

```{r}
result = ctassoc(X, W, Y, C = C,
                 test = "full")
result$coefficients
```

We recommend nonlinear regression with ridge regularization.
For DNA methylation, we use the **log** function for normalization,
and the test option is `nls.log`

```{r}
result = ctassoc(X, W, Y, C = C,
                 test = "nls.log",
                 regularize = TRUE)
print(result$coefficients, n = 15)
```

The first 13 lines show the result for transcript ENSG00000059804.
Line 1 shows that the basal expression level in Granulocytes is
exp(8.847) = 6955.
Line 7 shows that the Granulocytes-specific effect of age
is 4.62e-04	(in log scale).
Since the p.value is 0.82, this is not significant.
Line 13 shows that the effect of sexF (female compared to male) is
-2.57e-03 with p.value 0.97
Since sexF is a covariate that has uniform effect across cell types,
the celltype column is NA.

## Analyzing mQTL

For QTL analyses, we use `ctcisQTL` function instead of `ctassoc`.
To speed up computation, we perform linear ridge regression,
thus the statistical test is almost identical to `ctassoc(test = "nls.identity", regularize = TRUE)`.
We analyze only in the linear scale.
Association analysis is performed between each row of Y and each row of X.
See description.

```{r}
?ctcisQTL
```

Let's load a sample data.

```{r}
data(GSE79262small)
X    = GSE79262small$X
Xpos = GSE79262small$Xpos
W    = GSE79262small$W
Y    = GSE79262small$Y
Ypos = GSE79262small$Ypos
C    = GSE79262small$C
X    = X[seq(1, 3001, 100), ] # for brevity
Xpos = Xpos[seq(1, 3001, 100)]
Y    = Y[seq(1, 501, 100), ]
Ypos = Ypos[seq(1, 501, 100)]
```

See description.

```{r}
?GSE79262small
```

Analyze mQTL.

```{r}
ctcisQTL(X, Xpos, W, Y, Ypos, C = C)
```

The result is stored in a file.

```{r}
head(
  read.table(file.path(tempdir(), "ctcisQTL.out.txt"),
             header = TRUE,
             sep ="\t"))
```

The first 3 lines show the result for the association of
SNP rs6678176 with CpG site cg19251656.
Line 1 shows that the CD4T-specific effect of the SNP is -0.003.
Since the p.value is 0.75, this is not significant.
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GSE79262small.R
\docType{data}
\name{GSE79262small}
\alias{GSE79262small}
\title{Small Subset of GSE79262 Dataset From GEO}
\format{
An object of class \code{list} of length 6.
}
\source{
\href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE79262}{GEO}
}
\usage{
data(GSE79262small)
}
\description{
The dataset includes 53 samples.
A subset of 737 CpG sites and 3624 SNPs within Chr1:100,000,000-110,000,000
were selected from the original EWAS dataset.
DNA methylation was measured in T cells.
The estimated proportion of CD4T, CD8T, NK cells are saved in W.
}
\examples{
data(GSE79262small)
X    = GSE79262small$X
Xpos = GSE79262small$Xpos
W    = GSE79262small$W
Y    = GSE79262small$Y
Ypos = GSE79262small$Ypos
C    = GSE79262small$C
X    = X[seq(1, 3001, 100), ] # for brevity
Xpos = Xpos[seq(1, 3001, 100)]
Y    = Y[seq(1, 501, 100), ]
Ypos = Ypos[seq(1, 501, 100)]
ctcisQTL(X, Xpos, W, Y, Ypos, C = C)
}
\seealso{
ctcisQTL
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GSE42861small.R
\docType{data}
\name{GSE42861small}
\alias{GSE42861small}
\title{Small Subset of GSE42861 Dataset From GEO}
\format{
An object of class \code{list} of length 4.
}
\source{
\href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42861}{GEO}
}
\usage{
data(GSE42861small)
}
\description{
The dataset includes 336 rheumatoid arthritis cases and 322 controls.
A subset of 500 CpG sites were randomly selected from the original EWAS dataset.
}
\examples{
data(GSE42861small)
X = GSE42861small$X
W = GSE42861small$W
Y = GSE42861small$Y
Y = Y[seq(1, 20), ] # for brevity
C = GSE42861small$C
result = ctassoc(X, W, Y, C = C)
result$coefficients
}
\seealso{
ctassoc
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rrs.R
\name{rrs.fit}
\alias{rrs.fit}
\title{Fitting reduced-rank ridge regression with given rank and shrinkage penalty}
\usage{
rrs.fit(Y, X, nrank = min(ncol(Y), ncol(X)), lambda = 1, coefSVD = FALSE)
}
\arguments{
\item{Y}{a matrix of response (n by q)}

\item{X}{a matrix of covariate (n by p)}

\item{nrank}{an integer specifying the desired rank}

\item{lambda}{tunging parameter for the ridge penalty}

\item{coefSVD}{logical indicating the need for SVD for the
coeffient matrix int the output}
}
\value{
S3 \code{rrr} object, a list consisting of
  \item{coef}{coefficient of rrs}
  \item{coef.ls}{coefficient of least square}
  \item{fitted}{fitted value of rrs}
  \item{fitted.ls}{fitted value of least square}
  \item{A}{right singular matrix}
  \item{Ad}{sigular value vector}
  \item{nrank}{rank of the fitted rrr}
}
\description{
Fitting reduced-rank ridge regression with given rank and shrinkage penalty
This is a modification of rrs.fit in rrpack version 0.1-6.
In order to handle extremely large q = ncol(Y),
generation of a q by q matrix is avoided.
}
\examples{
Y <- matrix(rnorm(400), 100, 4)
X <- matrix(rnorm(800), 100, 8)
rfit <- rrs.fit(Y, X)
}
\references{
Mukherjee, A. and Zhu, J. (2011) Reduced rank ridge regression and its
kernal extensions.

Mukherjee, A., Chen, K., Wang, N. and Zhu, J. (2015) On the degrees of
freedom of reduced-rank estimators in multivariate
regression. \emph{Biometrika}, 102, 457--477.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GTExsmall.R
\docType{data}
\name{GTExsmall}
\alias{GTExsmall}
\title{Small Subset of GTEx Dataset}
\format{
An object of class \code{list} of length 4.
}
\source{
\href{https://gtexportal.org}{GTEx}
}
\usage{
data(GTExsmall)
}
\description{
The dataset includes gene expression measured in whole blood for 389 samples.
A subset of 500 genes were randomly selected from the original dataset.
}
\examples{
data(GTExsmall)
X = GTExsmall$X
W = GTExsmall$W
Y = GTExsmall$Y + 1
Y = Y[seq(1, 20), ] # for brevity
C = GTExsmall$C
result = ctassoc(X, W, Y, C = C)
result$coefficients
}
\seealso{
ctassoc
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ctcisQTL.R
\name{ctcisQTL}
\alias{ctcisQTL}
\title{Cell-Type-Specific QTL analysis}
\usage{
ctcisQTL(
  X,
  Xpos,
  W,
  Y,
  Ypos,
  C = NULL,
  max.pos.diff = 1e+06,
  outdir = tempdir(),
  outfile = "ctcisQTL.out.txt"
)
}
\arguments{
\item{X}{Matrix (or vector) of SNP genotypes; SNPs x samples.}

\item{Xpos}{Vector of the physical position of X}

\item{W}{Matrix of cell type composition; samples x cell types.}

\item{Y}{Matrix (or vector) of bulk omics measurements; markers x samples.}

\item{Ypos}{Vector of the physical position of Y}

\item{C}{Matrix (or vector) of covariates; samples x covariates.
X, Xpos, W, Y, Ypos, C should be numeric.}

\item{max.pos.diff}{Maximum positional difference to compute cis-QTL.
Association analysis is performed between a row of X and a row of Y,
only when they are within this limit.
Since the limiting is only by position, the function needs to be run
separately for each chromosome.}

\item{outdir}{Output directory.}

\item{outfile}{Output file.}
}
\value{
The estimate, statistic, p.value are written to the specified file.
}
\description{
Cell-Type-Specific QTL analysis
}
\details{
A function for analyses of QTL, such as eQTL, mQTL, pQTL.
The statistical test is almost identical to
\code{ctassoc(test =  "nls.identity", regularize = "TRUE")}.
Association analysis is performed between each row of Y and each row of X.
Usually, the former will be a methylation/expression marker,
and the latter will be a SNP.
To cope with the large number of combinations,
the testing is limited to pairs whose position is within
the difference specified by \code{max.pos.diff}; i.e., limited to cis-QTL.
In detail, this function performs linear ridge regression,
whereas \code{ctassoc(test =  "nls.identity", regularize = "TRUE")}
actually is nonlinear regression
but with \eqn{f} = identity as normalizing transformation.
In order to speed up computation, first, the parameters \eqn{\alpha_{h j}} and
\eqn{\gamma_{j l}} are fit by ordinary linear regression assuming \eqn{\beta_{h j k} = 0}.
Next, \eqn{\beta_{h j k}} are fit and tested by
linear ridge regression (see documentation for \link[omicwas]{ctassoc}).
}
\examples{
\donttest{
data(GSE79262small)
X    = GSE79262small$X
Xpos = GSE79262small$Xpos
W    = GSE79262small$W
Y    = GSE79262small$Y
Ypos = GSE79262small$Ypos
C    = GSE79262small$C
X    = X[seq(1, 3601, 100), ] # for brevity
Xpos = Xpos[seq(1, 3601, 100)]
ctcisQTL(X, Xpos, W, Y, Ypos, C = C)
}

}
\seealso{
ctassoc
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/omicwas.R
\name{ctassoc}
\alias{ctassoc}
\title{Cell-Type-Specific Association Testing}
\usage{
ctassoc(
  X,
  W,
  Y,
  C = NULL,
  test = "full",
  regularize = FALSE,
  num.cores = 1,
  chunk.size = 1000,
  seed = 123
)
}
\arguments{
\item{X}{Matrix (or vector) of traits; samples x traits.}

\item{W}{Matrix of cell type composition; samples x cell types.}

\item{Y}{Matrix (or vector) of bulk omics measurements; markers x samples.}

\item{C}{Matrix (or vector) of covariates; samples x covariates.
X, W, Y, C should be numeric.}

\item{test}{Statistical test to apply; either \code{"full"}, \code{"marginal"},
\code{"nls.identity"}, \code{"nls.log"}, \code{"nls.logit"},
\code{"propdiff.identity"}, \code{"propdiff.log"}, \code{"propdiff.logit"}
or \code{"reducedrankridge"}.}

\item{regularize}{Whether to apply Tikhonov (ie ridge) regularization
to \eqn{\beta_{h j k}}.
The regularization parameter is chosen automatically according to
an unbiased version of (Lawless & Wang, 1976).
Effective for \code{nls.*} and \code{propdiff.*} tests.}

\item{num.cores}{Number of CPU cores to use.
Full, marginal and propdiff tests are run in serial, thus num.cores is ignored.}

\item{chunk.size}{The size of job for a CPU core in one batch.
If you have many cores but limited memory, and there is a memory failure,
decrease num.cores and/or chunk.size.}

\item{seed}{Seed for random number generation.}
}
\value{
A list with one element, which is named "coefficients".
The element gives the estimate, statistic, p.value in tibble format.
In order to transform the estimate for \eqn{\alpha_{h j}} to the original scale,
apply \code{plogis} for \code{test = nls.logit} and
\code{exp} for \code{test = nls.log}.
The estimate for \eqn{\beta_{h j k}} by \code{test = nls.log} is
the natural logarithm of fold-change, not the log2.
If numerical convergence fails, \code{NA} is returned for that marker.
}
\description{
Cell-Type-Specific Association Testing
}
\details{
Let the indexes be
\eqn{h} for cell type, \eqn{i} for sample,
\eqn{j} for marker (CpG site or gene),
\eqn{k} for each trait that has cell-type-specific effect,
and \eqn{l} for each trait that has a uniform effect across cell types.
The input data are \eqn{X_{i k}}, \eqn{C_{i l}}, \eqn{W_{i h}} and \eqn{Y_{j i}},
where \eqn{C_{i l}} can be omitted.
\eqn{X_{i k}} and \eqn{C_{i l}} are the values for two types of traits,
showing effects that are cell-type-specific or not, respectively.
Thus, calling \eqn{X_{i k}} and \eqn{C_{i l}} as "traits" and "covariates"
gives a rough idea, but is not strictly correct.
\eqn{W_{i h}} represents the cell type composition and
\eqn{Y_{j i}} represents the marker level,
such as methylation or gene expression.
For each tissue sample, the cell type proportion \eqn{W_{i h}}
is the proportion of each cell type in the bulk tissue,
which is measured or imputed beforehand.
The marker level \eqn{Y_{j i}} in bulk tissue is measured and provided as input.

The parameters we estimate are
the cell-type-specific trait effect \eqn{\beta_{h j k}},
the tissue-uniform trait effect \eqn{\gamma_{j l}},
and the basal marker level \eqn{\alpha_{h j}} in each cell type.

We first describe the conventional linear regression models.
For marker \eqn{j} in sample \eqn{i},
the maker level specific to cell type \eqn{h} is
\deqn{\alpha_{h j} + \sum_k \beta_{h j k} * X_{i k}.}
This is a representative value rather than a mean, because we do not model
a probability distribution for cell-type-specific expression.
The bulk tissue marker level is the average weighted by \eqn{W_{i h}},
\deqn{\mu_{j i} = \sum_h W_{i h} [ \alpha_{h j} + \sum_k \beta_{h j k} * X_{i k} ] +
                  \sum_l \gamma_{j l} C_{i l}.}
The statistical model is
\deqn{Y_{j i} = \mu_{j i} + \epsilon_{j i},}
\deqn{\epsilon_{j i} ~ N(0, \sigma^2_j).}
The error of the marker level is is noramlly distributed with variance
\eqn{\sigma^2_j}, independently among samples.

The \code{full} model is the linear regression
\deqn{Y_{j i} = (\sum_h \alpha_{h j} * W_{i h}) +
                (\sum_{h k} \beta_{h j k} * W_{i h} * X_{i k}) +
                (\sum_l \gamma_{j l} * C_{i l}) +
                error.}
The \code{marginal} model tests the trait association only in one
cell type \eqn{h}, under the linear regression,
\deqn{Y_{j i} = (\sum_{h'} \alpha_{h' j} * W_{i h'}) +
                (\sum_k \beta_{h j k} * W_{i h} * X_{i k}) +
                (\sum_l \gamma_{j l} * C_{i l}) +
                error.}

The nonlinear model simultaneously analyze cell type composition in
linear scale and differential expression/methylation in log/logit scale.
The normalizing function is the natural logarithm \eqn{f} = log for gene
expression, and \eqn{f} = logit for methylation. Conventional linear regression
can be formulated by defining \eqn{f} as the identity function. The three models
are named \code{nls.log}, \code{nls.logit} and \code{nls.identity}.
We denote the inverse function of \eqn{f} by \eqn{g}; \eqn{g} = exp for
gene expression, and \eqn{g} = logistic for methylation.
The mean normalized marker level of marker \eqn{j} in sample \eqn{i} becomes
\deqn{\mu_{j i} = f(\sum_h W_{i h} g( \alpha_{h j} + \sum_k \beta_{h j k} * X_{i k} )) +
                  \sum_l \gamma_{j l} C_{i l}.}
The statistical model is
\deqn{f(Y_{j i}) = \mu_{j i} + \epsilon_{j i},}
\deqn{\epsilon_{j i} ~ N(0, \sigma^2_j).}
The error of the marker level is is noramlly distributed with variance
\eqn{\sigma^2_j}, independently among samples.

The ridge regression aims to cope with multicollinearity of
the interacting terms \eqn{W_{i h} * X_{i k}}.
Ridge regression is fit by minimizing the residual sum of squares (RSS) plus
\eqn{\lambda \sum_{h k} \beta_{h j k}^2}, where \eqn{\lambda > 0} is the
regularization parameter.

The propdiff tests try to cope with multicollinearity by, roughly speaking,
using mean-centered \eqn{W_{i h}}.
We obtain, instead of \eqn{\beta_{h j k}}, the deviation of
\eqn{\beta_{h j k}} from the average across cell types.
Accordingly, the null hypothesis changes.
The original null hypothesis was \eqn{\beta_{h j k} = 0}.
The null hypothesis when centered is
\eqn{\beta_{h j k} - (\sum_{i h'} W_{i h'} \beta_{h' j k}) / (\sum_{i h'} W_{i h'}) = 0}.
It becomes difficult to detect a signal for a major cell type,
because \eqn{\beta_{h j k}} would be close to the average across cell types.
The tests \code{propdiff.log} and \code{propdiff.logit} include
an additional preprocessing step that converts \eqn{Y_{j i}} to \eqn{f(Y_{j i})}.
Apart from the preprocessing, the computations are performed in linear scale.
As the preprocessing distorts the linearity between the dependent variable
and (the centered) \eqn{W_{i h}},
I actually think \code{propdiff.identity} is better.
}
\examples{
\donttest{
data(GSE42861small)
X = GSE42861small$X
W = GSE42861small$W
Y = GSE42861small$Y
C = GSE42861small$C
result = ctassoc(X, W, Y, C = C)
result$coefficients
}

}
\references{
Lawless, J. F., & Wang, P. (1976). A simulation study of ridge and other
regression estimators.
Communications in Statistics - Theory and Methods, 5(4), 307–323.
\url{https://doi.org/10.1080/03610927608827353}
}
\seealso{
ctcisQTL
}
