Robust backfitting
================
Matias Salibian
2021-04-07

## A robust backfitting algorithm

The `R` package `RBF` (available on CRAN
[here](https://cran.r-project.org/package=RBF)) implements the robust
back-fitting algorithm as proposed by Boente, Martinez and
Salibian-Barrera in

> Boente G, Martinez A, Salibian-Barrera M. (2017) Robust estimators for
> additive models using backfitting. Journal of Nonparametric
> Statistics. Taylor & Francis; 29, 744-767. [DOI:
> 10.1080/10485252.2017.1369077](https://doi.org/10.1080/10485252.2017.1369077)

A paper about this software package appeared in the [The Journal of Open
Source Software](https://joss.theoj.org/) and can be found here:
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02992/status.svg)](https://doi.org/10.21105/joss.02992)

This repository contains a development version of `RBF` which may differ
slightly from the one available on CRAN (until the CRAN version is
updated appropriately).

The package in this repository can be installed from within `R` by using
the following code (assuming the
[devtools](https://cran.r-project.org/package=devtools)) package is
available:

``` r
devtools::install_github("msalibian/RBF")
```

### An example

Here is a (longish) example on how `RBF` works. We use the Air Quality
data. The interest is in predicting `Ozone` in terms of three
explanatory variables: `Solar.R`, `Wind` and `Temp`:

``` r
library(RBF)
data(airquality)
pairs(airquality[, c('Ozone', 'Solar.R', 'Wind', 'Temp')], 
      pch=19, col='gray30', cex=1.5)
```

![](man/figures/intro-1.png)<!-- -->

The following bandwidths were obtained via a robust leave-one-out
cross-validation procedure (described in the paper). As a robust
prediction error measure we use `mu^2 + sigma^2` where `mu` and `sigma`
are M-estimators of location and scale of the prediction errors,
respectively. The code is copied below:

``` r
# Bandwidth selection with leave-one-out cross-validation
## Without outliers
# This takes a long time to compute (approx 380 minutes running
# R 3.6.1 on an Intel(R) Core(TM) i7-4790 CPU @ 3.60GHz)
ccs <- complete.cases(airquality)
x <- as.matrix( airquality[ccs, c('Solar.R', 'Wind', 'Temp')] )
y <- as.vector( airquality[ccs, 'Ozone'] )
a <- c(1/2, 1, 1.5, 2, 2.5, 3)
h1 <- a * sd(x[,1])
h2 <- a * sd(x[,2])
h3 <- a * sd(x[,3])
hh <- expand.grid(h1, h2, h3)
nh <- nrow(hh)
rmspe <- rep(NA, nh)
jbest <- 0
cvbest <- +Inf
# leave-one-out
n <- nrow(x)
for(i in 1:nh) {
  # leave-one-out CV loop
  preds <- rep(NA, n)
  for(j in 1:n) {
    tmp <- try( backf.rob(y ~ x, point = x[j, ],
                          windows = hh[i, ], epsilon = 1e-6,
                          degree = 1, type = 'Tukey', subset = c(-j) ))
    if (class(tmp)[1] != "try-error") {
      preds[j] <- rowSums(tmp$prediction) + tmp$alpha
    }
  }
  pred.res <- preds - y
  tmp.re <- RobStatTM::locScaleM(pred.res, na.rm=TRUE)
  rmspe[i] <- tmp.re$mu^2 + tmp.re$disper^2
  if( rmspe[i] < cvbest ) {
    jbest <- i
    cvbest <- rmspe[i]
    print('Record')
  }
  print(c(i, rmspe[i]))
}
bandw <- hh[jbest,]
```

Here we just set them to their optimal values:

``` r
bandw <- c(136.7285, 10.67314, 4.764985)
```

Now we use the robust backfitting algorithm to fit an additive model
using Tukey’s bisquare loss (the default tuning constant for this loss
function is 4.685). We remove cases with missing entries.

``` r
ccs <- complete.cases(airquality)
fit.full <- backf.rob(Ozone ~ Solar.R + Wind + Temp, data=airquality,
                subset=ccs, windows=bandw, degree=1, type='Tukey')
```

We display the 3 fits (one per additive component), being careful with
the axis limits (which are stored to use them later):

``` r
lim.cl <- lim.rob <- matrix(0, 2, 3)
x0 <- fit.full$Xp
for(j in 1:3) {
  re <- fit.full$y - fit.full$alpha - rowSums(fit.full$g.matrix[,-j])
  lim.rob[,j] <- c(min(re), max(re))
  plot(re ~ x0[,j], type='p', pch=19, col='gray30', 
       xlab=colnames(x0)[j], ylab='', cex=1.5)
  oo <- order(x0[,j])
  lines(x0[oo,j], fit.full$g.matrix[oo,j], lwd=5, col='blue')
}
```

<img src="man/figures/showfits-1.png" width="33%" /><img src="man/figures/showfits-2.png" width="33%" /><img src="man/figures/showfits-3.png" width="33%" />

NOTE: These plots could also be obtained using the `plot` method:

``` r
# Plot each component
plot(fit.full, which=1:3)
```

We now compute and display the classical backfitting fits, with
bandwidths chosen via leave-one-out CV:

``` r
library(gam)
x <- aircomplete[, c('Ozone', 'Solar.R', 'Wind', 'Temp')]
a <- c(.3, .4, .5, .6, .7, .8, .9)
hh <- expand.grid(a, a, a)
nh <- nrow(hh)
jbest <- 0
cvbest <- +Inf
n <- nrow(x)
for(i in 1:nh) {
  fi <- rep(0, n)
  for(j in 1:n) {
    tmp <- gam(Ozone ~ lo(Solar.R, span=hh[i,1]) + lo(Wind, span=hh[i,2])
               + lo(Temp, span=hh[i,3]), data=x, subset=c(-j))
    fi[j] <- as.numeric(predict(tmp, newdata=x[j,], type='response'))
  }
  ss <- mean((x$Ozone - fi)^2)
  if(ss < cvbest) {
    jbest <- i
    cvbest <- ss
  }
  print(c(i, ss))
}
(hh[jbest,])
```

The optimal bandwidths are `0.7`, `0.7` and `0.5` for `Solar.R`, `Wind`
and `Temp`, respectively. Below are plots of partial residuals with the
classical and robust fits overlaid:

``` r
aircomplete <- airquality[ complete.cases(airquality), ]
library(gam)
fit.gam <- gam(Ozone ~ lo(Solar.R, span=.7) + lo(Wind, span=.7)+
                 lo(Temp, span=.5), data=aircomplete)
# Plot both fits (robust and classical) 
x <- as.matrix( aircomplete[ , c('Solar.R', 'Wind', 'Temp')] )
y <- as.vector( aircomplete[ , 'Ozone'] )
fits <- predict(fit.gam, type='terms')
# alpha.gam <- attr(fits, 'constant')
for(j in 1:3) {
  re <- fit.full$yp - fit.full$alpha - rowSums(fit.full$g.matrix[,-j])
  plot(re ~ x[,j], type='p', pch=20, col='gray45', xlab=colnames(x)[j], ylab='')
  oo <- order(x[,j])
  lines(x[oo,j], fit.full$g.matrix[oo,j], lwd=5, col='blue')
  lines(x[oo,j], fits[oo,j], lwd=5, col='magenta')
}
```

<img src="man/figures/classicfits-1.png" width="33%" /><img src="man/figures/classicfits-2.png" width="33%" /><img src="man/figures/classicfits-3.png" width="33%" />

To identify potential outiers we look at the residuals from the robust
fit, and use the function `boxplot`:

``` r
re.ro <- residuals(fit.full)
ou.ro <- boxplot(re.ro, col='gray80')$out
in.ro <- (1:length(re.ro))[ re.ro %in% ou.ro ]
points(rep(1, length(in.ro)), re.ro[in.ro], pch=20, cex=1.5, col='red')
```

![](man/figures/outliers-1.png)<!-- -->

We highlight these suspicious observations on the scatter plot

``` r
cs <- rep('gray30', nrow(aircomplete))
cs[in.ro] <- 'red'
os <- 1:nrow(aircomplete)
os2 <- c(os[-in.ro], os[in.ro])
pairs(aircomplete[os2, c('Ozone', 'Solar.R', 'Wind', 'Temp')], 
      pch=19, col=cs[os2], cex=1.5)
```

![](man/figures/showouts-1.png)<!-- -->

and on the partial residuals plots

``` r
for(j in 1:3) {
  re <- fit.full$yp - fit.full$alpha - rowSums(fit.full$g.matrix[,-j])
  plot(re ~ x[,j], type='p', pch=20, col='gray45', xlab=colnames(x)[j], ylab='')
  points(re[in.ro] ~ x0[in.ro,j], pch=19, col='red', cex=1.5)
  oo <- order(x[,j])
  lines(x[oo,j], fit.full$g.matrix[oo,j], lwd=5, col='blue')
  lines(x[oo,j], fits[oo,j], lwd=5, col='magenta')
}
```

<img src="man/figures/showouts2-1.png" width="33%" /><img src="man/figures/showouts2-2.png" width="33%" /><img src="man/figures/showouts2-3.png" width="33%" />

We now compute the classical backfitting algorithm on the data without
the potential outliers identified by the robust fit (the optimal
smoothing parameters for the non-robust fit were re-computed using
leave-one-out cross-validation on the “clean” data set). Note that now
both fits (robust and non-robust) are  
almost identical.

``` r
# Run the classical backfitting algorithm without outliers
airclean <- aircomplete[-in.ro, ]
fit.gam2 <- gam(Ozone ~ lo(Solar.R, span=.7) + lo(Wind, span=.8)+
                  lo(Temp, span=.3), data=airclean)
fits2 <- predict(fit.gam2, type='terms')
# alpha.gam2 <- attr(fits2, 'constant')
dd2 <- aircomplete[-in.ro, c('Solar.R', 'Wind', 'Temp')]
for(j in 1:3) {
  re <- fit.full$yp - fit.full$alpha - rowSums(fit.full$g.matrix[,-j])
  plot(re ~ x[,j], type='p', pch=20, col='gray45', xlab=colnames(x)[j], ylab='')
  points(re[in.ro] ~ x[in.ro,j], pch=20, col='red', cex=1.5)
  oo <- order(dd2[,j])
  lines(dd2[oo,j], fits2[oo,j], lwd=5, col='magenta')
  oo <- order(x[,j])
  lines(x[oo,j], fit.full$g.matrix[oo,j], lwd=5, col='blue')
}
```

<img src="man/figures/bothonclean-1.png" width="33%" /><img src="man/figures/bothonclean-2.png" width="33%" /><img src="man/figures/bothonclean-3.png" width="33%" />

Finally, we compare the prediction accuracy obtained with each of the
fits. Because we are not interested in predicting well any possible
outliers in the data, we evaluate the quality of the predictions using a
5%-trimmed mean squared prediction error (effectively measuring the
prediction accuracy on 95% of the data). We use this alpha-trimmed mean
squared function:

``` r
tms <- function(a, alpha=.1) {
  # alpha is the proportion to trim
  a2 <- sort(a^2, na.last=NA)
  n0 <- floor( length(a) * (1 - alpha) )
  return( mean(a2[1:n0], na.rm=TRUE) )
}
```

We use 100 runs of 5-fold CV to compare the 5%-trimmed mean squared
prediction error of the robust fit and the classical one. Note that the
bandwidths are kept fixed at their optimal value estimated above.

``` r
dd <- airquality
dd <- dd[complete.cases(dd), c('Ozone', 'Solar.R', 'Wind', 'Temp')]
# 100 runs of K-fold CV
M <- 100
# 5-fold
K <- 5
n <- nrow(dd)
# store (trimmed) TMSPE for robust and gam, and also
tmspe.ro <- tmspe.gam <- vector('numeric', M)
set.seed(123)
ii <- (1:n)%%K + 1
for(runs in 1:M) {
  tmpro <- tmpgam <- vector('numeric', n)
  ii <- sample(ii)
  for(j in 1:K) {
    fit.full <- backf.rob(Ozone ~ Solar.R + Wind + Temp, 
                           point=dd[ii==j, -1], windows=bandw, 
                           epsilon=1e-6, degree=1, type='Tukey', 
                           subset = (ii!=j), data = dd)
    tmpro[ ii == j ] <- rowSums(fit.full$prediction) + fit.full$alpha
    fit.gam <- gam(Ozone ~ lo(Solar.R, span=.7) + lo(Wind, span=.7)+
                     lo(Temp, span=.5), data=dd[ii!=j, ])
    tmpgam[ ii == j ] <- predict(fit.gam, newdata=dd[ii==j, ], type='response')
  }
  tmspe.ro[runs] <- tms( dd$Ozone - tmpro, alpha=0.05)
  tmspe.gam[runs] <- tms( dd$Ozone - tmpgam, alpha=0.05)
}
```

These are the boxplots. We see that the robust fit consistently fits the
vast majority (95%) of the data better than the classical one.

``` r
boxplot(tmspe.ro, tmspe.gam, names=c('Robust', 'Classical'), 
        col=c('tomato3', 'gray80')) #, main='', ylim=c(130, 210))
```

![](man/figures/pred3-1.png)<!-- -->

As a sanity check, we compare the prediction accuracy of the robust and
non-robust fits using only the “clean” data set. We re-compute the
optimal bandwidths for the robust fit using leave-one-out cross
validation as above, and as above, note that these bandwidths are kept
fixed. We use 100 runs of 5-fold cross-validation, and compute both the
trimmed and the regular mean squared prediction errors of each fit.

``` r
aq <- airquality
aq2 <- aq[complete.cases(aq), c('Ozone', 'Solar.R', 'Wind', 'Temp')]
airclean <- aq2[ -in.ro, ]
bandw <- c(138.2699, 10.46753, 4.828436)
M <- 100 
K <- 5
n <- nrow(airclean)
mspe.ro <- mspe.gam <- tmspe.ro <- tmspe.gam <- vector('numeric', M)
set.seed(17)
ii <- (1:n)%%K + 1
for(runs in 1:M) {
  tmpro <- tmpgam <- vector('numeric', n)
  ii <- sample(ii)
  for(j in 1:K) {
    fit.full <- try( backf.rob(Ozone ~ Solar.R + Wind + Temp, 
                          point=airclean[ii==j, -1], windows=bandw, 
                          epsilon=1e-6, degree=1, type='Tukey', 
                          subset = (ii!=j), data = airclean) )
    if (class(fit.full)[1] != "try-error") {
      tmpro[ ii == j ] <- rowSums(fit.full$prediction) + fit.full$alpha
    }
    fit.gam <- gam(Ozone ~ lo(Solar.R, span=.7) + lo(Wind, span=.8)+
                     lo(Temp, span=.3), data=airclean, subset = (ii!=j) )
    tmpgam[ ii == j ] <- predict(fit.gam, newdata=airclean[ii==j, ], 
                                 type='response')
  }
  tmspe.ro[runs] <- tms( airclean$Ozone - tmpro, alpha=0.05)
  mspe.ro[runs] <- mean( ( airclean$Ozone - tmpro)^2, na.rm=TRUE)
  tmspe.gam[runs] <- tms( airclean$Ozone - tmpgam, alpha=0.05)
  mspe.gam[runs] <- mean( ( airclean$Ozone - tmpgam)^2, na.rm=TRUE)
}
```

The boxplots of the trimmed and regular mean squared prediction errors
over the 100 cross-validation runs are below. We see that for the
majority of the runs both estimators provide very similar prediction
errors. Note that we naturally expect a robust method to perform
slightly worse than the classical one when no model deviations occur.
The boxplots below show that for this robust backfitting estimator, this
loss in prediction accuracy is in fact very small.

``` r
boxplot(tmspe.ro, tmspe.gam, names=c('Robust', 'Classical'), 
        col=rep(c('tomato3', 'gray80'), 2), main='Trimmed MSPE On "clean" data')
boxplot(mspe.ro, mspe.gam, names=c('Robust', 'Classical'), 
        col=rep(c('tomato3', 'gray80'), 2), main='Non-Trimmed MSPE On "clean" data')
```

<img src="man/figures/boxplot.clean.predictions-1.png" width="33%" /><img src="man/figures/boxplot.clean.predictions-2.png" width="33%" />
## RBF 2.1.0

 - Add a vignette.
 - Fix outputs in the S3 methods `summary` and `print`.
 - Add `formula` to output of print methods.
 - Fix problem for predicting with the argument `point` of function `backf.cl`.
 - Updated README file with link to JOSS paper.
 ---
title: "Robust backfitting"
author: "Matias Salibian"
date: "`r format(Sys.Date())`"
output: github_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, warning=FALSE, message=FALSE,
                      fig.width=6, fig.height=6, fig.path = "man/figures/")
```

## A robust backfitting algorithm

The `R` package `RBF` (available on CRAN [here](https://cran.r-project.org/package=RBF)) 
implements the robust back-fitting algorithm as proposed by
Boente, Martinez and Salibian-Barrera  in 

> Boente G, Martinez A, Salibian-Barrera M. (2017) Robust estimators for 
> additive models using backfitting. Journal of Nonparametric Statistics. Taylor 
> & Francis; 29, 744-767.
> [DOI: 10.1080/10485252.2017.1369077](https://doi.org/10.1080/10485252.2017.1369077)

A paper about this software package appeared in the 
[The Journal of Open Source Software](https://joss.theoj.org/) and
can be found here:
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02992/status.svg)](https://doi.org/10.21105/joss.02992)

This repository contains a development version of `RBF`
which may differ slightly from the one available on CRAN
(until the CRAN version is updated appropriately). 

The package in this repository can be installed from within `R` by using the following code (assuming the [devtools](https://cran.r-project.org/package=devtools)) package is available:
```R
devtools::install_github("msalibian/RBF")
```

### An example

Here is a (longish) example on how `RBF` works. 
We use the Air Quality data. The interest is in
predicting `Ozone` in terms of three explanatory 
variables: `Solar.R`, `Wind` and `Temp`:
```{r intro}
library(RBF)
data(airquality)
pairs(airquality[, c('Ozone', 'Solar.R', 'Wind', 'Temp')], 
      pch=19, col='gray30', cex=1.5)
```

The following bandwidths were obtained via a robust 
leave-one-out cross-validation procedure (described in the paper).
As a robust prediction error measure we use `mu^2 + sigma^2` where
`mu` and `sigma` are M-estimators of location and scale of
the prediction errors, respectively. The code is copied below:
```{r robust.leaveoneout, eval=FALSE}
# Bandwidth selection with leave-one-out cross-validation
## Without outliers
# This takes a long time to compute (approx 380 minutes running
# R 3.6.1 on an Intel(R) Core(TM) i7-4790 CPU @ 3.60GHz)
ccs <- complete.cases(airquality)
x <- as.matrix( airquality[ccs, c('Solar.R', 'Wind', 'Temp')] )
y <- as.vector( airquality[ccs, 'Ozone'] )
a <- c(1/2, 1, 1.5, 2, 2.5, 3)
h1 <- a * sd(x[,1])
h2 <- a * sd(x[,2])
h3 <- a * sd(x[,3])
hh <- expand.grid(h1, h2, h3)
nh <- nrow(hh)
rmspe <- rep(NA, nh)
jbest <- 0
cvbest <- +Inf
# leave-one-out
n <- nrow(x)
for(i in 1:nh) {
  # leave-one-out CV loop
  preds <- rep(NA, n)
  for(j in 1:n) {
    tmp <- try( backf.rob(y ~ x, point = x[j, ],
                          windows = hh[i, ], epsilon = 1e-6,
                          degree = 1, type = 'Tukey', subset = c(-j) ))
    if (class(tmp)[1] != "try-error") {
      preds[j] <- rowSums(tmp$prediction) + tmp$alpha
    }
  }
  pred.res <- preds - y
  tmp.re <- RobStatTM::locScaleM(pred.res, na.rm=TRUE)
  rmspe[i] <- tmp.re$mu^2 + tmp.re$disper^2
  if( rmspe[i] < cvbest ) {
    jbest <- i
    cvbest <- rmspe[i]
    print('Record')
  }
  print(c(i, rmspe[i]))
}
bandw <- hh[jbest,]
```
Here we just set them to their optimal values:
```{r bandw}
bandw <- c(136.7285, 10.67314, 4.764985)
```
Now we use the robust backfitting algorithm to fit an additive
model using Tukey's bisquare loss (the default tuning
constant for this loss function is 4.685). We remove 
cases with missing entries. 
```{r rbfone}
ccs <- complete.cases(airquality)
fit.full <- backf.rob(Ozone ~ Solar.R + Wind + Temp, data=airquality,
                subset=ccs, windows=bandw, degree=1, type='Tukey')
```

We display the 3 fits (one per additive component), being
careful with the axis limits (which are stored to use them later):
```{r showfits, fig.show="hold", out.width="33%"}
lim.cl <- lim.rob <- matrix(0, 2, 3)
x0 <- fit.full$Xp
for(j in 1:3) {
  re <- fit.full$y - fit.full$alpha - rowSums(fit.full$g.matrix[,-j])
  lim.rob[,j] <- c(min(re), max(re))
  plot(re ~ x0[,j], type='p', pch=19, col='gray30', 
       xlab=colnames(x0)[j], ylab='', cex=1.5)
  oo <- order(x0[,j])
  lines(x0[oo,j], fit.full$g.matrix[oo,j], lwd=5, col='blue')
}
```

NOTE: These plots could also be obtained using the `plot` method: 
```{r plot.method, eval=FALSE}
# Plot each component
plot(fit.full, which=1:3)
```

We now compute and display the 
classical backfitting fits, with 
bandwidths chosen via leave-one-out CV:
```{r gam.loo, eval=FALSE}
library(gam)
x <- aircomplete[, c('Ozone', 'Solar.R', 'Wind', 'Temp')]
a <- c(.3, .4, .5, .6, .7, .8, .9)
hh <- expand.grid(a, a, a)
nh <- nrow(hh)
jbest <- 0
cvbest <- +Inf
n <- nrow(x)
for(i in 1:nh) {
  fi <- rep(0, n)
  for(j in 1:n) {
    tmp <- gam(Ozone ~ lo(Solar.R, span=hh[i,1]) + lo(Wind, span=hh[i,2])
               + lo(Temp, span=hh[i,3]), data=x, subset=c(-j))
    fi[j] <- as.numeric(predict(tmp, newdata=x[j,], type='response'))
  }
  ss <- mean((x$Ozone - fi)^2)
  if(ss < cvbest) {
    jbest <- i
    cvbest <- ss
  }
  print(c(i, ss))
}
(hh[jbest,])
```
The optimal bandwidths are `0.7`, `0.7` and `0.5` for `Solar.R`, 
`Wind` and `Temp`, respectively.
Below are plots of  partial residuals with
the classical and robust fits overlaid:
```{r classicfits, fig.show="hold", out.width="33%"}
aircomplete <- airquality[ complete.cases(airquality), ]
library(gam)
fit.gam <- gam(Ozone ~ lo(Solar.R, span=.7) + lo(Wind, span=.7)+
                 lo(Temp, span=.5), data=aircomplete)
# Plot both fits (robust and classical) 
x <- as.matrix( aircomplete[ , c('Solar.R', 'Wind', 'Temp')] )
y <- as.vector( aircomplete[ , 'Ozone'] )
fits <- predict(fit.gam, type='terms')
# alpha.gam <- attr(fits, 'constant')
for(j in 1:3) {
  re <- fit.full$yp - fit.full$alpha - rowSums(fit.full$g.matrix[,-j])
  plot(re ~ x[,j], type='p', pch=20, col='gray45', xlab=colnames(x)[j], ylab='')
  oo <- order(x[,j])
  lines(x[oo,j], fit.full$g.matrix[oo,j], lwd=5, col='blue')
  lines(x[oo,j], fits[oo,j], lwd=5, col='magenta')
}
```

To identify potential outiers 
we look at the residuals from the robust fit, and use
the function `boxplot`:
```{r outliers, fig.width=4, fig.height=4}
re.ro <- residuals(fit.full)
ou.ro <- boxplot(re.ro, col='gray80')$out
in.ro <- (1:length(re.ro))[ re.ro %in% ou.ro ]
points(rep(1, length(in.ro)), re.ro[in.ro], pch=20, cex=1.5, col='red')
```

We highlight these suspicious observations on the
scatter plot 
```{r showouts}
cs <- rep('gray30', nrow(aircomplete))
cs[in.ro] <- 'red'
os <- 1:nrow(aircomplete)
os2 <- c(os[-in.ro], os[in.ro])
pairs(aircomplete[os2, c('Ozone', 'Solar.R', 'Wind', 'Temp')], 
      pch=19, col=cs[os2], cex=1.5)
```

and on the partial residuals plots
```{r showouts2, fig.show="hold", out.width="33%"}
for(j in 1:3) {
  re <- fit.full$yp - fit.full$alpha - rowSums(fit.full$g.matrix[,-j])
  plot(re ~ x[,j], type='p', pch=20, col='gray45', xlab=colnames(x)[j], ylab='')
  points(re[in.ro] ~ x0[in.ro,j], pch=19, col='red', cex=1.5)
  oo <- order(x[,j])
  lines(x[oo,j], fit.full$g.matrix[oo,j], lwd=5, col='blue')
  lines(x[oo,j], fits[oo,j], lwd=5, col='magenta')
}
```

We now compute the classical backfitting algorithm on
the data without the potential outliers identified by
the robust fit (the optimal smoothing 
parameters for the non-robust fit were 
re-computed using leave-one-out cross-validation on the "clean" data set). 
Note that now both fits (robust and non-robust) are  
almost identical.
```{r bothonclean, fig.show="hold", out.width="33%"}
# Run the classical backfitting algorithm without outliers
airclean <- aircomplete[-in.ro, ]
fit.gam2 <- gam(Ozone ~ lo(Solar.R, span=.7) + lo(Wind, span=.8)+
                  lo(Temp, span=.3), data=airclean)
fits2 <- predict(fit.gam2, type='terms')
# alpha.gam2 <- attr(fits2, 'constant')
dd2 <- aircomplete[-in.ro, c('Solar.R', 'Wind', 'Temp')]
for(j in 1:3) {
  re <- fit.full$yp - fit.full$alpha - rowSums(fit.full$g.matrix[,-j])
  plot(re ~ x[,j], type='p', pch=20, col='gray45', xlab=colnames(x)[j], ylab='')
  points(re[in.ro] ~ x[in.ro,j], pch=20, col='red', cex=1.5)
  oo <- order(dd2[,j])
  lines(dd2[oo,j], fits2[oo,j], lwd=5, col='magenta')
  oo <- order(x[,j])
  lines(x[oo,j], fit.full$g.matrix[oo,j], lwd=5, col='blue')
}
```

Finally, we compare the prediction accuracy obtained with each
of the fits. Because we are not interested in predicting well any
possible outliers in the data, we evaluate the quality of
the predictions using a 5%-trimmed mean 
squared prediction error (effectively measuring the prediction 
accuracy on 95% of the data). We use this alpha-trimmed mean 
squared function:
```{r pred1}
tms <- function(a, alpha=.1) {
  # alpha is the proportion to trim
  a2 <- sort(a^2, na.last=NA)
  n0 <- floor( length(a) * (1 - alpha) )
  return( mean(a2[1:n0], na.rm=TRUE) )
}
```

We use 100 runs of 5-fold CV to compare the 
5%-trimmed mean squared prediction error of
the robust fit and the classical one. 
Note that the bandwidths are kept fixed 
at their optimal value estimated above. 
```{r pred2, cache=TRUE}
dd <- airquality
dd <- dd[complete.cases(dd), c('Ozone', 'Solar.R', 'Wind', 'Temp')]
# 100 runs of K-fold CV
M <- 100
# 5-fold
K <- 5
n <- nrow(dd)
# store (trimmed) TMSPE for robust and gam, and also
tmspe.ro <- tmspe.gam <- vector('numeric', M)
set.seed(123)
ii <- (1:n)%%K + 1
for(runs in 1:M) {
  tmpro <- tmpgam <- vector('numeric', n)
  ii <- sample(ii)
  for(j in 1:K) {
    fit.full <- backf.rob(Ozone ~ Solar.R + Wind + Temp, 
                           point=dd[ii==j, -1], windows=bandw, 
                           epsilon=1e-6, degree=1, type='Tukey', 
                           subset = (ii!=j), data = dd)
    tmpro[ ii == j ] <- rowSums(fit.full$prediction) + fit.full$alpha
    fit.gam <- gam(Ozone ~ lo(Solar.R, span=.7) + lo(Wind, span=.7)+
                     lo(Temp, span=.5), data=dd[ii!=j, ])
    tmpgam[ ii == j ] <- predict(fit.gam, newdata=dd[ii==j, ], type='response')
  }
  tmspe.ro[runs] <- tms( dd$Ozone - tmpro, alpha=0.05)
  tmspe.gam[runs] <- tms( dd$Ozone - tmpgam, alpha=0.05)
}
```
These are the boxplots. We see that the robust fit consistently
fits the vast majority (95%) of the data better than the 
classical one. 
```{r pred3}
boxplot(tmspe.ro, tmspe.gam, names=c('Robust', 'Classical'), 
        col=c('tomato3', 'gray80')) #, main='', ylim=c(130, 210))
```

As a sanity check, we compare the prediction accuracy of the robust and 
non-robust fits using only the "clean" data set. We re-compute the 
optimal bandwidths for the robust fit using leave-one-out cross validation
as above, and as above, note that these bandwidths are kept fixed. 
We use 100 runs of 5-fold cross-validation, and
compute both the trimmed and the regular mean squared prediction 
errors of each fit.
```{r pred.clean, cache=TRUE}
aq <- airquality
aq2 <- aq[complete.cases(aq), c('Ozone', 'Solar.R', 'Wind', 'Temp')]
airclean <- aq2[ -in.ro, ]
bandw <- c(138.2699, 10.46753, 4.828436)
M <- 100 
K <- 5
n <- nrow(airclean)
mspe.ro <- mspe.gam <- tmspe.ro <- tmspe.gam <- vector('numeric', M)
set.seed(17)
ii <- (1:n)%%K + 1
for(runs in 1:M) {
  tmpro <- tmpgam <- vector('numeric', n)
  ii <- sample(ii)
  for(j in 1:K) {
    fit.full <- try( backf.rob(Ozone ~ Solar.R + Wind + Temp, 
                          point=airclean[ii==j, -1], windows=bandw, 
                          epsilon=1e-6, degree=1, type='Tukey', 
                          subset = (ii!=j), data = airclean) )
    if (class(fit.full)[1] != "try-error") {
      tmpro[ ii == j ] <- rowSums(fit.full$prediction) + fit.full$alpha
    }
    fit.gam <- gam(Ozone ~ lo(Solar.R, span=.7) + lo(Wind, span=.8)+
                     lo(Temp, span=.3), data=airclean, subset = (ii!=j) )
    tmpgam[ ii == j ] <- predict(fit.gam, newdata=airclean[ii==j, ], 
                                 type='response')
  }
  tmspe.ro[runs] <- tms( airclean$Ozone - tmpro, alpha=0.05)
  mspe.ro[runs] <- mean( ( airclean$Ozone - tmpro)^2, na.rm=TRUE)
  tmspe.gam[runs] <- tms( airclean$Ozone - tmpgam, alpha=0.05)
  mspe.gam[runs] <- mean( ( airclean$Ozone - tmpgam)^2, na.rm=TRUE)
}
```
The boxplots of the trimmed and regular mean squared prediction 
errors over the 100 cross-validation runs are below. We see that
for the majority of the runs both estimators provide very similar
prediction errors. Note that we naturally expect a robust method to 
perform slightly worse than the classical one when no model
deviations occur. The boxplots below show that for this 
robust backfitting estimator, this loss in prediction accuracy is
in fact very small. 
```{r boxplot.clean.predictions, fig.show="hold", out.width="33%"} 
boxplot(tmspe.ro, tmspe.gam, names=c('Robust', 'Classical'), 
        col=rep(c('tomato3', 'gray80'), 2), main='Trimmed MSPE On "clean" data')
boxplot(mspe.ro, mspe.gam, names=c('Robust', 'Classical'), 
        col=rep(c('tomato3', 'gray80'), 2), main='Non-Trimmed MSPE On "clean" data')
```
---
title: "Examples"
author: "Martínez and Salibian"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Examples}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```
# About this vignette

In this vignette we discuss some properties of a
robust backfitting estimator for additive models, and 
illustrate the use of the package `RBF` that implements it. 
These estimators were originally proposed in 
Boente G, Martinez A, Salibian-Barrera M. (2017). 
See also Martinez A. and Salibian-Barrera M. (2021).

Below we analyze two data sets. The first one
shows the robustness of the robust backfitting estimators 
when a small 
proportion of very large outliers are present in the training
set, and compares them with those obtained with the standard backfitting 
algorithm. With the second example, we illustrate how these robust estimators
can be interpreted as automatically detecting and downweighting 
potential outliers in the training set. We also compare the 
prediction accuracy of the robust and classical backfitting 
algorithms. 

<!-- Robust estimators for additive models using backfitting. Journal of Nonparametric Statistics. Taylor & Francis; 29, 744-767. [DOI: 10.1080/10485252.2017.1369077](https://doi.org/10.1080/10485252.2017.1369077) -->


<!-- The `R` package `RBF` (available on CRAN [here](https://cran.r-project.org/package=RBF))  -->
<!-- implements the robust back-fitting algorithm as proposed by -->
<!-- Boente, Martinez and Salibian-Barrera  in  -->

<!-- This repository contains a development version of `RBF` -->
<!-- which may differ slightly from the one available on CRAN -->
<!-- (until the CRAN version is updated appropriately).  -->

<!-- The package in this repository can be installed from within `R` by using the following code (assuming the [devtools](https://cran.r-project.org/package=devtools)) package is available: -->
<!-- ```R -->
<!-- devtools::install_github("msalibian/RBF") -->
<!-- ``` -->
<!-- and charged in the `R` session by the following command -->
<!-- ```{r library} -->
<!-- library("RBF") -->
<!-- ``` -->
<!-- Now that the `R` package is downloaded, we can now start to see how this procedure works. -->

## Boston example

Consider the well-known `Boston` house price data of Harrinson and Rubinfeld (1978). This dataset was used as an example to model an additive model by Härdle et a. (2004). The data are available in the `MASS` package. It contains $n=506$ observations and 14 variables measured on the census districts of the Boston metropolitan area. Following the analysis in Härdle et a. (2004) we use the following 10 explanatory variables:

- `crim`: per capita crime rate by town ($X_1$),
- `indus`: proportion of non-retail business acres per town ($X_2$),
- `nox`: nitric oxides concentration (parts per 10 million) ($X_3$),
- `rm`: average number of rooms per dwelling ($X_4$),
- `age`: proportion of owner-occupied units built prior to 1940 ($X_5$),
- `dis`: weighted distances to five Boston employment centers ($X_6$),
- `tax`: full-value property tax rate per 10,000 ($X_7$),
- `ptratio`: pupil-teacher ratio by town ($X_8$),
- `black`: $1,000(Bk-0.63)^2$ where $Bk$ is the proportion of people of Afrom American descent by town ($X_9$),
- `lstat`: percent lower status of population ($X_{10}$).

The response variable $Y$ is `medv`, the median value of the owner-occupied homes in 1,000 USD, and the proposed  additive model is
$$Y= \mu+ \sum_{j=1}^{10} g_j(\log(X_j))+ \epsilon.$$
where $\mu \in \mathbb{R}$, $g_j$ are unknown smooth functions,
and $\epsilon$ are random errors. 

First we load the data, and transform the explanatory variables as
required by the model above:
```{r read the dataset}
data(Boston, package='MASS')
dd <- Boston[, c(1, 3, 5:8, 10:14)]
dd[, names(dd) != 'medv'] <- log( dd[, names(dd) != 'medv'] )
```
Next, we load the `RBF` package
```{r loadpckg}
library(RBF)
```
The robust backfitting estimators for each additive component are 
computed using robust 
kernel-based local polynomial regression.
The model to be fit is specified using the
standard `formula` notation in `R`. We also need to specify
the following arguments:

- `windows`: the bandwidths for the kernel estimators,
- `degree`: the degree of the  polynomial used for the kernel local regression, defaults to `0`,
- `type`: specifies the robust loss function, options are `Huber` or `Tukey`, defaults to `Huber`.

As with all kernel-based estimators, bandwidth selection is an important step. In this example
we follow 
Härdle et al. (2004) and select bandwidths $h_j$, $1 \le j \le 10$, 
proportional to the standard deviation of the corresponding
explanatory variables $\log(X_j)$. Specifically we set 
$h_j = \hat{\sigma}_j / 2$: 
```{r bandwidths}
bandw <- apply(dd[, names(dd) != 'medv'], 2, sd) / 2
```
We are now ready to compute the robust backfitting estimators: 
```{r robustfit}
robust.fit <- backf.rob(medv ~ ., data = dd, degree = 0, type = 'Huber', 
                        windows = bandw)
```
Information about the fit can be obtained using the `summary` method: 
```{r summary}
summary(robust.fit)
```
Note that the summary output includes a robust 
version of the R-squared coefficient, which is 
computed as 
$$
R^2_{rob}=\frac{\sum_{i=1}^n \rho\left((Y_i-\widehat{a})/\widehat{\sigma}\right)-\sum_{i=1}^n \rho\left(R_i/\widehat{\sigma}\right)}{\sum_{i=1}^n \rho\left((Y_i-\widehat{a})/\widehat{\sigma}\right)} \, ,
$$
where $\rho$ is the loss function used by the M-kernel smoothers (as determined by the argument `type` in the call to `backf.rob`), and $\widehat{a}$ is a robust location estimator for a model without explanatory variables: $Y=a+\epsilon$.
<!-- are present. $R_i=Y_i-\widehat{\mu}-\sum_{j=1}^p \widehat{g}_j(X_{ij})$ and $\rho$ is the loss function in the `type` argument. -->


The plot method can be used to visualize the estimated 
additive components $\hat{g}_j$, displayed over the corresponding partial
residuals: 
$$
R_{ij}=Y_i-\hat{\mu}-\sum_{k\neq i}\hat{g}_k(X_{ik}) \, , 
\quad 1 \le i \le 506\, , \quad 1 \le j \le 10 \, . 
$$ 
```{r plot, out.width  = "45%"}
plot(robust.fit)
```

By default, `backf.rob` computes fitted values on the training
set. If predictions at a different specific point are desired, 
we can pass those points using the argumen `point`. For example, 
to obtain predicted values at a point `po` given by the average of 
the (log transformed) explanatory variables, we can use the following 
command (note that this step implies re-fitting the whole model):
```{r preds}
po <- colMeans(dd[, names(dd) != 'medv'])
robust.fit1 <- backf.rob(medv ~ ., data = dd, degree = 0, type = 'Huber', 
                         windows = bandw, point = po)
```
The values of the estimated components evaluated at the 
corresponding coordinates of `po` are 
returned in the `$prediction` element: 
```{r showpred}
robust.fit1$prediction
```
In order to illustrate the behaviour of the robust fit when 
outliers are present in the data, we artifically introduce
1% of atypical values in the response variable: 
```{r outliers}
dd2 <- dd
dd2$medv[1:5]<- rep(400, 5)
```
We now calculate the robust estimators using  the contaminated
data set. Note that we expect them  to be very similar to 
those we obtained before with the original 
training set. 
```{r robustplotswithoutliers}
robust.fit.new <- backf.rob(medv ~ ., data = dd2, degree = 0, type = 'Huber', 
                            windows = bandw, point = po)
summary(robust.fit.new)
robust.fit.new$prediction
```
From the output above we verify that the predictions at the point 
`po` with both fits are very similar to each other.  
Because the magnitude of the outliers affects the scale of the 
partial residuals plots, to compare both fits below we plot
the robust estimators for each additive component
trained on the original and 
contaminated data sets (without including the partial residuals). 
In green and dashed lines the robust estimator computed with the original data set, and in blue and solid lines the robust estimator computed with the contaminated data set. We see that, indeed, both sets of 
estimated additive components are very similar to 
each other. 
```{r robustplots2, warning=FALSE, out.width  = "45%"}
for(j in 1:10) {
  name.x <- names(dd)[j] 
  name.y <- bquote(paste(hat('g')[.(j)]))
  oo <- order(dd2[,j])
  plot(dd2[oo,j], robust.fit.new$g.matrix[oo,j], type="l", lwd=2, col='blue', lty=1, 
       xlab=name.x, ylab=name.y)
  lines(dd2[oo,j], robust.fit$g.matrix[oo,j], lwd=2, col='green', lty=2)
}
```

It is easy to see that when we use the classical
backfitting estimator the fits obtained 
when the training set contains outliers
are dramatically different from the ones obtained with the original
data set. We use the package `gam` to compute the 
standard backfitting algorithm (based on local kernel regression smoothers).
The bandwidths used to compute the classical estimates are again proportional to the standard deviations of the explanatory variables,
but we use a slightly larger coefficient to avoid numerical issues 
with the local fits. We set $h_j=(3/4) \, \widehat{\sigma}_j$, 
for $1 \le j \le 10$: 
```{r gam, warning=FALSE}
library(gam)
fit.gam <- gam(medv ~ lo(crim, span=1.62) + lo(indus, span=0.58) + 
                 lo(nox, span=0.15) + lo(rm, span=0.08) +
                 lo(age, span=0.46) + lo(dis, span=0.40) + 
                 lo(tax, span=0.30) + lo(ptratio, span=0.09) +
                 lo(black, span=0.58) + lo(lstat, span=0.45), data=dd)
fits <- predict(fit.gam, type='terms')
fit.gam.new <- gam(medv ~ lo(crim, span=1.62) + lo(indus, span=0.58) + 
                 lo(nox, span=0.15) + lo(rm, span=0.08) +
                 lo(age, span=0.46) + lo(dis, span=0.40) + 
                 lo(tax, span=0.30) + lo(ptratio, span=0.09) +
                 lo(black, span=0.58) + lo(lstat, span=0.45), data=dd2)
fits.new <- predict(fit.gam.new, type='terms')
```
In the plots below, the standard backfitting estimates 
calculated with the original 
`Boston` data set are shown with 
orange and dashed lines, while those obtained with the 
contaminated training set are shown with 
purple and  solid lines. 
```{r gamplots, out.width  = "45%"}
for(j in 1:10) {
  oo <- order(dd2[,j])
  name.x <- names(dd)[j] 
  name.y <- bquote(paste(hat('g')[.(j)]))
  plot(dd2[oo,j], fits.new[oo,j], type="l", lwd=2, col='purple', lty=1, 
       xlab=name.x, ylab=name.y)
  lines(dd2[oo,j], fits[oo,j], lwd=2, col='darkorange2', lty=2)
}
```

## Airquality example

The ``airquality`` data set contains 153 daily air quality measurements in the New York region between May and September, 1973 (Chambers et al., 1983). The interest is in modeling the mean Ozone (\lq\lq $\mbox{O}_3$\rq\rq) concentration as a function of 3 potential
explanatory variables: solar radiance in the frequency band
4000-7700 (\lq\lq Solar.R\rq\rq), wind speed (\lq\lq Wind\rq\rq) and temperature (\lq\lq Temp\rq\rq). We focus on the 111 complete entries in the data set. 
```{r scatterplot, out.width  = "75%"}
data(airquality)
ccs <- complete.cases(airquality)
aircomplete <- airquality[ccs, c('Ozone', 'Solar.R', 'Wind', 'Temp')]
pairs(aircomplete[, c('Ozone', 'Solar.R', 'Wind', 'Temp')], pch=19, col='gray30')
```

The scatter plot suggests that the relationship between ozone and the other variables in not linear and so we propose using an additive regression model of the form 
\begin{equation} \label{eq:ozone-model}
\mbox{Ozone}=\mu+g_{1}(\mbox{Solar.R})+g_{2}(\mbox{Wind})+g_{3}(\mbox{Temp}) + \varepsilon \, .
\end{equation} 
To fit this model above we use robust local linear kernel M-estimators and Tukey's bisquare loss function. These choices are set using the arguments ``degree = 1`` and ``type='Tukey'`` in the call to the function ``backf.rob``.  The model is specified with the standard formula notation in R. The argument ``windows`` is a vector with the bandwidths to be used with each kernel smoother. To estimate optimal values we used a robust leave-one-out cross-validation approach (see Boente et al., 2017). As a robust prediction error measure we use `mu^2 + sigma^2` where `mu` and `sigma` are M-estimators of location and scale of the prediction errors, respectively. 

```{r robustcv, warning=FALSE, eval=FALSE}
library(RBF)

# Bandwidth selection with leave-one-out cross-validation
## Without outliers
# This takes a long time to compute (approx 380 minutes running
# R 3.6.1 on an Intel(R) Core(TM) i7-4790 CPU @ 3.60GHz)
a <- c(1/2, 1, 1.5, 2, 2.5, 3)
h1 <- a * sd(aircomplete[,2])
h2 <- a * sd(aircomplete[,3])
h3 <- a * sd(aircomplete[,4])
hh <- expand.grid(h1, h2, h3)
nh <- nrow(hh)
rmspe <- rep(NA, nh)
jbest <- 0
cvbest <- +Inf
n <- nrow(aircomplete)
for(i in 1:nh) {
  # leave-one-out CV loop
  preds <- rep(NA, n)
  for(j in 1:n) {
    tmp <- try( backf.rob(Ozone ~ Solar.R + Wind + Temp, point = aircomplete[j, -1],
                          windows = hh[i, ], epsilon = 1e-6, data = aircomplete,
                          degree = 1, type = 'Tukey', subset = c(-j) ))
    if (class(tmp)[1] != "try-error") {
      preds[j] <- rowSums(tmp$prediction) + tmp$alpha
    }
  }
  tmp.re <- RobStatTM::locScaleM(preds - aircomplete$Ozone, na.rm=TRUE)
  rmspe[i] <- tmp.re$mu^2 + tmp.re$disper^2
  if( rmspe[i] < cvbest ) {
    jbest <- i
    cvbest <- rmspe[i]
  }
}
(bandw <- hh[jbest,])
```
The resulting bandwidths are:
```{r bandw}
bandw <- c(136.7285, 10.67314, 4.764985)
```
Now we use the robust backfitting algorithm to fit an additive model using Tukey's bisquare loss (the default tuning constant for this loss function is 4.685) and the optimal bandwidths.

```{r fitfull}
fit.full <- backf.rob(Ozone ~ Solar.R + Wind + Temp, windows = bandw, 
                      epsilon = 1e-6, degree = 1, type = 'Tukey', 
                      subset = ccs, data = airquality)
```
We can visually explore the estimated additive functions plotted 
over the corresponding partial residuals using the method `plot`:
```{r plotfitfull, out.width  = "45%"}
plot(fit.full)
```

As before, we use the R package `gam` to compute the 
classical additive model estimators for this model. 
Optimal bandwidths were calculated using 
leave-one-out cross-validation as before:
```{r gamcv, warning=FALSE, eval=FALSE}
library(gam)
a <- c(.3, .4, .5, .6, .7, .8, .9)
hh <- expand.grid(a, a, a)
nh <- nrow(hh)
jbest <- 0
cvbest <- +Inf
n <- nrow(aircomplete)
for(i in 1:nh) {
  fi <- rep(0, n)
  for(j in 1:n) {
    tmp <- gam(Ozone ~ lo(Solar.R, span=hh[i,1]) + lo(Wind, span=hh[i,2])
               + lo(Temp, span=hh[i,3]), data = aircomplete, subset=c(-j))
    fi[j] <- as.numeric(predict(tmp, newdata=aircomplete[j, -1], type='response'))
  }
  ss <- mean((aircomplete$Ozone - fi)^2)
  if(ss < cvbest) {
    jbest <- i
    cvbest <- ss
  }
}
(hh[jbest,])
# Var1 Var2 Var3
# 0.7  0.7  0.5
```
The optimal bandwidths are `0.7`, `0.7` and `0.5` for `Solar.R`, `Wind` and `Temp`, respectively, and we use them to compute the backfitting estimators: 
```{r fitgam}
fit.gam <- gam(Ozone ~ lo(Solar.R, span=.7) + lo(Wind, span=.7)+
                 lo(Temp, span=.5), data = aircomplete)
```
Both classical (in magenta and dashed lines) and robust (in blue and solid lines) fits are shown in the following plot together with the partial residuals obtained by the robust fit.
```{r plotrobgam, out.width  = "45%"}
x <- as.matrix( aircomplete[ , c('Solar.R', 'Wind', 'Temp')] )
y <- as.vector( aircomplete[ , 'Ozone'] )
fits <- predict(fit.gam, type='terms')
for(j in 1:3) {
  re <- fit.full$yp - fit.full$alpha - rowSums(fit.full$g.matrix[,-j])
  plot(re ~ x[,j], type='p', pch=20, col='gray45', xlab=colnames(x)[j], ylab='')
  oo <- order(x[,j])
  lines(x[oo,j], fit.full$g.matrix[oo,j], lwd=2, col='blue', lty=1)
  lines(x[oo,j], fits[oo,j], lwd=2, col='magenta', lty=2)
}
```

The two fits differ mainly on the estimated effects of wind speed and temperature. The classical estimate for $g_1(\mbox{Temp})$ is consistently lower than the robust counterpart for $\mbox{Temp} \ge 85$. For wind speed,
the non-robust estimate $\hat{g}_2(\mbox{Wind})$ suggests a higher effect over Ozone concentrations for low wind speeds than the one given by the robust estimate, and the opposite difference for higher speeds.

Since residuals from a robust fit can generally be used to detect the presence of atypical observations in the training data, we plot the boxplot of the residuals obtained by the robust fit and 4 possible outlying points (indicated with red circles) can be observed.

```{r boxplot, out.width  = "50%"}
re.ro <- residuals(fit.full)
ou.ro <- boxplot(re.ro, col='gray80')$out
in.ro <- (1:length(re.ro))[ re.ro %in% ou.ro ]
points(rep(1, length(in.ro)), re.ro[in.ro], pch=20, col='red')
(in.ro)
```

We highlight these suspicious observations on the
scatter plot.

```{r scatterplotpoints, out.width  = "75%"}
cs <- rep('gray30', nrow(aircomplete))
cs[in.ro] <- 'red'
os <- 1:nrow(aircomplete)
os2 <- c(os[-in.ro], os[in.ro])
pairs(aircomplete[os2, c('Ozone', 'Solar.R', 'Wind', 'Temp')], 
      pch=19, col=cs[os2])
```


Note that not all these suspected atypical observations are particularly extreme, or directly evident on the scatter plot. However, as we will show below, they do have an important effect on the estimates of the components of the additive model. 

The partial residuals corresponding to these points can be also visualized in red in the plot of the estimated curves.
```{r plotoutred, out.width  = "45%"}
# Plot both fits (robust and classical) 
x <- as.matrix( aircomplete[ , c('Solar.R', 'Wind', 'Temp')] )
y <- as.vector( aircomplete[ , 'Ozone'] )
fits <- predict(fit.gam, type='terms')
for(j in 1:3) {
  re <- fit.full$yp - fit.full$alpha - rowSums(fit.full$g.matrix[,-j])
  plot(re ~ x[,j], type='p', pch=20, col='gray45', xlab=colnames(x)[j], ylab='')
  points(re[in.ro] ~ x[in.ro,j], pch=20, col='red')
  oo <- order(x[,j])
  lines(x[oo,j], fit.full$g.matrix[oo,j], lwd=2, col='blue', lty=1)
  lines(x[oo,j], fits[oo,j], lwd=2, col='magenta', lty=2)
}
```

To investigate whether the differences between the robust and non-robust estimators are due to the outliers, we recomputed the classical fit after removing them.

We ran a similar leave-one-out cross-validation experiment to select the spans for each the 3 univariate smoothers.

```{r cvgamclean, eval=FALSE}
airclean <- aircomplete[-in.ro, c('Ozone', 'Solar.R', 'Wind', 'Temp')]
a <- c(.3, .4, .5, .6, .7, .8, .9)
hh <- expand.grid(a, a, a)
nh <- nrow(hh)
jbest <- 0
cvbest <- +Inf
n <- nrow(airclean)
for(i in 1:nh) {
  fi <- rep(0, n)
  for(j in 1:n) {
    tmp <- gam(Ozone ~ lo(Solar.R, span=hh[i,1]) + lo(Wind, span=hh[i,2])
               + lo(Temp, span=hh[i,3]), data=airclean, subset=c(-j))
    fi[j] <- as.numeric(predict(tmp, newdata=airclean[j,], type='response'))
  }
  ss <- mean((airclean$Ozone - fi)^2)
  if(ss < cvbest) {
    jbest <- i
    cvbest <- ss
  }
}
(hh[jbest,])
# Var1 Var2 Var3
# 0.7  0.8  0.3
```

We use the optimal bandwidths to compute non-robust fit.
```{r fitgam2}
airclean <- aircomplete[-in.ro, c('Ozone', 'Solar.R', 'Wind', 'Temp')]
fit.gam2 <- gam(Ozone ~ lo(Solar.R, span=.7) + lo(Wind, span=.8)+
                  lo(Temp, span=.3), data=airclean) 
```

The following plot shows the estimated curves obtained with the classical estimator using the \lq\lq clean\rq\rq\ data together with the robust ones (computed on the whole data set).  Outliers are highlighted in red. 

```{r finalplot, out.width  = "45%"}
fits2 <- predict(fit.gam2, type='terms')
dd2 <- aircomplete[-in.ro, c('Solar.R', 'Wind', 'Temp')]
for(j in 1:3) {
  re <- fit.full$yp - fit.full$alpha - rowSums(fit.full$g.matrix[,-j])
  plot(re ~ x[,j], type='p', pch=20, col='gray45', xlab=colnames(x)[j], ylab='')
  points(re[in.ro] ~ x[in.ro,j], pch=20, col='red')
  oo <- order(dd2[,j])
  lines(dd2[oo,j], fits2[oo,j], lwd=2, col='magenta', lty=2)
  oo <- order(x[,j])
  lines(x[oo,j], fit.full$g.matrix[oo,j], lwd=2, col='blue', lty=1)
}
```

Note that both fits are now very close. An intuitive interpretation is that the robust fit has automatically down-weighted potential outliers and produced estimates very similar to the classical ones applied to the \lq\lq clean\rq\rq\ observations.


### Prediction comparison

Finally, we compare the prediction accuracy obtained with each
of the fits. Because we are not interested in predicting well any
possible outliers in the data, we evaluate the quality of
the predictions using a 5%-trimmed mean 
squared prediction error (effectively measuring the prediction 
accuracy on 95% of the data). We use this alpha-trimmed mean 
squared function:
```{r pred1}
tms <- function(a, alpha=.1) {
  # alpha is the proportion to trim
  a2 <- sort(a^2, na.last=NA)
  n0 <- floor( length(a) * (1 - alpha) )
  return( mean(a2[1:n0], na.rm=TRUE) )
}
```

We use 100 runs of 5-fold CV to compare the 
5%-trimmed mean squared prediction error of
the robust fit and the classical one. 
Note that the bandwidths are kept fixed 
at their optimal value estimated above. 
```{r pred2, warning=FALSE, eval=FALSE}
dd <- airquality
dd <- dd[complete.cases(dd), c('Ozone', 'Solar.R', 'Wind', 'Temp')]
# 100 runs of K-fold CV
M <- 100
# 5-fold
K <- 5
n <- nrow(dd)
# store (trimmed) TMSPE for robust and gam, and also
tmspe.ro <- tmspe.gam <- vector('numeric', M)
set.seed(123)
ii <- (1:n)%%K + 1
for(runs in 1:M) {
  tmpro <- tmpgam <- vector('numeric', n)
  ii <- sample(ii)
  for(j in 1:K) {
    fit.full <- backf.rob(Ozone ~ Solar.R + Wind + Temp, 
                           point=dd[ii==j, -1], windows = bandw, 
                           epsilon = 1e-6, degree = 1, type = 'Tukey', 
                           subset = (ii!=j), data = dd)
    tmpro[ ii == j ] <- rowSums(fit.full$prediction) + fit.full$alpha
    fit.gam <- gam(Ozone ~ lo(Solar.R, span=.7) + lo(Wind, span=.7)+
                     lo(Temp, span=.5), data = dd[ii!=j, ])
    tmpgam[ ii == j ] <- predict(fit.gam, newdata=dd[ii==j, ], type='response')
  }
  tmspe.ro[runs] <- tms( dd$Ozone - tmpro, alpha=0.05)
  tmspe.gam[runs] <- tms( dd$Ozone - tmpgam, alpha=0.05)
}
```
These are the boxplots. We see that the robust fit consistently fits the vast majority (95%) of the data better than the classical one. 

![Boxplots of the prediction errors \label{fig:pred3}](pred3.png){ width=50% }

<!-- As a sanity check, we compare the prediction accuracy of the robust and non-robust fits using only the "clean" data set. We re-compute the optimal bandwidths for the robust fit using leave-one-out cross validation as above, and as above, note that these bandwidths are kept fixed.  We use 100 runs of 5-fold cross-validation, and compute both the trimmed and the regular mean squared prediction errors of each fit. -->
<!-- ```{r pred.clean, cache=TRUE, warning=FALSE, eval=FALSE} -->
<!-- aq <- airquality -->
<!-- aq2 <- aq[complete.cases(aq), c('Ozone', 'Solar.R', 'Wind', 'Temp')] -->
<!-- airclean <- aq2[ -in.ro, ] -->
<!-- bandw <- c(138.2699, 10.46753, 4.828436) -->
<!-- M <- 100  -->
<!-- K <- 5 -->
<!-- n <- nrow(airclean) -->
<!-- mspe.ro <- mspe.gam <- tmspe.ro <- tmspe.gam <- vector('numeric', M) -->
<!-- set.seed(17) -->
<!-- ii <- (1:n)%%K + 1 -->
<!-- for(runs in 1:M) { -->
<!--   tmpro <- tmpgam <- vector('numeric', n) -->
<!--   ii <- sample(ii) -->
<!--   for(j in 1:K) { -->
<!--     fit.full <- try( backf.rob(Ozone ~ Solar.R + Wind + Temp,  -->
<!--                           point=airclean[ii==j, -1], windows = bandw,  -->
<!--                           epsilon = 1e-6, degree = 1, type = 'Tukey',  -->
<!--                           subset = (ii!=j), data = airclean) ) -->
<!--     if (class(fit.full)[1] != "try-error") { -->
<!--       tmpro[ ii == j ] <- rowSums(fit.full$prediction) + fit.full$alpha -->
<!--     } -->
<!--     fit.gam <- gam(Ozone ~ lo(Solar.R, span=.7) + lo(Wind, span=.8)+ -->
<!--                      lo(Temp, span=.3), data = airclean, subset = (ii!=j) ) -->
<!--     tmpgam[ ii == j ] <- predict(fit.gam, newdata=airclean[ii==j, ],  -->
<!--                                  type='response') -->
<!--   } -->
<!--   tmspe.ro[runs] <- tms( airclean$Ozone - tmpro, alpha=0.05) -->
<!--   mspe.ro[runs] <- mean( ( airclean$Ozone - tmpro)^2, na.rm=TRUE) -->
<!--   tmspe.gam[runs] <- tms( airclean$Ozone - tmpgam, alpha=0.05) -->
<!--   mspe.gam[runs] <- mean( ( airclean$Ozone - tmpgam)^2, na.rm=TRUE) -->
<!-- } -->
<!-- ``` -->
<!-- The boxplots of the trimmed and regular mean squared prediction  -->
<!-- errors over the 100 cross-validation runs are below. We see that -->
<!-- for the majority of the runs both estimators provide very similar -->
<!-- prediction errors. Note that we naturally expect a robust method to  -->
<!-- perform slightly worse than the classical one when no model -->
<!-- deviations occur. The boxplots below show that for this  -->
<!-- robust backfitting estimator, this loss in prediction accuracy is -->
<!-- in fact very small.  -->


<!-- ![Trimmed MSPE On \lq\lq clean\rq\rq data \label{fig:tmspeclean}](tmspeclean.png){ width=50% } -->





<!-- ![Non-Trimmed MSPE On \lq\lq clean\rq\rq data \label{fig:nontmspeclean}](nontmspeclean.png){ width=50% } -->



## Bibliography

Boente G, Martinez A, and Salibian-Barrera M. (2017) Robust estimators for additive models using backfitting. *Journal of Nonparametric Statistics*, **29**, 744-767. [DOI: 10.1080/10485252.2017.1369077](https://doi.org/10.1080/10485252.2017.1369077)

Chambers, J. M., Cleveland W. S., Kleiner B. and Tukey A. (1983). *Graphical Methods for Data Analysis*. 2nd edition. Chapman \& Hall. [DOI: 10.1201/9781351072304](https://doi.org/10.1201/9781351072304)

Härdle, W., Müller, M., Sperlich, S. and Werwatz, A. (2004). *Nonparametric and Semiparametric Models*. Springer. [DOI: 10.1007/978-3-642-17146-8](https://doi.org/10.1007/978-3-642-17146-8)

Harrinson, D. and Rubinfeld, D. L. (1978). Hedonic prices and the demand for clean air. *J. Environ. Economics and Management*, **5**, 81-102. [DOI: 10.1016/0095-0696(78)90006-2](https://doi.org/10.1016/0095-0696(78)90006-2)

Martinez A, and Salibian-Barrera M. (2021) RBF: An R package to compute a robust backfitting estimator for additive models. *Journal of Open Source Software*, **6(60)**. [DOI: 10.21105/joss.02992](https://doi.org/10.21105/joss.02992)
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RBF-fn.R
\name{backf.cl}
\alias{backf.cl}
\title{Classic Backfitting}
\usage{
backf.cl(
  formula,
  data,
  subset,
  point = NULL,
  windows,
  epsilon = 1e-06,
  degree = 0,
  prob = NULL,
  max.it = 100
)
}
\arguments{
\item{formula}{an object of class \code{formula} (or one that can be coerced to 
that class): a symbolic description of the model to be fitted.}

\item{data}{an optional data frame, list or environment (or object coercible 
by \link{as.data.frame} to a data frame) containing the variables in the model. 
If not found in \code{data}, the variables are taken from \code{environment(formula)}, 
typically the environment from which the function was called.}

\item{subset}{an optional vector specifying a subset of observations to be used in 
the fitting process.}

\item{point}{matrix of points where predictions will be computed and returned.}

\item{windows}{vector of bandwidths for the local polynomial smoother,
one per explanatory variable.}

\item{epsilon}{convergence criterion. Maximum allowed relative difference between
consecutive estimates}

\item{degree}{degree of the local polynomial smoother. Defaults to \code{0} (local constant).}

\item{prob}{vector of probabilities of observing each response (length n).
Defaults to \code{NULL} and in that case it is ignored.}

\item{max.it}{Maximum number of iterations for the algorithm.}
}
\value{
A list with the following components:
\item{alpha}{Estimate for the intercept.}
\item{g.matrix }{Matrix of estimated additive components (n by p).}
\item{prediction }{Matrix of estimated additive components for the points listed in
the argument \code{point}.}
}
\description{
This function computes the standard backfitting algorithm for additive models.
}
\details{
This function computes the standard backfitting algorithm for additive models,
using a squared loss function and local polynomial smoothers.
}
\examples{
data(airquality)
tmp <- backf.cl(Ozone ~ Solar.R + Wind + Temp, data=airquality, 
subset=complete.cases(airquality), windows=c(130, 9, 10), degree=1)

}
\references{
Hasie, TJ and Tibshirani, RJ. Generalized Additive Models, 1990. Chapman
and Hall, London.
}
\author{
Matias Salibian-Barrera, \email{matias@stat.ubc.ca}, Alejandra Martinez
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RBF-fn.R
\name{summary.backf}
\alias{summary.backf}
\alias{summary.backf.cl}
\alias{summary.backf.rob}
\title{Summary for additive models fits using backfitting}
\usage{
\method{summary}{backf}(object, ...)
}
\arguments{
\item{object}{an object of class \code{backf}, a result of a call to
\code{\link{backf.cl}} or \code{\link{backf.rob}}.}

\item{...}{additional other arguments. Currently ignored.}
}
\description{
Summary method for class \code{backf}.
}
\details{
This function returns the estimation of the intercept and also the
five-number summary and the mean of the residuals for both classical and
robust estimators. For the classical estimator, it also returns the R-squared.
For the robust estimator it returns a robust version of the R-squared and 
the estimate of the residual standard error.
}
\author{
Alejandra Mercedes Martinez \email{ale_m_martinez@hotmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RBF-fn.R
\name{formula.backf}
\alias{formula.backf}
\title{Additive model formula}
\usage{
\method{formula}{backf}(x, ...)
}
\arguments{
\item{x}{an object of class \code{backf}, a result of a call to \code{\link{backf.cl}} or \code{\link{backf.rob}}.}

\item{...}{additional other arguments. Currently ignored.}
}
\value{
A model formula.
}
\description{
Description of the additive model formula extracted from an object of class \code{backf}.
}
\author{
Alejandra Mercedes Martinez \email{ale_m_martinez@hotmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RBF-fn.R
\name{predict.backf}
\alias{predict.backf}
\title{Fitted values for objects of class \code{backf}.}
\usage{
\method{predict}{backf}(object, ...)
}
\arguments{
\item{object}{an object of class \code{backf}, a result of a call to \code{\link{backf.cl}} or \code{\link{backf.rob}}.}

\item{...}{additional other arguments. Currently ignored.}
}
\value{
A vector of fitted values.
}
\description{
This function returns the fitted values given the covariates of
the original sample under an additive model using the classical or
robust backfitting approach computed with \code{\link{backf.cl}} or
\code{\link{backf.rob}}.
}
\author{
Alejandra Mercedes Martinez \email{ale_m_martinez@hotmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RBF-fn.R
\name{fitted.values.backf}
\alias{fitted.values.backf}
\title{Fitted values for objects of class \code{backf}}
\usage{
fitted.values.backf(object, ...)
}
\arguments{
\item{object}{an object of class \code{backf}, a result of a call to \code{\link{backf.cl}} or \code{\link{backf.rob}}.}

\item{...}{additional other arguments. Currently ignored.}
}
\value{
A vector of fitted values.
}
\description{
This function returns the fitted values given the covariates of the original sample under an additive model using a classical or robust marginal integration procedure estimator computed with \code{backf.cl} or \code{backf.rob}.
}
\author{
Alejandra Mercedes Martinez \email{ale_m_martinez@hotmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RBF-fn.R
\name{residuals.backf}
\alias{residuals.backf}
\title{Residuals for objects of class \code{backf}}
\usage{
\method{residuals}{backf}(object, ...)
}
\arguments{
\item{object}{an object of class \code{backf}, a result of a call to \code{\link{backf.cl}} or \code{\link{backf.rob}}.}

\item{...}{additional other arguments. Currently ignored.}
}
\value{
A vector of residuals.
}
\description{
This function returns the residuals of the fitted additive model using
the classical or robust backfitting estimators, as computed with \code{\link{backf.cl}} or
\code{\link{backf.rob}}.
}
\author{
Alejandra Mercedes Martinez \email{ale_m_martinez@hotmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RBF-fn.R
\name{psi.huber}
\alias{psi.huber}
\title{Derivative of Huber's loss function.}
\usage{
psi.huber(r, k = 1.345)
}
\arguments{
\item{r}{a vector of real numbers}

\item{k}{a positive tuning constant.}
}
\value{
A vector of the same length as \code{x}.
}
\description{
This function evaluates the first derivative of Huber's loss function.
}
\details{
This function evaluates the first derivative of Huber's loss function.
}
\examples{
x <- seq(-2, 2, length=10)
psi.huber(r=x, k = 1.5)

}
\author{
Matias Salibian-Barrera, \email{matias@stat.ubc.ca}, Alejandra Martinez
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RBF-fn.R
\name{backf.rob}
\alias{backf.rob}
\title{Robust Backfitting}
\usage{
backf.rob(
  formula,
  data,
  subset,
  windows,
  point = NULL,
  epsilon = 1e-06,
  degree = 0,
  sigma.hat = NULL,
  prob = NULL,
  max.it = 50,
  k.h = 1.345,
  k.t = 4.685,
  type = "Huber"
)
}
\arguments{
\item{formula}{an object of class \code{formula} (or one that can be coerced to 
that class): a symbolic description of the model to be fitted.}

\item{data}{an optional data frame, list or environment (or object coercible 
by \link{as.data.frame} to a data frame) containing the variables in the model. 
If not found in \code{data}, the variables are taken from \code{environment(formula)}, 
typically the environment from which the function was called.}

\item{subset}{an optional vector specifying a subset of observations to be used in 
the fitting process.}

\item{windows}{vector of bandwidths for the local polynomial smoother,
one per explanatory variable.}

\item{point}{matrix of points where predictions will be computed and returned.}

\item{epsilon}{convergence criterion. Maximum allowed relative difference between
consecutive estimates}

\item{degree}{degree of the local polynomial smoother. Defaults to \code{0} (local constant).}

\item{sigma.hat}{estimate of the residual standard error. If \code{NULL} (default) we use the
\link{mad} of the residuals obtained with local medians.}

\item{prob}{vector of probabilities of observing each response (length n).
Defaults to \code{NULL} and in that case it is ignored.}

\item{max.it}{Maximum number of iterations for the algorithm.}

\item{k.h}{tuning constant for a Huber-type loss function.}

\item{k.t}{tuning constant for a Tukey-type loss function.}

\item{type}{one of either \code{'Tukey'} or \code{'Huber'}.}
}
\value{
A list with the following components:
\item{alpha}{Estimate for the intercept.}
\item{g.matrix }{Matrix of estimated additive components (n by p).}
\item{prediction }{Matrix of estimated additive components for the points listed in
the argument \code{point}.}
\item{sigma.hat }{Estimate of the residual standard error.}
}
\description{
This function computes a robust backfitting algorithm for additive models
}
\details{
This function computes a robust backfitting algorithm for additive models
using robust local polynomial smoothers.
}
\examples{
data(airquality)
tmp <- backf.rob(Ozone ~ Solar.R + Wind + Temp, data=airquality, 
subset=complete.cases(airquality), windows=c(136.7, 8.9, 4.8), degree=1)

}
\references{
Boente G, Martinez A, Salibian-Barrera M. Robust estimators
for additive models using backfitting. Journal of Nonparametric Statistics,
2017; 29:744-767. https://doi.org/10.1080/10485252.2017.1369077
}
\author{
Matias Salibian-Barrera, \email{matias@stat.ubc.ca}, Alejandra Martinez
}
\name{RBF-package}
\alias{RBF-package}
\alias{RBF}
\docType{package}
\title{
A robust backfitting algorithm for additive models.
}
\description{
A robust backfitting algorithm for additive models.
}
\details{
\tabular{ll}{
Package: \tab RBF\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-01-19\cr
License: \tab GPL 3.0\cr
}
}
\author{
Matias Salibian-Barrera, Alejandra Martinez

Maintainer: Matias Salibian-Barrera <matias@stat.ubc.ca>
}
\references{
Boente G, Martinez A, Salibian-Barrera M. Robust estimators for additive models using backfitting. Journal of Nonparametric Statistics, 2017; 29:744-767. https://doi.org/10.1080/10485252.2017.1369077
}
\keyword{ Robust Backfitting }
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RBF-fn.R
\name{plot.backf}
\alias{plot.backf}
\title{Diagnostic plots for objects of class \code{backf}}
\usage{
\method{plot}{backf}(x, ask = FALSE, which = 1:np, ...)
}
\arguments{
\item{x}{an object of class \code{backf}, a result of a call to \code{\link{backf.cl}} or \code{\link{backf.rob}}.}

\item{ask}{logical value. If \code{TRUE}, the graphical device will prompt for confirmation before
going to the next page/screen of output.}

\item{which}{vector of indices of explanatory variables for which partial residuals plots will
be generaetd. Defaults to all available explanatory variables.}

\item{...}{additional other arguments. Currently ignored.}
}
\description{
Plot method for objects of class \code{backf}.
}
\examples{
tmp <- backf.rob(Ozone ~ Solar.R + Wind + Temp, data=airquality, 
subset=complete.cases(airquality), windows=c(136.7, 8.9, 4.8), degree=1)
plot(tmp, which=1:2)

}
\author{
Alejandra Mercedes Martinez \email{ale_m_martinez@hotmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RBF-fn.R
\name{psi.tukey}
\alias{psi.tukey}
\title{Derivative of Tukey's bi-square loss function.}
\usage{
psi.tukey(r, k = 4.685)
}
\arguments{
\item{r}{a vector of real numbers}

\item{k}{a positive tuning constant.}
}
\value{
A vector of the same length as \code{x}.
}
\description{
This function evaluates the first derivative of Tukey's bi-square loss function.
}
\details{
This function evaluates the first derivative of Tukey's bi-square loss function.
}
\examples{
x <- seq(-2, 2, length=10)
psi.tukey(r=x, k = 1.5)

}
\author{
Matias Salibian-Barrera, \email{matias@stat.ubc.ca}, Alejandra Martinez
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RBF-fn.R
\name{deviance.backf}
\alias{deviance.backf}
\title{Deviance for objects of class \code{backf}}
\usage{
\method{deviance}{backf}(object, ...)
}
\arguments{
\item{object}{an object of class \code{backf}, a result of a call to \code{\link{backf.cl}} or \code{\link{backf.rob}}.}

\item{...}{additional other arguments. Currently ignored.}
}
\value{
A real number.
}
\description{
This function returns the deviance of the fitted additive model using one of the three
classical or robust marginal integration estimators, as computed with \code{\link{backf.cl}} or
\code{\link{backf.rob}}.
}
\author{
Alejandra Mercedes Martinez \email{ale_m_martinez@hotmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RBF-fn.R
\name{print.backf}
\alias{print.backf}
\title{Print a Marginal Integration procedure}
\usage{
\method{print}{backf}(x, ...)
}
\arguments{
\item{x}{an object of class \code{backf}, a result of a call to \code{\link{backf.cl}} or \code{\link{backf.rob}}.}

\item{...}{additional other arguments. Currently ignored.}
}
\value{
A real number.
}
\description{
The default print method for a \code{backf} object.
}
\author{
Alejandra Mercedes Martinez \email{ale_m_martinez@hotmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RBF-fn.R
\name{k.epan}
\alias{k.epan}
\title{Epanechnikov kernel}
\usage{
k.epan(x)
}
\arguments{
\item{x}{a vector of real numbers}
}
\value{
A vector of the same length as \code{x} where each entry is
\code{0.75 * (1 - x^2)} if \code{x < 1} and 0 otherwise.
}
\description{
This function evaluates an Epanechnikov kernel
}
\details{
This function evaluates an Epanechnikov kernel
}
\examples{
x <- seq(-2, 2, length=10)
k.epan(x)

}
\author{
Matias Salibian-Barrera, \email{matias@stat.ubc.ca}, Alejandra Martinez
}
