---
title: "figure1"
output: html_document
---

# setup
```{r knitr_setup, cache=FALSE, echo=FALSE, results='hide'}
library(knitr)
cachedir <- 'cache/'
opts_chunk$set(cache=TRUE, autodep=TRUE, dev=c('png', 'cairo_pdf', 'svglite'), dpi=200, fig.path='plots/', fig.width=5, fig.height=2.5, dev.args=list(system_fonts=list(sans='Helvetica')), cache.path=cachedir, fig.keep='high', tidy='styler', warning=FALSE, message=FALSE)
options(width=110)
```

```{r setup, cache=FALSE, message=FALSE}
library(tibble)
library(dplyr)
library(tidyr)
library(magrittr)
library(purrr)
library(ggplot2)
library(ggsci)
library(reticulate)
library(arrow)

th <- theme_bw()
th$line$size <- 0.5
th$rect$size <- 0.5
th$panel.grid.minor <- element_blank()
th$panel.border$colour <- 'black'
th$axis.ticks$colour <- "black"
th$axis.text <- element_text()
th$legend.text <- element_text()
th$strip.text['colour'] <- list(NULL)
th$strip.text['size'] <- list(NULL)
th$plot.title$hjust <- 0.5
th$plot.title['size'] <- list(NULL)
th$aspect.ratio <- 1
theme_set(th)

guide_fcolorbar <- function(...)guide_colourbar(frame.colour='black', ticks.colour='black')
scale_color_discrete <- scale_colour_discrete <- function(...) scale_color_startrek(...)# scale_color_brewer(..., palette='Set1')
scale_fill_discrete <- function(...) scale_fill_startrek(...)#scale_fill_brewer(..., palette='Set1')

update_geom_defaults("point", list(shape=16, size=1.5))
update_geom_defaults("path", list(size=0.5))
update_geom_defaults('smooth', list(size=0.5, color='black'))
update_geom_defaults('rect', list(size=0.5))
update_geom_defaults("text", list(size=th$text$size / .pt))
update_geom_defaults("violin", list(color='black'))
update_geom_defaults("boxplot", list(color='black'))
#update_geom_defaults("text_repel", list(size=th$text$size / .pt))

set.seed(42)
```
```{python score_test_simulation, message=FALSE, results="hide"}
import sys
sys.path.append("/home/kats/SpatialDE")
import SpatialDE

import numpy as np
import scipy
import pandas as pd
import anndata as ad

import gpflow

rng = np.random.default_rng(seed=42)

npoints_per_dim = 25
npoints = npoints_per_dim ** 2
coords = (
    np.stack(tuple(i.reshape(-1) for i in np.meshgrid(np.arange(npoints_per_dim), np.arange(npoints_per_dim), indexing="ij")), axis=-1) +
    rng.normal(scale=0.05, size=(npoints, 2))
)

nsamples = 10000
dispersion = 0.5
m = 5
sigma = 0.1

sizefactors = 1 + rng.random(npoints) * 50
means = m * sizefactors

var = means + dispersion * means**2
r = np.repeat(1/dispersion, npoints)
p = 1 - (var - means) / var

samples_null = rng.negative_binomial(np.repeat(r[:, np.newaxis], nsamples, axis=1), p[:, np.newaxis])

adata_null = ad.AnnData(samples_null)
adata_null.obsm["spatial"] = coords

results_null_cauchy = SpatialDE.test(adata_null)
results_null_omnibus = SpatialDE.test(adata_null, omnibus=True)
```
```{python score_test_simulation_alt, echo=FALSE, eval=FALSE}
distances = scipy.spatial.distance_matrix(coords, coords)
max_distance = np.max(distances)
np.fill_diagonal(distances, np.inf)
min_distance = np.min(distances)

lengthscales = np.logspace(np.log10(2 * min_distance), np.log10(max_distance), 5)
kernels = [gpflow.kernels.SquaredExponential(lengthscales=l) for l in lengthscales]

kernelmats = [k.K(coords).numpy() for k in kernels]

samples_alt = np.full((npoints, nsamples), np.log(m))
for i in range(nsamples):
samples_alt[:,i] += rng.multivariate_normal(np.zeros(npoints), sigma * kernelmats[rng.choice(len(kernels))])
samples_alt = np.exp(samples_alt) * sizefactors[:, np.newaxis]

var = samples_alt + dispersion * samples_alt**2
r = samples_alt**2 / (var - samples_alt)
p = 1 - (var - samples_alt) / var

samples_alt = rng.negative_binomial(r, p)

adata_alt = ad.AnnData(samples_alt)
adata_alt.obsm["spatial"] = coords

results_alt_cauchy = SpatialDE.test(adata_alt, kernel_space={"SE": lengthscales})
results_alt_omnibus = SpatialDE.test(adata_alt, omnibus=True, kernel_space={"SE": lengthscales})
```

```{r get_python_data, dependson="score_test_simulation"}
simulation <- bind_rows(mutate(py$results_null_cauchy[[1]], method="cauchy"),
                        mutate(py$results_null_omnibus[[1]], method="omnibus"))
```
```{r figure_2a, fig.height=2.2, fig.width=5}
group_by(simulation, method) %>%
    mutate(expected=row_number(pval) / n()) %>%
    ggplot(aes(expected, pval, color=method)) +
        geom_abline(slope=1, color="lightgray") +
        geom_point() +
        labs(x="expected p-value", y="observed p-value", color=NULL) +
        scale_x_log10(limits=c(6e-6, 1), expand=expansion(), labels=scales::label_math(format=log10)) +
        scale_y_log10(limits=c(6e-6, 1), expand=expansion(), labels=scales::label_math(format=log10)) +
        guides(color=guide_legend(override.aes=list(size=3))) +
        scale_color_startrek(labels=c(cauchy="Cauchy combination", omnibus="omnibus"))
```


```{r read_data}
ad <- import("anndata")

make_benchmark_summary <- function(bench) {
    group_by_at(bench, vars(-pval, -p.adj, -gene)) %>%
        group_modify(function(.x, .y) {
            map_dfr(seq(0, 1, 0.01), function(nfdr) {
                idx <- which(.x$gene <= 1000)
                power <- sum(.x$p.adj[idx] < nfdr) / length(idx)
                fdr <- sum(.x$p.adj[-idx] < nfdr) / sum(.x$p.adj < nfdr)
                tibble(nominal_fdr=nfdr, fdr=fdr, power=power)
            })
        }) %>%
        ungroup()
}

method_names <- . %>%
    mutate(method_name=recode(method, omnibus="omnibus",
                                      cauchy="cauchy combination",
                                      spark="SPARK",
                                      spark_nocheckpositive="SPARK (no eigenvalue clipping)",
                                      spatialde1="SpatialDE1"))
bench <- readRDS("benchmarks/sparsim_benchmark_rep_1.rds") %>%
    method_names() %>%
    mutate(foldchange=2L * as.integer(foldchange), lengthscaleidx=as.integer(lengthscale), p.adj=if_else(is.na(padj), p.adj, padj)) %>%
    select(-padj, -lengthscale) %>%
    as_tibble()

bench <- pmap_dfr(list(unique(bench$lengthscaleidx)), function(lengthscaleidx) {
    adata <- ad$read_h5ad(file.path("simulated", paste0("c2l_mouse_brain_lengthscale_", lengthscaleidx, "_fc_2_rep_1.h5ad")), backed="r")
    tibble(lengthscale=adata$uns[["simulation_params"]][["lengthscale"]])
}, .id="lengthscaleidx") %>%
    mutate(lengthscaleidx=as.integer(lengthscaleidx)) %>%
    inner_join(bench)

bench_summary_default <- make_benchmark_summary(bench)
```
```{r figure2b, fig.width=6, fig.height=2.15}
filter(bench_summary_default, nominal_fdr == 0.05, foldchange == 4) %>%
    ggplot(aes(lengthscale, power, color=method_name)) +
        geom_line() +
        #geom_point() +
        scale_y_continuous(limits=c(0,1), expand=expansion()) +
        guides(color=guide_legend(title=NULL, override.aes=list(size=2)))
```

```{r figure2c, fig.width=11, fig.height=2.4}
single <- readRDS("tests_runtime_1cores.rds") %>%
    mutate(cores="1 core")
multi <- readRDS("tests_runtime_10cores.rds") %>%
    mutate(cores="10 cores")
gpu <- readRDS("tests_runtime_gpu.rds") %>%
    mutate(cores="GPU")
speed_bench <- bind_rows(single, multi, gpu) %>%
    method_names()
ggplot(speed_bench, aes(npoints, color=method_name)) +
    geom_line(aes(y=time/60)) +
    facet_wrap(~cores, nrow=1) +
    scale_y_log10(name="time / min", labels=scales::label_math(format=log10)) +
    geom_line(aes(y=time), color="gray", linetype="dashed", data=data.frame(npoints=500:20000) %>% mutate(time=0.000002 * npoints^2 + 0.5)) +
    geom_line(aes(y=time), color="gray", linetype="twodash", data=data.frame(npoints=500:20000) %>% mutate(time=0.0000000004 * npoints^3 + 0.5)) +
    annotate(geom="text", x=7500, y=1000, label="N^3", parse=TRUE) +
    annotate(geom="segment", x=7500, xend=8000, yend=0.0000000004 * 8000^3 + 0.5 + 20, y=800, color="darkgray") +
    annotate(geom="text", x=10500, y=1500, label="N^2", parse=TRUE) +
    annotate(geom="segment", x=10500, y=1300, xend=11000, yend=0.000002 * 11000^2 + 0.5 + 20, color="darkgray") +
    labs(x="number of locations", color=NULL) +
    scale_x_continuous(limits=c(0, NA), expand=expansion()) +
    guides(color=guide_legend(override.aes=list(size=2)))
```
```{r figure 1e, fig.width=8, fig.height=2.96}
filter(speed_bench, method %in% c("omnibus", "spark", "spatialde1")) %>%
    mutate(method_name=recode(method_name, omnibus="SpatialDE2")) %>%
    ggplot(aes(npoints, color=method_name)) +
        geom_line(aes(y=time/60)) +
        facet_wrap(~cores, nrow=1) +
        scale_y_log10(name="time / min", labels=scales::label_math(format=log10)) +
        geom_line(aes(y=time), color="gray", linetype="dashed", data=data.frame(npoints=500:20000) %>% mutate(time=0.00000012 * npoints^2 + 0.1)) +
        geom_line(aes(y=time), color="gray", linetype="twodash", data=data.frame(npoints=500:20000) %>% mutate(time=0.0000000004 * npoints^3 + 0.5)) +
        annotate(geom="text", x=7500, y=1000, label="N^3", parse=TRUE) +
        annotate(geom="segment", x=7500, xend=8000, yend=0.0000000004 * 8000^3 + 0.5 + 20, y=800, color="darkgray") +
        annotate(geom="text", x=7500, y=30, label="N^2", parse=TRUE) +
        annotate(geom="segment", x=7500, y=20, xend=8000, yend=0.00000012 * 8000^2 + 0.1 + 1, color="darkgray") +
        labs(x="number of locations", color=NULL) +
        scale_x_continuous(limits=c(0, NA), expand=expansion()) +
        guides(color=guide_legend(override.aes=list(size=2))) +
        theme(legend.position="top")
```
```{r figure 1d, fig.width=8, fig.height=2.96}
clusterspeed <- read_feather("clustering_speed_benchmark.feather")
mutate(clusterspeed, method=recode(method, leiden="Leiden", spatialde2="SpatialDE2"), cores=recode(ncores, `1`="1 core", `10`="10 cores", gpu="GPU")) %>%
    ggplot(aes(npoints, time_seconds / 60, color=method)) +
        geom_line() +
        facet_wrap(~cores, nrow=1) +
        labs(x="number of locations", y="time / min", color=NULL) +
        scale_x_continuous(limits=c(0, NA), expand=expansion()) +
        guides(color=guide_legend(override.aes=list(size=2))) +
        theme(legend.position="top") 
```

