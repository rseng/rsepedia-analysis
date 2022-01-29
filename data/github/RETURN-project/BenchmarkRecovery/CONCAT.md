[![License](https://img.shields.io/github/license/RETURN-project/BenchmarkRecovery)](https://choosealicense.com/licenses/apache-2.0/)
[![Build Status](https://github.com/RETURN-project/BenchmarkRecovery/workflows/R-CMD-check/badge.svg?branch=master)](https://github.com/RETURN-project/BenchmarkRecovery/actions)
[![codecov](https://codecov.io/gh/RETURN-project/BenchmarkRecovery/graph/badge.svg)](https://codecov.io/gh/RETURN-project/BenchmarkRecovery)
[![codecov](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4320502.svg)](https://doi.org/10.5281/zenodo.4320502)
[![fair-software.eu](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B%20%20%E2%97%8F%20%20%E2%97%8B-orange)](https://fair-software.eu)

# Benchmarking recovery indicators derived from remote sensing time series

This project simulates Landsat data and evaluates the performance of recovery indicators with respect to data and disturbance characteristics.

## Background

The context of this project is the study of the recovery of tropical forests after an abrupt disturbance (typically a forest fire) using satellite images as a data source.

The speed of recovery after a disturbance is known to be correlated with the concept of resilience. This is true not only for forests, but for many dynamical systems. To put it simply: forests that recover fast are more resilient. Forests that recover slowly may be in danger of permanent disappearance.

The specialized literature proposes different metrics for measuring the recovery speed. The performance of these metrics depends on many factors. Some of them are natural, such as the intensity of the perturbation or the seasonality. Others are technical, such as the sampling frequency or the spatial resolution.

## Purpose


The purpose of this project is to **efficiently** **compare** the  **reliability** of different post-disturbance recovery **metrics**.

## Mechanics

1. **Infers** time series' **parameters** and characteristics from optical satellite image data​
2. Uses those parameters to **create** a large collection of synthetic (but realistic) **time series**
3. **Calculates** several state-of-the-art recovery **metrics**

## Simplified workflow

![Simplified workflow](./img/flow.png)

## Install

You can install the master version from R via:

```r
library(devtools)
install_github("RETURN-project/BenchmarkRecovery")
```

## Citation info

Please cite as:
```
Wanda De Keersmaecker, & Pablo Rodríguez-Sánchez. (2020, December 14). RETURN-project/BenchmarkRecovery: Benchmarking recovery metrics derived from remote sensing time series (Version v1.0). Zenodo. http://doi.org/10.5281/zenodo.432050
```
Citation can be exported to different formats (BibTeX, JSON, ...) with our [Zenodo link](https://zenodo.org/record/4320503#.X9dFkFOYWhc).

## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## License

[Apache](https://choosealicense.com/licenses/apache/)
