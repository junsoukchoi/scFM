scFM: Bayesian segmented Gaussian copula factor model for single-cell
sequencing data
================

- [Installation](#installation)
- [Example](#example)

The R package `scFM` implements an efficient data-augmented Gibbs
sampler for the scFM model, a novel Bayesian segmented Gaussian copula
factor model. This model is designed for analyzing single-cell
sequencing data, specifically addressing challenges like excessive low
counts and high skewness, while automatically determining the number of
latent factors.

## Installation

To install the latest version of the R package `scFM` from Github, use

``` r
library(devtools)
devtools::install_github("junsoukchoi/scFM")
```

## Example
