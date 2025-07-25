scFM: Bayesian segmented Gaussian copula factor model for single-cell
sequencing data
================

- [Installation](#installation)
- [Example](#example)

The R package `scFM` implements an efficient data-augmented Gibbs
sampler for scFM (Choi et al., 2025+), a novel Bayesian segmented
Gaussian copula factor model for single-cell sequencing data. This
method enables factor analysis of single-cell sequencing data by
effectively addressing challenges such as inflated low counts (dropout)
and high skewness in the data, while also automatically determining the
number of latent factors.

**Reference:**

Choi, J., Chung, H. C., Gaynanova, I., & Ni, Y. (2025+). Bayesian
segmented Gaussian copula factor model for single-cell sequencing data.

## Installation

To install the latest version of the R package `scFM` from Github, use

``` r
library(devtools)
devtools::install_github("junsoukchoi/scFM")
```

## Example

The following example demonstrates how to apply our proposed scFM to the
scRNA-seq dataset analyzed in Section 5 of the paper.

``` r
# load the scFM package along with necessary libraries
library(scFM)
library(GIGrvg)
library(statmod)
library(ggplot2)

# set a random seed
set.seed(7)

# load the built-in LCL dataset, corresponding to the scRNA-seq data analyzed in Section 5
data(LCL)
x = LCL

#################### set hyperparameters ####################
m = 1
k_max = 8
priors = list()
priors$Sigma = c(0.1, 0.1)
priors$phi   = 0.5

#################### initialize parameters ####################
# initialize Z using the scaled empirical CDF of each variable
starting = list()
n = nrow(x)
p = ncol(x)
ecdf.scale = n / (n + 1)
eFx   = apply(x, 2, ecdf)
eFx   = lapply(eFx, function(x){  function(y) ecdf.scale * x(y)  })
eFxx  = Map(function(f, x) do.call(f, list(x)), eFx, plyr::alply(x, 2))
starting$Z     = matrix(unlist(lapply(eFxx, function(pr) qnorm(pr))), n, p)

# initialize delta based on the proportion of low integer values (0 and 1, assuming m = 1)
starting$delta = matrix(0, p, 2)
starting$delta[ , 1] = qnorm(colMeans(x == 0))
starting$delta[ , 2] = qnorm(colMeans(x <= 1))

# initialize Lambda and U using typical factor analysis
out0 = factanal(starting$Z, factors = k_max, scores = "regression")
starting$Lambda = as.matrix(out0$loadings)
starting$U = out0$scores

# initialize the remaining parameters by sampling from their full conditional distributions, 
# given the initialized values of Z, delta, Lambda, and U.
Sigma = diag(1, p)
phi = rep(1 / k_max, k_max)
tau = 1
Xi = matrix(rexp(p * k_max, rate = 0.5), p, k_max)
for (iter in 1 : 1000) 
{
    Sigma[diag(rep(TRUE, p))] = 1 / rgamma(p, shape = priors$Sigma[1] + rep(n / 2, p), 
                                           rate = priors$Sigma[2] + 
                                             0.5 * colSums((starting$Z - tcrossprod(starting$U, starting$Lambda))^2))
  
    H = rep(NA, k_max)
    for (l in 1 : k_max) 
    {
        H[l] = rgig(1, lambda = priors$phi - p, chi = 2 * sum(abs(starting$Lambda[ , l])), psi = 1)
    }
  
    phi = H / sum(H)
  
    tau = rgig(1, lambda = k_max * (priors$phi - p), chi = 2 * sum(t(abs(starting$Lambda)) / phi), psi = 1)
  
    Xi_tilde = rinvgauss(p * k_max, mean = c(t(phi * t(tau / abs(starting$Lambda)))))
    Xi = matrix(1 / Xi_tilde, p, k_max)
}

starting$Sigma = Sigma
starting$phi   = phi
starting$tau   = tau
starting$Xi    = Xi

#################### run our proposed data-augmented Gibbs sampler for scFM ####################
nmcmc   = 10000
nburnin = 5000
p.time = proc.time()
out    = gibbs_scFM(x, k_max, m, starting, priors, nmcmc = nmcmc, nburnin = nburnin)
proc.time() - p.time

#################### conduct posterior inference ####################
# compute the posterior mean of the factor loadings and scores
loading_est = apply(out$Lambda, c(1, 2), mean) 
score_est = apply(out$U, c(1, 2), mean) 

# estimate the number of signicant factors
sig_facto = matrix(0, dim(out$Lambda)[2], dim(out$Lambda)[3])
for (i in 1 : dim(out$Lambda)[3])
{
  kmeans_lambda = kmeans(t(out$Lambda[ , , i]), centers = 2, nstart = 10)
  sig_facto[which(kmeans_lambda$cluster == which.max(rowSums(kmeans_lambda$centers^2))), i] = 1
}

k_hat = as.integer(names(which.max(table(colSums(sig_facto)))))
k_hat

# the posterior inclusion probability-based method also yields k_hat = 2
#calculate_eFDR = function(c, pip)
#{
#  return(sum((pip > c) * (1 - pip)) / (sum(pip > c) + 1e-8))
#}
#
#cutoffs = seq(0.005, 0.995, by = 0.005)
#eFDR = sapply(cutoffs, function(c) calculate_eFDR(c, rowMeans(sig_facto)))
#k_hat = sum(rowMeans(sig_facto) > cutoffs[which(eFDR < 0.01)[1]])
#k_hat

# given k_hat, identify significant factors 
ids_sig_factors = order(sqrt(colSums(loading_est^2)), decreasing = TRUE)[1:k_hat]
score_est_sig   = score_est[ , ids_sig_factors]

# apply the k-means clustering with k = 3 (chosen based on the elbow method);
# this reveals three distinct sell subpopulations
out_cell_types = kmeans(score_est_sig, centers = 3, nstart = 100, iter.max = 10000)
df = data.frame(facto1 = score_est_sig[ , 1], facto2 = score_est_sig[ , 2], Cluster = factor(out_cell_types$cluster))
plot_cell_types = ggplot(df[-4148, ], aes(facto1, facto2, color = Cluster)) + 
  geom_point(alpha = 0.3) + labs(x = "Factor 1", y = "Factor 2")
plot_cell_types

# plot the loadings of each signifcant factor, sorted by their absolute magnitude;
# this provides interpretation for the factors
rownames(loading_est) = colnames(x)
order1 = order(abs(loading_est[ , ids_sig_factors[1]]), decreasing = TRUE)
df = data.frame(gene = factor(colnames(x), levels = colnames(x)[order1]), value = abs(loading_est[ , 1]))
loading1 = ggplot(df, aes(x = gene, y = value)) + geom_bar(stat = "identity") + 
  scale_x_discrete(guide = guide_axis(n.dodge = 5)) + labs(x = "Genes", y = "Absolute loadings on factor 1")
loading1
order2 = order(abs(loading_est[ , ids_sig_factors[2]]), decreasing = TRUE)
df = data.frame(gene = factor(colnames(x), levels = colnames(x)[order2]), value = abs(loading_est[ , 2]))
loading2 = ggplot(df, aes(x = gene, y = value)) + geom_bar(stat = "identity") + 
  scale_x_discrete(guide = guide_axis(n.dodge = 5)) + labs(x = "Genes", y = "Absolute loadings on factor 2")
loading2
```
