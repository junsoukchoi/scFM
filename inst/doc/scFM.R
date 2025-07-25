## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  message = FALSE,
  warning = FALSE,
  collapse = TRUE,
  comment = NA,
  eval = TRUE,
  cache = TRUE,
  tidy = TRUE,
  tidy.opts = list(width.cutoff = 60),
  fig.width = 7,
  fig.height = 5.25,
  fig.align = "center"
)

## ----load-data, cache=TRUE----------------------------------------------------
# load necessary packages including scFM
library(scFM)
library(GIGrvg)
library(statmod)
library(ggplot2)
library(tibble)
library(SingleR)
library(celldex)

# set random seed
set.seed(7)

# load the built-in LCL dataset
data(LCL)
x = LCL
dim(x)      # check dimensions of the dataset

## ----set-hyperparameters, cache=TRUE------------------------------------------
# set k_{max} and m
k = 8
m = 1

# set hyperparameters for priors
priors = list()
priors$Sigma = c(0.1, 0.1)
priors$phi   = 0.5

## ----init-Z-delta, cache=TRUE-------------------------------------------------
# initialize parameters Z and delta
starting = list()
n = nrow(x)
p = ncol(x)
ecdf.scale = n / (n + 1)
eFx   = apply(x, 2, ecdf)
eFx   = lapply(eFx, function(x){  function(y) ecdf.scale * x(y)  })
eFxx  = Map(function(f, x) do.call(f, list(x)), eFx, plyr::alply(x, 2))
starting$Z     = matrix(unlist(lapply(eFxx, function(pr) qnorm(pr))), n, p)
starting$delta = matrix(0, p,  m + 1)
for (d in 0 : m)
{
   starting$delta[ , d + 1] = qnorm(colMeans(x <= d))
}

## ----init-Lambda-U, cache=TRUE------------------------------------------------
# initialize Lambda and U using standard factor analysis
out0 = factanal(starting$Z, factors = k, scores = "regression")
starting$Lambda = as.matrix(out0$loadings)
starting$U = out0$scores

## ----init-remainders, cache=TRUE----------------------------------------------
# initialize Sigma, phi, tau, and Xi via a short Gibbs sampling
Sigma = diag(1, p)
phi = rep(1 / k, k)
tau = 1
Xi = matrix(rexp(p * k, rate = 0.5), p, k)

for (iter in 1 : 1000) {
   Sigma[diag(rep(TRUE, p))] = 1 / rgamma(p, shape = priors$Sigma[1] + rep(n / 2, p), 
                                          rate = priors$Sigma[2] + 0.5 * colSums((starting$Z - tcrossprod(starting$U, starting$Lambda))^2))
   
   H = rep(NA, k)
   for (l in 1 : k) {
      H[l] = rgig(1, lambda = priors$phi - p, chi = 2 * sum(abs(starting$Lambda[ , l])), psi = 1)
   }
   phi = H / sum(H)
   tau = rgig(1, lambda = k * (priors$phi - p), chi = 2 * sum(t(abs(starting$Lambda)) / phi), psi = 1)
   Xi_tilde = rinvgauss(p * k, mean = c(t(phi * t(tau / abs(starting$Lambda)))))
   Xi = matrix(1 / Xi_tilde, p, k)
}

starting$Sigma = Sigma
starting$phi   = phi
starting$tau   = tau
starting$Xi    = Xi

## ----fit-scFM, cache=TRUE-----------------------------------------------------
# run the data-augmented Gibbs sampling for scFM
nmcmc   = 10000
nburnin = 5000
out    = gibbs_scFM(x, k, m = 1, starting, priors, nmcmc = nmcmc, nburnin = nburnin)

# compute posterior means of the factor loadings and scores
loading_est = apply(out$Lambda, c(1, 2), mean) 
score_est = apply(out$U, c(1, 2), mean) 

## ----identify-significant-factors, cache=TRUE---------------------------------
# infer the number of significant latent factors
sig_facto = matrix(0, dim(out$Lambda)[2], dim(out$Lambda)[3])
for (i in 1 : dim(out$Lambda)[3])
{
   kmeans_lambda = kmeans(t(out$Lambda[ , , i]), centers = 2, nstart = 10)
   sig_facto[which(kmeans_lambda$cluster == which.max(rowSums(kmeans_lambda$centers^2))), i] = 1
}

table(colSums(sig_facto))         # indicates the existence of two significant latent factors

# select the significant latent factors
sqrt(colSums(loading_est^2))      # the first two factors show the largest L2-norm
score_est = score_est[ , 1 : 2]

## ----interpret-factors, cache=TRUE--------------------------------------------
# display absolute values of the estimated factor loadings on the significant factors
rownames(loading_est) = colnames(x)

order1 = order(abs(loading_est[ , 1]), decreasing = TRUE)
df = data.frame(gene = factor(colnames(x), levels = colnames(x)[order1]), value = abs(loading_est[ , 1]))
loading1 = ggplot(df, aes(x = gene, y = value)) + geom_bar(stat = "identity") + scale_x_discrete(guide = guide_axis(n.dodge = 8)) + labs(x = "Genes", y = "Absolute loadings on Factor 1")
loading1

order2 = order(abs(loading_est[ , 2]), decreasing = TRUE)
df = data.frame(gene = factor(colnames(x), levels = colnames(x)[order2]), value = abs(loading_est[ , 2]))
loading2 = ggplot(df, aes(x = gene, y = value)) + geom_bar(stat = "identity") + scale_x_discrete(guide = guide_axis(n.dodge = 8)) + labs(x = "Genes", y = "Absolute loadings on Factor 2")
loading2

# check genes with the 10 largest abosolute loadings for each significant factor
sort(loading_est[rank(abs(loading_est[ , 1])) > 90, 1], decreasing = TRUE)
sort(loading_est[rank(abs(loading_est[ , 2])) > 90, 2], decreasing = TRUE)

## ----clustering, cache=TRUE---------------------------------------------------
# determine the number of clusters via the elbow method
n_clusters = 10
wss = numeric(n_clusters)
for (i in 1 : n_clusters)
{
   # fit the model for the given number of clusters
   km.out = kmeans(score_est, centers = i, nstart = 100, iter.max = 10000)
   # save the within cluster sum of squares
   wss[i] = km.out$tot.withinss
}

# produce a scree plot to determine the number of the underlying cell clusters
wss_df = tibble(clusters = 1 : n_clusters, wss = wss)
scree_plot = ggplot(wss_df, aes(x = clusters, y = wss, group = 1)) + 
   geom_point(size = 4) + 
   geom_line() +
   scale_x_continuous(breaks = c(2, 4, 6, 8, 10)) +
   xlab('Number of clusters')
scree_plot

# apply k-means clustering with the chosen k = 3
out_cell_types = kmeans(score_est, centers = 3, nstart = 100, iter.max = 10000)

## ----visualize-clustering, cache=TRUE-----------------------------------------
# visualize cell clustering (excluding an outlier)
df = data.frame(facto1 = score_est[ , 1], facto2 = score_est[ , 2], Cluster = factor(out_cell_types$cluster))
plot_cell_types = ggplot(df[-4148, ], aes(facto1, facto2, color = Cluster)) + geom_point(alpha = 0.3) + labs(x = "Factor 1", y = "Factor 2")
plot_cell_types

## ----singleR, cache=TRUE------------------------------------------------------
# annotate clusters using SingleR
hpca.se = HumanPrimaryCellAtlasData()
hpca.se
annotate_clusters = SingleR(test = t(x), ref = hpca.se, clusters = out_cell_types$cluster, labels = hpca.se$label.main)
annotate_clusters

