#' Data-augmented Gibbs sampler for Bayesian segmented Gaussian copula factor model (scFM)
#'
#' @param x An n-by-p matrix where each row represents the read counts for a single cell in single-cell sequencing data.
#' @param k An integer indicating the maximum allowable number of latent factors.
#' @param m An integer specifying a model hyperparameter in scFM that serves as a threshold to distinguish between inflated low counts and typical higher counts in single-cell sequencing data.
#' @param starting A list of initial values for the data-augmented Gibbs sampler for scFM. The list should include the tags \code{Z}, \code{delta}, \code{Lambda}, \code{U}, \code{Sigma}, \code{Phi}, \code{tau}, and \code{Psi}, each representing the initial MCMC value for the corresponding scFM parameter.
#' @param priors A list of hyperparameter values. Valid tags are \code{Sigma} and \code{Phi}, specifying the hyperparameters for the inverse-gamma prior on error variances (i.e., \eqn{a_{\sigma}}, \eqn{b_{\sigma}}) and the Dirichlet-Laplace prior on factor loadings (i.e., \eqn{\alpha}), respectively.
#' @param nmcmc An integer specifying the total number of MCMC iterations.
#' @param nburnin An integer indicating the number of burn-in samples.
#'
#' @return A list containing MCMC samples of the factor loadings and factor scores: \code{Lambda} holds the MCMC samples for factor loadings and \code{U} holds the MCMC samples for factor scores.
#' @export
#'
#' @examples 
#' # load the scFM package along with necessary libraries
#' library(scFM)
#' library(GIGrvg)
#' library(statmod)
#' library(ggplot2)
#' 
#' # set a random seed
#' set.seed(7)
#' 
#' # load the built-in LCL dataset, corresponding to the scRNA-seq data analyzed in Section 5
#' data(LCL)
#' x = LCL
#' 
#' #################### set hyperparameters ####################
#' m = 1
#' k_max = 8
#' priors = list()
#' priors$Sigma = c(0.1, 0.1)
#' priors$Phi   = 0.5
#' 
#' #################### initialize parameters ####################
#' # initialize Z using the scaled empirical CDF of each variable
#' starting = list()
#' n = nrow(x)
#' p = ncol(x)
#' ecdf.scale = n / (n + 1)
#' eFx   = apply(x, 2, ecdf)
#' eFx   = lapply(eFx, function(x){  function(y) ecdf.scale * x(y)  })
#' eFxx  = Map(function(f, x) do.call(f, list(x)), eFx, plyr::alply(x, 2))
#' starting$Z     = matrix(unlist(lapply(eFxx, function(pr) qnorm(pr))), n, p)
#' 
#' # initialize delta based on the proportion of low integer values (0 and 1, assuming m = 1)
#' starting$delta = matrix(0, p, 2)
#' starting$delta[ , 1] = qnorm(colMeans(x == 0))
#' starting$delta[ , 2] = qnorm(colMeans(x <= 1))
#' 
#' # initialize Lambda and U using typical factor analysis
#' out0 = factanal(starting$Z, factors = k_max, scores = "regression")
#' starting$Lambda = as.matrix(out0$loadings)
#' starting$U = out0$scores
#' 
#' # initialize the remaining parameters by sampling from their full conditional distributions, 
#' # given the initialized values of Z, delta, Lambda, and U.
#' Sigma = diag(1, p)
#' Phi = rep(1 / k_max, k_max)
#' tau = 1
#' Psi = matrix(rexp(p * k_max, rate = 0.5), p, k_max)
#' 
#' for (iter in 1 : 1000) 
#' {
#'    Sigma[diag(rep(TRUE, p))] = 1 / rgamma(p, shape = priors$Sigma[1] + rep(n / 2, p), 
#'                                           rate = priors$Sigma[2] + 
#'                                              0.5 * colSums((starting$Z - tcrossprod(starting$U, starting$Lambda))^2))
#'    
#'    H = rep(NA, k_max)
#'    for (l in 1 : k_max) 
#'    {
#'       H[l] = rgig(1, lambda = priors$Phi - p, chi = 2 * sum(abs(starting$Lambda[ , l])), psi = 1)
#'    }
#'    
#'    Phi = H / sum(H)
#'    
#'    tau = rgig(1, lambda = k_max * (priors$Phi - p), chi = 2 * sum(t(abs(starting$Lambda)) / Phi), psi = 1)
#'    
#'    Psi_tilde = rinvgauss(p * k_max, mean = c(t(Phi * t(tau / abs(starting$Lambda)))))
#'    Psi = matrix(1 / Psi_tilde, p, k_max)
#' }
#' 
#' starting$Sigma = Sigma
#' starting$Phi   = Phi
#' starting$tau   = tau
#' starting$Psi   = Psi
#' 
#' #################### run our proposed data-augmented Gibbs sampler for scFM ####################
#' nmcmc   = 10000
#' nburnin = 5000
#' p.time = proc.time()
#' out    = gibbs_scFM(x, k_max, m, starting, priors, nmcmc = nmcmc, nburnin = nburnin)
#' proc.time() - p.time
#' 
#' #################### conduct posterior inference ####################
#' # compute the posterior mean of the factor loadings and scores
#' loading_est = apply(out$Lambda, c(1, 2), mean) 
#' score_est = apply(out$U, c(1, 2), mean) 
#' 
#' # estimate the number of signicant factors
#' sig_facto = matrix(0, dim(out$Lambda)[2], dim(out$Lambda)[3])
#' for (i in 1 : dim(out$Lambda)[3])
#' {
#'    kmeans_lambda = kmeans(t(out$Lambda[ , , i]), centers = 2, nstart = 10)
#'    sig_facto[which(kmeans_lambda$cluster == which.max(rowSums(kmeans_lambda$centers^2))), i] = 1
#' }
#' 
#' k_hat = as.integer(names(which.max(table(colSums(sig_facto)))))
#' k_hat
#' 
#' # the posterior inclusion probability-based method also yields k_hat = 2
#' #calculate_eFDR = function(c, pip)
#' #{
#' #  return(sum((pip > c) * (1 - pip)) / (sum(pip > c) + 1e-8))
#' #}
#' #
#' #cutoffs = seq(0.005, 0.995, by = 0.005)
#' #eFDR = sapply(cutoffs, function(c) calculate_eFDR(c, rowMeans(sig_facto)))
#' #k_hat = sum(rowMeans(sig_facto) > cutoffs[which(eFDR < 0.01)[1]])
#' #k_hat
#' 
#' # given k_hat, identify significant factors 
#' ids_sig_factors = order(sqrt(colSums(loading_est^2)), decreasing = TRUE)[1:k_hat]
#' score_est_sig   = score_est[ , ids_sig_factors]
#' 
#' # apply the k-means clustering with k = 3 (chosen based on the elbow method);
#' # this reveals three distinct sell subpopulations
#' out_cell_types = kmeans(score_est_sig, centers = 3, nstart = 100, iter.max = 10000)
#' df = data.frame(facto1 = score_est_sig[ , 1], facto2 = score_est_sig[ , 2], Cluster = factor(out_cell_types$cluster))
#' plot_cell_types = ggplot(df[-4148, ], aes(facto1, facto2, color = Cluster)) + 
#'    geom_point(alpha = 0.3) + labs(x = "Factor 1", y = "Factor 2")
#' plot_cell_types
#' 
#' # plot the loadings of each signifcant factor, sorted by their absolute magnitude;
#' # this provides interpretation for the factors
#' rownames(loading_est) = colnames(x)
#' order1 = order(abs(loading_est[ , ids_sig_factors[1]]), decreasing = TRUE)
#' df = data.frame(gene = factor(colnames(x), levels = colnames(x)[order1]), value = abs(loading_est[ , 1]))
#' loading1 = ggplot(df, aes(x = gene, y = value)) + geom_bar(stat = "identity") + 
#'    scale_x_discrete(guide = guide_axis(n.dodge = 5)) + labs(x = "Genes", y = "Absolute loadings on factor 1")
#' loading1
#' order2 = order(abs(loading_est[ , ids_sig_factors[2]]), decreasing = TRUE)
#' df = data.frame(gene = factor(colnames(x), levels = colnames(x)[order2]), value = abs(loading_est[ , 2]))
#' loading2 = ggplot(df, aes(x = gene, y = value)) + geom_bar(stat = "identity") + 
#'    scale_x_discrete(guide = guide_axis(n.dodge = 5)) + labs(x = "Genes", y = "Absolute loadings on factor 2")
#' loading2
gibbs_scFM = function(x, k, m, starting, priors, nmcmc, nburnin)
{
   n = nrow(x)
   p = ncol(x)
   
   # initialize parameters
   Z      = starting$Z
   delta  = starting$delta
   Lambda = starting$Lambda
   U      = starting$U
   Sigma  = starting$Sigma
   Phi    = starting$Phi
   tau    = starting$tau
   Psi    = starting$Psi
   
   # initialize MCMC samples
   #Z_gibbs      = array(0, dim = c(n, p, nmcmc - nburnin))
   #delta_gibbs  = array(0, dim = c(p, 2, nmcmc - nburnin))
   Lambda_gibbs = array(0, dim = c(p, k, nmcmc - nburnin))
   U_gibbs      = array(0, dim = c(n, k, nmcmc - nburnin))
   #Sigma_gibbs  = array(0, dim = c(p, p, nmcmc - nburnin))
   #Phi_gibbs    = matrix(0, k, nmcmc - nburnin)
   #tau_gibbs    = rep(0, nmcmc - nburnin)
   #Psi_gibbs    = array(0, dim = c(p, k, nmcmc - nburnin))
   
   # iterate
   for (iter in 1 : nmcmc)
   {
      # sample Z from its full conditional
      D = sqrt(diag(tcrossprod(Lambda) + Sigma))
      #for (i in 1 : n)
      #{
      #   ind0 = which(x[i, ] == 0)
      #   ind1 = which(x[i, ] == 1)
      #   tind = c(ind0, ind1)
      #   
      #   cmu    = (tcrossprod(Lambda, U[i, , drop = FALSE]) / D)[tind, ]
      #   cSigma = (Sigma / D^2)[tind, tind]
      #   
      #   lower = rep(-Inf, length(tind))
      #   if (length(ind1) > 0) lower[length(ind0) + (1 : length(ind1))] = delta[ind1, 1]
      #   upper = rep(Inf, length(tind))
      #   if (length(ind0) > 0) upper[1 : length(ind0)] = delta[ind0, 1]
      #   if (length(ind1) > 0) upper[length(ind0) + (1 : length(ind1))] = delta[ind1, 2]
      #   
      #   #zt_new = tmvtnorm::rtmvnorm(1, mean = cmu, sigma = cSigma, lower = lower, upper = upper, 
      #   #                            algorithm = "gibbs", burn.in.samples = 50, thinning = 1)
      #   #Z[i, tind] = zt_new
      #   Z[i, tind] = tmvnsim::tmvnsim(1, length(tind), lower = lower, upper = upper, means = cmu, sigma = cSigma)$samp
      #   if (any(Z[i, tind] < lower | Z[i, tind] > upper))
      #      Z[i, tind] = TruncatedNormal::rtmvnorm(1, mu = cmu, sigma = cSigma, lb = lower, ub = upper)
      #}
      Z = sample_Z(Z, delta, Lambda, U, Sigma, diag(1 / D), x, m)
      
      # sample delta from its full conditional
      #z0_max = rep(NA, p)
      #z1_min = rep(NA, p)
      #z1_max = rep(NA, p)
      #zo_min = rep(NA, p)
      #for (j in 1 : p)
      #{
      #   z0_max[j] = max(Z[x[ , j] == 0, j])
      #   z1_min[j] = min(Z[x[ , j] == 1, j])
      #   z1_max[j] = max(Z[x[ , j] == 1, j])
      #   zo_min[j] = min(Z[x[ , j] != 0 & x[ , j] != 1, j])
      #}
      #delta[ , 1] = runif(p, min = z0_max, max = z1_min)
      #delta[ , 2] = runif(p, min = z1_max, max = zo_min)
      delta = sample_delta(delta, Z, x, m)    # Rcpp function to sample delta
      
      # sample Lambda from its full conditional
      #for (j in 1 : p)
      #{
      #   invXi = crossprod(U) / Sigma[j, j] + diag(1 / (Psi[j, ] * tau^2 * Phi^2), k)
      #   invXi = (invXi + t(invXi)) / 2   # guarantee symmetry
      #   invXi = as(invXi, "sparseMatrix")
      #   chol_invXi = Matrix::Cholesky(invXi, LDL = FALSE)
      #   e_Xi = as.vector(Matrix::solve(chol_invXi, rnorm(k), system = "Lt"))
      #   invXi_eta = crossprod(U, Z[ , j]) / Sigma[j, j]
      #   eta = as.vector(Matrix::solve(chol_invXi, Matrix::solve(chol_invXi, invXi_eta, system = "L"), system = "Lt"))
      #   Lambda[j, ] = eta + e_Xi
      #}
      Lambda = sample_Lambda(Z, U, Sigma, Phi, tau, Psi)
      
      # sample U from its full conditional
      #invOmega_j    = crossprod(Lambda, Lambda / diag(Sigma)) + diag(rep(1, k))
      #invOmega_j    = (invOmega_j + t(invOmega_j)) / 2   # guarantee symmetry
      #invOmega      = do.call(Matrix::bdiag, rep(list(invOmega_j), n))
      #chol_invOmega = Matrix::Cholesky(invOmega, LDL = FALSE)
      #e_Omega       = as.vector(Matrix::solve(chol_invOmega, rnorm(n * k), system = "Lt"))
      #invOmega_nu   = c(crossprod(Lambda, t(Z) / diag(Sigma)))
      #nu            = as.vector(Matrix::solve(chol_invOmega, Matrix::solve(chol_invOmega, invOmega_nu, system = "L"), system = "Lt"))
      #U             = matrix(nu + e_Omega, n, k, byrow = TRUE)
      U = sample_U(Z, Lambda, diag(1 / diag(Sigma)))
      
      # sample Sigma from its full conditional
      Sigma[diag(rep(TRUE, p))] = 1 / rgamma(p, shape = priors$Sigma[1] + rep(n / 2, p), 
                                             rate = priors$Sigma[2] + 0.5 * colSums((Z - tcrossprod(U, Lambda))^2))
      
      # sample Phi from its full conditional
      H = rep(NA, k)
      for (l in 1 : k)
      {
         H[l] = rgig(1, lambda = priors$Phi - p, chi = 2 * sum(abs(Lambda[ , l])), psi = 1)
      }
      
      Phi = H / sum(H)
      
      # sample tau from its full conditional
      tau = rgig(1, lambda = k * (priors$Phi - p), chi = 2 * sum(t(abs(Lambda)) / Phi), psi = 1)
      
      # sample Psi from its full conditional
      Psi_tilde = rinvgauss(p * k, mean = c(t(Phi * t(tau / abs(Lambda)))))
      Psi = matrix(1 / Psi_tilde, p, k)
      
      # print progress
      if (iter %% 100 == 0) cat("iter =", iter, "\n")
      
      # store MCMC samples
      if (iter > nburnin)
      {
         #Z_gibbs[ , , iter - nburnin]      = Z
         #delta_gibbs[ , , iter - nburnin]  = delta
         Lambda_gibbs[ , , iter - nburnin] = Lambda
         U_gibbs[ , , iter - nburnin]      = U
         #Sigma_gibbs[ , , iter - nburnin]  = Sigma
         #Phi_gibbs[ , iter - nburnin]      = Phi
         #tau_gibbs[iter - nburnin]         = tau
         #Psi_gibbs[ , , iter - nburnin]    = Psi
      }
   }
   
   # return MCMC samples
   out = list()
   #out$Z      = Z_gibbs
   #out$delta  = delta_gibbs
   out$Lambda = Lambda_gibbs
   out$U      = U_gibbs
   #out$Sigma  = Sigma_gibbs
   #out$Phi    = Phi_gibbs
   #out$tau    = tau_gibbs
   #out$Psi    = Psi_gibbs
   return(out)
}
