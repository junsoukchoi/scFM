#' Gibbs sampler for Bayesian segmented Gaussian copula factor model
#'
#' @param x Data matrix (n × p), typically with count-type entries.
#' @param k Number of latent factors.
#' @param m To be written...
#' @param starting A list of starting values for the parameters (e.g., Z, delta, Lambda, U, etc.).
#' @param priors A list of prior hyperparameters (e.g., for Sigma and Phi).
#' @param nmcmc Number of MCMC iterations.
#' @param nburnin Burn-in period.
#'
#' @return A list containing MCMC samples for Lambda and U.
#' @export
#'
#' @examples
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
