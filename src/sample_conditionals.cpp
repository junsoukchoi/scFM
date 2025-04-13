#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat sample_Z(const arma::mat& Z, const arma::mat& delta, const arma::mat& Lambda, const arma::mat& U, const arma::mat& Sigma, const arma::mat& invD, const arma::mat& X, double beta)
{
   int n = X.n_rows, p = X.n_cols;
   arma::mat Z_new(Z.begin(), n, p);
   Environment pkg1 = Environment::namespace_env("tmvnsim");
   Function f1 = pkg1["tmvnsim"];
   Environment pkg2 = Environment::namespace_env("TruncatedNormal");
   Function f2 = pkg2["rtmvnorm"];
   
   arma::mat Mean = invD * Lambda * U.t(); 
   arma::mat Var  = invD * Sigma * invD;
   for (int i = 0; i < n; i++)
   {
      double n_trunc = sum(X.row(i) <= beta);
      if (n_trunc > 0)
      {
         arma::uvec ind_trunc(n_trunc, arma::fill::none);
         arma::vec  lower(n_trunc, arma::fill::value(-arma::datum::inf));
         arma::vec  upper(n_trunc, arma::fill::value(arma::datum::inf));
         arma::uvec ind_trunc_k = arma::find(X.row(i) == 0);
         double  n_tunc_km1 = ind_trunc_k.n_elem, n_trunc_k = ind_trunc_k.n_elem;
         if (n_trunc_k > 0)
         {
            ind_trunc.subvec(0, n_trunc_k - 1) = ind_trunc_k;
            upper.subvec(0, n_trunc_k - 1)     = delta(ind_trunc_k, arma::linspace<arma::uvec>(0, 0, 1));
         }
         
         for (int k = 1; k < beta + 1; k++)
         {
            ind_trunc_k = arma::find(X.row(i) == k);
            n_trunc_k = n_tunc_km1 + ind_trunc_k.n_elem;
            if (n_trunc_k > n_tunc_km1)
            {
               ind_trunc.subvec(n_tunc_km1, n_trunc_k - 1) = ind_trunc_k;
               lower.subvec(n_tunc_km1, n_trunc_k - 1)     = delta(ind_trunc_k, arma::linspace<arma::uvec>(k - 1, k - 1, 1));
               upper.subvec(n_tunc_km1, n_trunc_k - 1)     = delta(ind_trunc_k, arma::linspace<arma::uvec>(k, k, 1));
               n_tunc_km1 = n_trunc_k;
            }
         }
         
         arma::vec cmean = Mean(ind_trunc, arma::linspace<arma::uvec>(i, i, 1));
         arma::mat cVar  = Var(ind_trunc, ind_trunc);
         
         NumericVector cmean_NumVec(cmean.begin(), cmean.end());
         NumericMatrix cVar_NumMat(cVar.n_rows, cVar.n_cols, cVar.begin());
         NumericVector lower_NumVec(lower.begin(), lower.end());
         NumericVector upper_NumVec(upper.begin(), upper.end());
         List          out = f1(1, n_trunc, lower_NumVec, upper_NumVec, Named("means") = cmean_NumVec, Named("sigma") = cVar_NumMat);
         NumericVector sample = out[0];
         if (any((sample < lower_NumVec) | (sample > upper_NumVec)).is_true())
         {
            sample = f2(1, cmean_NumVec, cVar_NumMat, lower_NumVec, upper_NumVec);
         }
         
         arma::rowvec sample_arma(sample.begin(), sample.size());
         Z_new(arma::linspace<arma::uvec>(i, i, 1), ind_trunc) = sample_arma;
      }
  }
   
   return Z_new;
}

// [[Rcpp::export]]
arma::mat sample_delta(const arma::mat& Z, const arma::mat& X, double beta)
{
   double p = Z.n_cols;
   arma::mat delta_new(p, beta + 1, arma::fill::none);
   for (int j = 0; j < p; j++)
   {
      arma::vec X_j = X.col(j);
      arma::vec Z_j = Z.col(j);
      for (int k = 0; k < beta; k++)
      {
         delta_new(j, k) = arma::randu(arma::distr_param(max(Z_j(arma::find(X_j == k))), min(Z_j(arma::find(X_j == k + 1)))));
      }
      delta_new(j, beta) = arma::randu(arma::distr_param(max(Z_j(arma::find(X_j == beta))), min(Z_j(arma::find(X_j > beta)))));
   }
   
   return delta_new;
}

// [[Rcpp::export]]
arma::mat sample_Lambda(const arma::mat& Z, const arma::mat& U, const arma::mat& Sigma, const arma::vec& Phi, double tau, const arma::mat& Psi)
{
   double p = Z.n_cols, k = U.n_cols;
   arma::mat Lambda_new(p, k, arma::fill::none);
   arma::mat invXi(k, k, arma::fill::none);
   arma::mat Lt_invXi(k, k, arma::fill::none);
   arma::vec invXi_eta(k, arma::fill::none);
   arma::vec lambda_j(k, arma::fill::none);
   for (int j = 0; j < p; j++)
   {
      invXi = U.t() * U / Sigma(j, j) + arma::diagmat(arma::pow(Psi.row(j).t() % arma::square(tau * Phi), -1));   
      invXi = (invXi + invXi.t()) / 2;
      invXi_eta = U.t() * Z.col(j) / Sigma(j, j);
      Lt_invXi  = arma::chol(invXi);
      lambda_j  = arma::solve(arma::trimatu(Lt_invXi), arma::randn(k), arma::solve_opts::fast);
      lambda_j += arma::solve(arma::trimatu(Lt_invXi), arma::solve(arma::trimatl(Lt_invXi.t()), invXi_eta, arma::solve_opts::fast), arma::solve_opts::fast);
      Lambda_new.row(j) = lambda_j.t();
   }
   
   return Lambda_new;
}

// [[Rcpp::export]]
arma::mat sample_U(const arma::mat& Z, const arma::mat& Lambda, const arma::mat& invSigma)
{
   double n = Z.n_rows, k = Lambda.n_cols;
   arma::mat invOmega = Lambda.t() * invSigma * Lambda + arma::diagmat(arma::vec(k, arma::fill::value(1)));
   arma::mat invOmega_Nu = Lambda.t() * invSigma * Z.t();
   arma::mat Lt_invOmega = arma::chol(invOmega);
   arma::mat U_new  = arma::solve(arma::trimatu(Lt_invOmega), arma::randn(k, n), arma::solve_opts::fast).t();
   U_new += arma::solve(arma::trimatu(Lt_invOmega), arma::solve(arma::trimatl(Lt_invOmega.t()), invOmega_Nu, arma::solve_opts::fast), arma::solve_opts::fast).t();   
   return U_new;
}

// [[Rcpp::export]]
arma::vec sample_Phi(const arma::mat& Lambda, double alpha)
{
   double p = Lambda.n_rows, k = Lambda.n_cols;
   arma::vec H(k, arma::fill::none);
   Environment pkg = Environment::namespace_env("GIGrvg");
   Function f = pkg["rgig"];
   for (int l = 0; l < k; l++)
   {
      NumericVector out = f(1, alpha - p, 2 * arma::accu(arma::abs(Lambda.col(l))), 1);
      H(l) = out[0];
   }
   arma::vec Phi_new = (1 / arma::accu(H)) * H;
   return Phi_new;
}

// [[Rcpp::export]]
arma::vec sample_c(const arma::mat& U, const arma::vec& pi, const arma::mat& mu)
{
   double n = U.n_rows, H = mu.n_cols;
   arma::vec c_new(n, arma::fill::none);
   arma::vec prob_c(n, arma::fill::none);
   for (int i = 0; i < n; i++)
   {
      prob_c  = pi % arma::exp(-0.5 * arma::sum(arma::square(mu.each_col() - U.row(i).t()), 0).t());
      prob_c /= sum(prob_c);
      NumericVector prob_c_NumVec(prob_c.begin(), prob_c.end());
      c_new(i) = Rcpp::sample(H, 1, false, prob_c_NumVec)[0];
   }
   
   return c_new;
}

// [[Rcpp::export]]
arma::mat sample_mu_pi(const arma::mat& U, const arma::vec& c, double phi, double H, double alpha)
{
   double n = U.n_rows, k = U.n_cols;
   arma::mat mu_new(k, H, arma::fill::zeros);
   arma::vec sizes(H, arma::fill::zeros);
   for (int h = 1; h < H; h++)
   {
      arma::uvec hind = arma::find(c == h + 1);
      sizes(h)  = hind.n_elem;
      mu_new.col(h) = arma::mvnrnd(arma::sum(U.rows(hind), 0).t() / (sizes(h) + 1 / phi) , arma::diagmat(arma::vec(k, arma::fill::value(1 / (sizes(h) + 1 / phi)))));
   }
   
   sizes(0) = n - sum(sizes.tail(H - 1));
   arma::rowvec Y(H, arma::fill::none);
   for (int h = 0; h < H; h++)
   {
      Y(h) = Rcpp::rgamma(1, alpha + sizes(h), 1)[0];
   }
   arma::rowvec pi_new = Y / sum(Y);
   return arma::join_cols(mu_new, pi_new);
}
