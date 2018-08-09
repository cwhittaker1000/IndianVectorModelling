#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Multivariate Normal Sampler
//'
//' Random number generation for a multivariate normal (Gaussian) distribution. Uses a Cholesky
//' decomposition to propose new values from a multivariate normal distribution with a given mean
//' vector (mu) and a specified covariance matrix (sigma).
//'
//' Note: Currently (for my own modelling purposes) won't return values < 0.
//'
//' @param mu Mean of the multivariate normal distribution (has to be specified as a matrix with
//' 1 row and number of columns equal to the number of variables in the multivariate distribution e.g. 2 for bivariate)
//' @param sigma A covariance matrix associated with mu
//' @export
//'
// [[Rcpp::export]]
arma::mat mvrnormArma(arma::mat mu, arma::mat sigma) {

  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(1, ncols);
  arma::mat MVN_samples = mu + Y * arma::chol(sigma);

  // this next bit is a placeholder until I can get the truncated normal going.
  // prevents it from returning negative numbers

  for (int i = 0; i < MVN_samples.size(); i++) {
    if (MVN_samples(0, i) < 0) {
      MVN_samples(0, i) = 0;
    }
  }

  return(MVN_samples);
}

