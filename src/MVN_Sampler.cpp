#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' @export
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

