#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' @export
// [[Rcpp::export]]
arma::mat mvrnormArma(arma::mat mu, arma::mat sigma) {

  int ncols = sigma.n_cols;

  arma::mat Y = arma::randn(1, ncols);

  return mu.t() + Y * arma::chol(sigma);
}

