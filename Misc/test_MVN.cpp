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
arma::mat testmvrnormArma(arma::mat mu, arma::mat sigma) {

  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(1, ncols);
  arma::vec eigenvalues_calc = arma::eig_sym(sigma);
  double minimum = eigenvalues_calc.min();

  if (minimum <= 0.0) {

    arma::vec eigenvalues;
    arma::mat eigenvectors;
    arma::eig_sym(eigenvalues, eigenvectors, sigma);
    //Rcpp::Rcout << "Eigenvalues are" << eigenvalues << std::endl;
    //Rcpp::Rcout << "Eigenvectors are" << eigenvectors << std::endl;

    arma::mat new_matrix(eigenvalues.size(), eigenvalues.size(), arma::fill::zeros);
    for (int i = 0; i < eigenvalues.size(); i++) {
      for (int j = 0; j < eigenvalues.size(); j++) {
        if(i == j) {
          if(eigenvalues[i] >= 0.0) {
            new_matrix(i, j) = eigenvalues[i];
          }
          else {
            new_matrix(i, j) = 0.001;
          }
        }
      }
    }
    //Rcpp::Rcout << "New Matrix is" << new_matrix << std::endl;


    // calculate diagonal scaling matrix
    arma::mat scaling_matrix(eigenvalues.size(), eigenvalues.size(), arma::fill::zeros);

    for (int i = 0; i < eigenvalues.size(); i++) {

      double res = 0;

      for (int j = 0; j < eigenvalues.size(); j++) {

        res = res + pow(eigenvectors(i, j), 2) * new_matrix(j, j);

      }

      scaling_matrix(i, i) = 1/res;

    }
    //Rcpp::Rcout << "Scaling Matrix is" << scaling_matrix << std::endl;


    // calculate the square root of the scaling matrix
    arma::mat scaling_sq_root = sqrt(scaling_matrix);
    //Rcpp::Rcout << "sqrtScaling Matrix is" << scaling_sq_root << std::endl;


    // calculate the squre root of the new_matrix
    arma::mat new_matrix_sq_root = sqrt(new_matrix);
    //Rcpp::Rcout << "sqrtNew Matrix is" << new_matrix_sq_root << std::endl;


    arma::mat decomposed_matrix = scaling_sq_root * eigenvectors * new_matrix_sq_root;
    arma::mat corrected_matrix = decomposed_matrix * arma::trans(decomposed_matrix);
    //Rcpp::Rcout << "decomposed_matrix is" << decomposed_matrix << std::endl;
    //Rcpp::Rcout << "corrected_matrix is" << corrected_matrix << std::endl;

    arma::mat variances(eigenvalues.size(), eigenvalues.size(), arma::fill::zeros);
    for (int k = 0; k < eigenvalues.size(); k++) {
      for (int l = 0; l < eigenvalues.size(); l++) {
        if (k == l) {
          variances(k, l) = sigma(k, l);
        }
      }
    }
    arma::mat stdev = sqrt(variances);

    arma::mat corrected_covariance = stdev * corrected_matrix * stdev;

    return(corrected_covariance);
  }

  // this next bit is a placeholder until I can get the truncated normal going.
  // prevents it from returning negative numbers
  else {

    arma::mat MVN_samples = mu + Y * arma::chol(sigma);
    for (int i = 0; i < MVN_samples.size(); i++) {
      if (MVN_samples(0, i) < 0) {
        MVN_samples(0, i) = 0;
      }
    }
    return(MVN_samples);
  }

}

