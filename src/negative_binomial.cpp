// Specifying all the includes and depends required to run the model
#include "Negative_Binomial.hpp"
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
double Negative_Binomial(double k, double n, double r, double p) {
  double m = p * n;
  double res = (lgamma(r + k) - (lgamma(k + 1) + lgamma(r))) + k * (log(m) - (log(r + m))) + r * (log(r) - log(r + m));
  return(res);
}
