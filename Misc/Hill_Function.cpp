#include <RcppArmadillo.h>
#include <vector>
#include <string>
#include <iostream>
// [[Rcpp::depends(RcppArmadillo)]]

//' @export
// [[Rcpp::export]]
double Hill_Function(double rainfall, double K_max, double a, double b) {
  double output = K_max / (1 + pow((rainfall / a), b));
  return(output);
}
