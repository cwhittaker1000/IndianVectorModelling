// Specifying all the includes and depends required to run the model
#include "Initial_State_Sample.hpp"
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]

// Linear
// [[Rcpp::export]]
Rcpp::NumericVector initial_state_sample(Rcpp::NumericVector fitted_parameters, Rcpp::NumericVector static_parameters) {

  double dE = fitted_parameters[0];
  double dL = fitted_parameters[1];
  double dP = fitted_parameters[2];
  double muE0 = fitted_parameters[3];
  double muL0 = fitted_parameters[4];
  double muP = fitted_parameters[5];
  double muM = fitted_parameters[6];
  double lambda = fitted_parameters[7];
  double tau = fitted_parameters[8];
  double beta = fitted_parameters[9];
  double dt = static_parameters[0];
  double z = 150; // z = initial E - L, but unclear how we calculate it
  double K = 150; // again need to work out how to calculate and incorporate this in

  double a = (0.5 * dL * dP)/(muM * (dP + muP));
  double muE = muE0 * exp(z/K);
  double muL = muL0 * exp(lambda * z/K);

  int E0 = ((beta * a) * z)/(dE + (beta * a) + muE);
  int L0 = (dE * z)/(dE + dL + muL);
  int M0 = (0.5 * dL * dP * L0)/muM*(dP+muP);
  int P0 = (2 * muM * M0) / dP;


  return Rcpp::NumericVector::create(E0, L0, M0, P0);

}
