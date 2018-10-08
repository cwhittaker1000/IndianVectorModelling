// Specifying all the includes and depends required to run the model
#include "Initial_State_Sample.hpp"
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]

// Linear
//' @export
// [[Rcpp::export]]
std::vector <int> initial_state_sample(Rcpp::NumericVector fitted_parameters, Rcpp::NumericVector static_parameters,
                                       double initial_K) {

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
  double K = initial_K;
  double dt = static_parameters[0];

  double x = (lambda * (muL0 / muE0) ) -  ((1 / dE) / (1/ dL)) - ((lambda - 1) * (muL0 * (1/ dE)));
  double y = (lambda * (beta * muL0 * (1/ dE))) / ((2 * muE0 * muM * (1/ dL))*(1 + (1/ dP)*muP));
  double W = (-0.5 * x) + sqrt(0.25 * pow(x, 2) + y);
  int M = initial_K * ((W /(muL0 * (1/dE))) - (1 / (muL0 * (1/dL))) - 1);

  int M0 = M;
  int E0 = ((2 * W * muM * (1/ dL)) * (1 + ((1/ dP) * muP))) * M0;
  int L0 = ((2 * muM * (1/ dL)) * (1 * ((1/ dP) * muP))) * M0;
  int P0 = 2 * (1/ dL) * muM * M0;

  std::vector<int> output = {E0, L0, P0, M0};
  return(output);

}
