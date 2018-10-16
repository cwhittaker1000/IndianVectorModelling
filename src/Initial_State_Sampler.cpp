#include "Initial_State_Sampler.hpp"
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]

// Linear Initital State Sampler
//' @export
// [[Rcpp::export]]
std::vector <int> test_initial_state_sample(Rcpp::NumericVector fitted_parameters,
                                            Rcpp::NumericVector static_parameters,
                                            Rcpp::String mortality_density_function) {

  double dE = 1/fitted_parameters[0]; double dL = 1/fitted_parameters[1]; double dP = 1/fitted_parameters[2];
  double muE0 = fitted_parameters[3]; double muL0 = fitted_parameters[4]; double muP = fitted_parameters[5]; double muM = fitted_parameters[6];
  double lambda = fitted_parameters[7]; double beta = fitted_parameters[8];
  double z = fitted_parameters[11];

  // Various Components Involved in the Initial States Calculation for Exponential Density Dependence Governing Mortality - NOT USING
  // double scaled_static_K = fitted_parameters[12] * fitted_parameters[14]; // fitted_parameters[14] = K_static, fitted_parameters[12] = scaling factor
  // double K = initial_K + scaled_static_K;
  // double dd_pow = static_parameters[1]; Not needed just yet but might be needed at some point
  // double muE = muE0 * exp((round(z) / K));
  // double muL = muL0 * exp((lambda * round(z) / K));

  // Initialising initial conditions for storage of output
  int E; int L; int P; int M;

  if (mortality_density_function == "power") { // needs to be modified to include dd_pow

    // Power Density Dependence
    double a = (beta * dP * dL) / ((2 * muM) * (muP + dP));
    double b = (muE0 / (lambda * muL0)) * (dL + muL0) - dE - muE0;
    double c = -(muE0 * dE) / (muL0 * lambda);
    double x = (-b + sqrt(pow(b, 2) - 4 * a * c)) / (2 * a);

    // Initial Conditions
    L = round(z);
    E = round(L / x);
    P = round((dL * L) / (muP + dP));
    M = round((dP * P) / (2 * muM));

  }

  // else if (mortality_density_function == "exponential") {
  //
  //   //Exponential Density Dependence
  //   double a = (0.5 * dL * dP) / (muM * (dP + muP));
  //
  //   // Initial Conditions
  //   E = round(((beta * a) * round(z)) / ((dE + (beta * a) + muE)));
  //   L = round((dE * round(z)) / ((dE + dL + muL)));
  //   M = round((0.5 * dL * dP * L) / (muM * (dP + muP)));
  //   P = round((2 * muM * M) / dP);
  //
  // }

  else if (mortality_density_function == "linear") {

    // Linear Density Dependence
    double a = (beta * dP * dL) / ((2 * muM) * (muP + dP));
    double b = (muE0 / (lambda * muL0)) * (dL + muL0) - dE - muE0;
    double c = -(muE0 * dE) / (muL0 * lambda);
    double x = (-b + sqrt(pow(b, 2) - 4 * a * c)) / (2 * a);

    // Initial Conditions
    L = round(z);
    E = round(L / x);
    P = round((dL * L) / (muP + dP));
    M = round((dP * P) / (2 * muM));

  }

  std::vector<int> output = {E, L, P, M};
  return(output);

}





