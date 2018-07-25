// Specifying all the includes and depends required to run the model
#include "Mosquito_Population_Model.hpp"
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector test_mosquito_population_model(int start_time, int end, Rcpp::NumericVector fitted_parameters, Rcpp::NumericVector static_parameters, std::vector<double> rainfall, Rcpp::String mortality_density_function) {

  // Setting the Start and Endtime
  int t = start_time;
  int end_time = end;

  // Setting Various Fitted Model Parameters
  double dE = fitted_parameters[0]; double dL = fitted_parameters[1]; double dP = fitted_parameters[2]; double muE0 = fitted_parameters[3];
  double muL0 = fitted_parameters[4]; double muP = fitted_parameters[5]; double muM = fitted_parameters[6];
  double lambda = fitted_parameters[7]; double tau = fitted_parameters[8]; double beta = fitted_parameters[9];

  // Setting the Static Parameters
  double dt = static_parameters[0];
  int dd_pow = static_parameters[1];

  // Setting the Initial State Variables
  std::vector<double> rF = rainfall;
  double E = fitted_parameters[12]; double L = fitted_parameters[13];
  double P = fitted_parameters[14]; double M = fitted_parameters[15];
  double Be = 0; double Bl = 0; double Bp = 0; double Bm = 0;

  // Declaring Additional Parameters
  double K; // ecological carrying capacity of the environment
  double muE; // mortality rate for Early larvae including density dependence
  double muL; // mortality rate for Late larvae including density dependence
  int tau_with_dt = tau / dt; // the number of timesteps of past rainfall that contribute to K, taking into account the timestep
  double rFsum; // rainfall summed over the tau days
  double rFsum_fudge; // rainfall summed over the tau days
  double mRan; // how many events out of the total pupal events (Bp) to assign as development into mosquitoes

  // Vectors to store the output of the model at each timestep
  Rcpp::NumericVector koutput(end);
  int i = 0;
  koutput[0] = 10;

  // Iterating the model through multiple timepoints
  while (t < end_time) {

    // Specifying the value of K for all timepoints
    if (t <= tau_with_dt) {

      std::vector<double>::const_iterator first_fudge = rF.begin();
      std::vector<double>::const_iterator last_fudge = rF.begin() + t;
      std::vector<double> rFx_fudge(first_fudge, last_fudge);
      rFsum_fudge = std::accumulate(rFx_fudge.begin(), rFx_fudge.end(), 0.0); // Don't forget to do 0.0!!! Otherwise accumulator will produce an int!!

      if (t == 0) {
        K = rainfall[t];
      }

      else {
        K = ((1.0 / t) * rFsum_fudge);
      }

      koutput[t] = K;

    }

    // if (t <= tau_with_dt) {
    //   K = (1.0 + (1.0 / tau_with_dt));
    // }

    else {

      std::vector<double>::const_iterator first = rF.begin() + t - tau_with_dt;
      std::vector<double>::const_iterator last = rF.begin() + t;
      std::vector<double> rFx(first, last);
      rFsum = std::accumulate(rFx.begin(), rFx.end(), 0.0); // Don't forget to do 0.0!!! Otherwise accumulator will produce an int!!
      K = (1.0 + ((1.0 / tau_with_dt) * rFsum));

      koutput[t] = K;
    }

  t = t + 1;
  }

  t = t + 1;
  return(koutput);
}
