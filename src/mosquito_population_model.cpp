// Specifying all the includes and depends required to run the model
#include "Mosquito_Population_Model.hpp"
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]

//' @export
// [[Rcpp::export]]
Rcpp::List mosquito_population_model(int start_time, int end, Rcpp::NumericVector fitted_parameters, Rcpp::NumericVector static_parameters, std::vector<double> rainfall, Rcpp::String mortality_density_function) {

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
  Rcpp::NumericVector E_output(end_time - start_time); // why can't I use vector<double> without stating size?
  Rcpp::NumericVector L_output(end_time - start_time); // also when I was initialising these without the size, it was crashing automatically
  Rcpp::NumericVector P_output(end_time - start_time);
  Rcpp::NumericVector M_output(end_time - start_time);
  Rcpp::NumericVector k_output(end_time - start_time);
  int i = 0;
  Rcpp::NumericVector koutput(end);

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

    // Setting the density dependent function regulating larval mortality
    if (mortality_density_function == "power") {
      muE = muE0 * (1 + pow(((E + L) / (K)), dd_pow));
      muL = muL0 * (1 + (lambda * pow(((E + L) / (K)), dd_pow)));
    }
    else if (mortality_density_function == "exponential") {
      muE = muE0 * exp(((E + L) / (K)));
      muL = muL0 * exp(lambda * (((E + L) / (K))));
    }
    else if (mortality_density_function == "linear") {
      muE = muE0 * (1 + ((E + L) / (K)));
      muL = muL0 * (1 + lambda * (((E + L) / (K))));
    }

    // Specifying the total number of events occurring for each compartment at each timepoint

    // Be
    if ((dE + muE)*dt < 1) {
      Be = R::rbinom(E, (dE + muE)*dt);
    }
    else {
      Be = R::rbinom(E, 1);
    }

    // Bl
    if ((dL + muL)*dt < 1) {
      Bl = R::rbinom(L, (dL + muL)*dt);
    }
    else {
      Bl = R::rbinom(L, 1);
    }

    // Bp
    if ((dP + muP)*dt < 1) {
      Bp = R::rbinom(P, (dP + muP)*dt);
    }
    else {
      Bp = R::rbinom(P, 1);
    }

    // Bm (dealing with instances of 0 mosquitoes) - Mosquito Deaths
    if (M >= 1) { // Aaron has >= 1, surely should be > 1??
      Bm = R::rbinom(M, muM * dt);
    }
    else {
      Bm = 0;
    }

    // Updating the State Variables

    // E
    E = round(E - Be + M * beta * dt);
    if (E <= 0) {
      E = 0;
    }

    // L
    if (Be >= 1) {
      L = round(L - Bl + R::rbinom(Be, (dE / (muE + dE))));
    }
    else {
      L = round(L - Bl);
    }

    // P
    if (Bl >= 1) {
      P = round(P - Bp + R::rbinom(Bl, (dL / (muL + dL))));
    }
    else {
      P = round(P - Bp);
    }

    // M
    if (Bp > 1) {
      mRan = R::rbinom(Bp, (dP / (muP + dP)));
    }
    else {
      mRan = 0;
    }
    M = round(M + (0.5 * mRan) - Bm);
    if (M < 1) {
      M = 1;
    }

    E_output[i] = E;
    L_output[i] = L;
    P_output[i] = P;
    M_output[i] = M;

    // Stepping forward another timestep
    t = t + 1;
    i = i + 1;

  }

  return Rcpp::List::create(Rcpp::Named("E_Output") = E_output,
                            Rcpp::Named("L_Output") = L_output,
                            Rcpp::Named("P_Output") = P_output,
                            Rcpp::Named("M_Output") = M_output,
                            Rcpp::Named("K") = koutput);
}
