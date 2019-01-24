// Specifying all the includes and depends required to run the model
#include "Mosquito_Population_Model.hpp"
#include "Hill_Function.hpp"
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]

//' @export
// [[Rcpp::export]]
Rcpp::List mosquito_population_model_cul(int start_time, int end, Rcpp::NumericVector fitted_parameters,
                                         Rcpp::NumericVector static_parameters,
                                         std::vector<double> rainfall,
                                         Rcpp::String mortality_density_function,
                                         Rcpp::String rainfall_relationship,
                                         Rcpp::String rainfall_effect) {

  // Setting the Start and Endtime
  int t = start_time;
  int end_time = end;

  // Setting Various Fitted Model Parameters
  double dE = 1/fitted_parameters[0]; double dL = 1/fitted_parameters[1]; double dP = 1/fitted_parameters[2]; double muE0 = fitted_parameters[3];
  double muL0 = fitted_parameters[4]; double muP = fitted_parameters[5]; double muM = fitted_parameters[6];
  double lambda = fitted_parameters[7]; double tau = fitted_parameters[8]; double beta = fitted_parameters[9];
  double scaling_factor = fitted_parameters[12]; double K_max = fitted_parameters[15]; double hill_1 = fitted_parameters[16]; double hill_2 = fitted_parameters[17];

  // Setting the Static Parameters
  double dt = static_parameters[0];
  int dd_pow = static_parameters[1];

  // Setting the Initial State Variables
  std::vector<double> rF = rainfall;
  double E = fitted_parameters[18]; double L = fitted_parameters[19];
  double P = fitted_parameters[20]; double M = fitted_parameters[21];
  double Be = 0; double Bl = 0; double Bp = 0; double Bm = 0;

  // Declaring Additional Parameters
  double K; // ecological carrying capacity of the environment
  double K_static = fitted_parameters[14] * scaling_factor;
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
  Rcpp::NumericVector k_rain_output(end);
  Rcpp::NumericVector k_total_output(end);

  // Iterating the model through multiple timepoints
  while (t < end_time) {

    // Specifying the value of K for all timepoints
    if (t <= tau_with_dt) { // note only mean and linear require separate rainfall calculations when t < tau_with_dt

      if (rainfall_relationship == "mean") {

        if (rainfall_effect == "raw") {

          std::vector<double>::const_iterator first_fudge = rF.begin();
          std::vector<double>::const_iterator last_fudge = rF.begin() + t;
          std::vector<double> rFx_fudge(first_fudge, last_fudge);
          rFsum_fudge = std::accumulate(rFx_fudge.begin(), rFx_fudge.end(), 0.0); // Don't forget to do 0.0!!! Otherwise accumulator will produce an int!!

          if (t == 0) {
            K = scaling_factor * rainfall[t];
          }

          else {
            K = scaling_factor * ((1.0 / t) * rFsum_fudge);
          }

          k_rain_output[t] = K;
          k_total_output[t] = K + K_static;

        }

        else if (rainfall_effect == "hill") {

          std::vector<double>::const_iterator first_fudge = rF.begin();
          std::vector<double>::const_iterator last_fudge = rF.begin() + t;
          std::vector<double> rFx_fudge(first_fudge, last_fudge);
          rFsum_fudge = std::accumulate(rFx_fudge.begin(), rFx_fudge.end(), 0.0); // Don't forget to do 0.0!!! Otherwise accumulator will produce an int!!

          if (t == 0) {

            double hill_output = Hill_Function(rainfall[t], K_max, hill_1, hill_2);
            K = scaling_factor * hill_output;

          }

          else {

            double mean_rainfall = (1.0 / t) * rFsum_fudge;
            double hill_output = Hill_Function(mean_rainfall, K_max, hill_1, hill_2);
            K = scaling_factor * hill_output;

          }

          k_rain_output[t] = K;
          k_total_output[t] = K + K_static;

        }

      }

      else if (rainfall_relationship == "linear") {

        if (rainfall_effect == "raw") {

          double rFsum = 0.0;
          for (int r = 0; r <= t; r++) {
            rFsum = rFsum + ((r + 1) * rF[r]);
          }

          if (t == 0) {
            K = scaling_factor * rF[t];
          }
          else {
            K = scaling_factor * rF[t]; // NOTE NEED TO STILL SORT THIS OUT
            // K = ((2.0 / pow(t + 1, 2)) * rFsum);
          }

          k_rain_output[t] = K;
          k_total_output[t] = K + K_static;

        }

        else if (rainfall_effect == "hill") {

          double rFsum = 0.0;
          for (int r = 0; r <= t; r++) {
            rFsum = rFsum + ((r + 1) * rF[r]);
          }

          if (t == 0) {
            double hill_output = Hill_Function(rF[t], K_max, hill_1, hill_2);
            K = scaling_factor * hill_output;
          }
          else { // NOTE NEED TO SORT THIS OUT STILL.
            double hill_output = Hill_Function(rF[t], K_max, hill_1, hill_2);
            K = scaling_factor * hill_output;
            // double linearly_weighted_rainfall = ((2.0 / pow(t + 1, 2)) * rFsum);
            // double hill_output = Hill_Function(linearly_weighted_rainfall, K_max, hill_1, hill_2);
            // K = scaling_factor * hill_output;
          }

          k_rain_output[t] = K;
          k_total_output[t] = K + K_static;

        }
      }

      else if (rainfall_relationship == "exponential") {

        if (rainfall_effect == "raw") {

          double temp_tau = tau_with_dt; // to convert tau_with_dt to a double and ensure downstream outputs are doubles
          double rFsum = 0.0;
          double calc;

          for (int r = 0; r <= t; r++) {

            calc = (-(t - r))/temp_tau;
            rFsum = rFsum + exp(calc) * rF[r];

          }

          K = scaling_factor * (1.0 / (temp_tau * (1 - exp(-(t + 1) / temp_tau)))) * rFsum;
          k_rain_output[t] = K;
          k_total_output[t] = K + K_static;

        }

        else if (rainfall_effect == "hill") {

          double temp_tau = tau_with_dt; // to convert tau_with_dt to a double and ensure downstream outputs are doubles
          double rFsum = 0.0;
          double calc;

          for (int r = 0; r <= t; r++) {

            calc = (-(t - r))/temp_tau;
            rFsum = rFsum + exp(calc) * rF[r];

          }

          double mean_rainfall = (1.0 / (temp_tau * (1 - exp(-(t + 1) / temp_tau)))) * rFsum;
          double hill_output = Hill_Function(mean_rainfall, K_max, hill_1, hill_2);
          K = scaling_factor * hill_output;
          k_rain_output[t] = K;
          k_total_output[t] = K + K_static;

        }
      }
    }

    else {

      if (rainfall_relationship == "mean") {

        if (rainfall_effect == "raw") {

          std::vector<double>::const_iterator first = rF.begin() + t - tau_with_dt;
          std::vector<double>::const_iterator last = rF.begin() + t;
          std::vector<double> rFx(first, last);
          rFsum = std::accumulate(rFx.begin(), rFx.end(), 0.0); // Don't forget to do 0.0!!! Otherwise accumulator will produce an int!!
          K = scaling_factor * ((1.0 / tau_with_dt) * rFsum);

          k_rain_output[t] = K;
          k_total_output[t] = K + K_static;

        }

        else if (rainfall_effect == "hill") {

          std::vector<double>::const_iterator first = rF.begin() + t - tau_with_dt;
          std::vector<double>::const_iterator last = rF.begin() + t;
          std::vector<double> rFx(first, last);
          rFsum = std::accumulate(rFx.begin(), rFx.end(), 0.0); // Don't forget to do 0.0!!! Otherwise accumulator will produce an int!!

          double mean_rainfall = ((1.0 / tau_with_dt) * rFsum);
          double hill_output = Hill_Function(mean_rainfall, K_max, hill_1, hill_2);
          K = scaling_factor * hill_output;

          k_rain_output[t] = K;
          k_total_output[t] = K + K_static;

        }

      }

      else if (rainfall_relationship == "linear") {

        if (rainfall_effect == "raw") {

          double rFsum = 0.0;
          int counter = 1;
          for (int r = (t - tau_with_dt + 1); r <= t; r++) {

            // THIS I think takes the rainfall up to the day before t: I start with t = 0, and
            // so t = 1 is actually day 2 etc and this needs to be accounted for given that
            // the model I'm exploring currently includes tau days previous INCLUDING present day (as in Michael's)
            // for (int r = (t - tau_with_dt); r < t; r++) {
            // rFsum = rFsum + (((t - tau_with_dt + counter) - t + tau_with_dt) * rF[r]);

            // This should be correct
            rFsum = rFsum + (((t - tau_with_dt + counter) - t + tau_with_dt) * rF[r]);
            counter = counter + 1;
          }

          K = scaling_factor * ((2.0 / pow(tau_with_dt, 2)) * rFsum);
          k_rain_output[t] = K;
          k_total_output[t] = K + K_static;

        }

        else if (rainfall_effect == "hill") {

          double rFsum = 0.0;
          int counter = 1;
          for (int r = (t - tau_with_dt + 1); r <= t; r++) {

            // THIS I think takes the rainfall up to the day before t: I start with t = 0, and
            // so t = 1 is actually day 2 etc and this needs to be accounted for given that
            // the model I'm exploring currently includes tau days previous INCLUDING present day (as in Michael's)
            // for (int r = (t - tau_with_dt); r < t; r++) {
            // rFsum = rFsum + (((t - tau_with_dt + counter) - t + tau_with_dt) * rF[r]);

            // This should be correct
            rFsum = rFsum + (((t - tau_with_dt + counter) - t + tau_with_dt) * rF[r]);
            counter = counter + 1;
          }

          double mean_rainfall = ((2.0 / pow(tau_with_dt, 2)) * rFsum);
          double hill_output = Hill_Function(mean_rainfall, K_max, hill_1, hill_2);
          K = scaling_factor * hill_output;

          k_rain_output[t] = K;
          k_total_output[t] = K + K_static;

        }
      }

      else if (rainfall_relationship == "exponential") {

        if (rainfall_effect == "raw") {

          double temp_tau = tau_with_dt;
          double rFsum = 0.0;
          double calc;

          for (int r = 0; r <= t; r++) {

            calc = (-(t - r))/temp_tau;
            rFsum = rFsum + exp(calc) * rF[r];

          }

          K = scaling_factor * (1.0 / (temp_tau * (1 - exp(-(t + 1) / temp_tau)))) * rFsum;
          k_rain_output[t] = K;
          k_total_output[t] = K + K_static;

        }

        else if (rainfall_effect == "hill") {

          double temp_tau = tau_with_dt;
          double rFsum = 0.0;
          double calc;

          for (int r = 0; r <= t; r++) {

            calc = (-(t - r))/temp_tau;
            rFsum = rFsum + exp(calc) * rF[r];

          }

          double mean_rainfall = (1.0 / (temp_tau * (1 - exp(-(t + 1) / temp_tau)))) * rFsum;
          double hill_output = Hill_Function(mean_rainfall, K_max, hill_1, hill_2);
          K = scaling_factor * hill_output;

          k_rain_output[t] = K;
          k_total_output[t] = K + K_static;

        }
      }
    }

    // Setting the density dependent function regulating larval mortality
    if (mortality_density_function == "power") {
      muE = muE0 * (1 + pow(((E + L) / (K + K_static)), dd_pow));
      muL = muL0 * (1 + (lambda * pow(((E + L) / (K + K_static)), dd_pow)));
    }
    else if (mortality_density_function == "exponential") {
      muE = muE0 * exp(((E + L) / (K + K_static)));
      muL = muL0 * exp(lambda * ((E + L) / (K + K_static)));
    }
    else if (mortality_density_function == "linear") {
      muE = muE0 * (1 + ((E + L) / (K + K_static)));
      muL = muL0 * (1 + lambda * ((E + L) / (K + K_static)));
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
                            Rcpp::Named("K_rain") = k_rain_output,
                            Rcpp::Named("K_total") = k_total_output);
}
