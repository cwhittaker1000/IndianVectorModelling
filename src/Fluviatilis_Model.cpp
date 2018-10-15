// Specifying all the includes and depends required to run the model
#include "Mosquito_Population_Model.hpp"
#include "Hill_Function.hpp"
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]

//' @export
// [[Rcpp::export]]
Rcpp::List mosquito_population_model_fluv(int start_time, int end, Rcpp::NumericVector fitted_parameters,
                                          Rcpp::NumericVector static_parameters,
                                          std::vector<double> rainfall,
                                          Rcpp::String mortality_density_function,
                                          Rcpp::String rainfall_relationship,
                                          Rcpp::String decline_type) {

  // Setting the Start and Endtime
  int t = start_time;
  int end_time = end;

  // Setting Various Fitted Model Parameters
  double dE = 1/fitted_parameters[0]; double dL = 1/fitted_parameters[1]; double dP = 1/fitted_parameters[2]; double muE0 = fitted_parameters[3];
  double muL0 = fitted_parameters[4]; double muP = fitted_parameters[5]; double muM = fitted_parameters[6];
  double lambda = fitted_parameters[7]; double tau = fitted_parameters[8]; double beta = fitted_parameters[9];
  double scaling_factor = fitted_parameters[12];

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

  // Fluviatilis Specific Parameters
  double Washout_Threshold = fitted_parameters[22];
  double washout_exp_scaling_factor = fitted_parameters[23];
  double washout_hill_one = fitted_parameters[24];
  double washout_hill_two = fitted_parameters[25];
  double marker = 0;
  double marker_stop;

  // Vectors to store the output of the model at each timestep
  Rcpp::NumericVector E_output(end_time - start_time); // why can't I use vector<double> without stating size?
  Rcpp::NumericVector L_output(end_time - start_time); // also when I was initialising these without the size, it was crashing automatically
  Rcpp::NumericVector P_output(end_time - start_time);
  Rcpp::NumericVector M_output(end_time - start_time);
  Rcpp::NumericVector k_output(end_time - start_time);
  Rcpp::NumericVector average_rainfall_K_Static(end);
  int i = 0;
  Rcpp::NumericVector k_total_output(end);

  // Iterating the model through multiple timepoints
  while (t < end_time) {

    // Specifying the value of K for all timepoints
    // Note: come back to this and change!!! Specifically, probably use rainfall data from beforehand to parameterise.
    if (t <= tau_with_dt) {

      if (rainfall_relationship == "mean") {
        K = 0;
        k_total_output[t] = 0;
      }
      else if (rainfall_relationship == "linear") {
        K = 0;
        k_total_output[t] = 0;
      }
      else if (rainfall_relationship == "exponential") {
        K = 0;
        k_total_output[t] = 0;
      }
    }

    else {

      if (rainfall_relationship == "mean") {

        if (decline_type == "exponential") {
          // Calculates Average Rainfall
          std::vector<double>::const_iterator first = rF.begin() + t - tau_with_dt;
          std::vector<double>::const_iterator last = rF.begin() + t;
          std::vector<double> rFx(first, last);
          rFsum = std::accumulate(rFx.begin(), rFx.end(), 0.0); // Don't forget to do 0.0!!! Otherwise accumulator will produce an int!!
          double rainfall_average = rFsum / tau_with_dt;
          average_rainfall_K_Static[t] = rainfall_average;
          // If Average Rainfall Exceeds Washout Threshold, Set Carrying Capacity to 0
          if (rainfall_average > Washout_Threshold) {
            K = 0;
            k_total_output[t] = K;
            marker = 0;
          }

          else {

            // First Time Washout Stops, Carrying Capacity Is At Max
            if (marker == 0) {
              K = K_static;
              k_total_output[t] = K;
              marker_stop = t;
              Rcpp::Rcout << "Washout stops at " << marker_stop << std::endl;
              marker = 1;
            }

            // Then Carrying Capacity Begins to Decline
            else {
              double calculation = -washout_exp_scaling_factor * (t - marker_stop);
              K = K_static * exp(calculation);
              k_total_output[t] = K;
            }
          }
        }

        else if (decline_type == "hill") {
          // Calculates Average Rainfall
          std::vector<double>::const_iterator first = rF.begin() + t - tau_with_dt;
          std::vector<double>::const_iterator last = rF.begin() + t;
          std::vector<double> rFx(first, last);
          rFsum = std::accumulate(rFx.begin(), rFx.end(), 0.0); // Don't forget to do 0.0!!! Otherwise accumulator will produce an int!!
          double rainfall_average = rFsum / tau_with_dt;

          // If Average Rainfall Exceeds Washout Threshold, Set Carrying Capacity to 0
          if (rainfall_average > Washout_Threshold) {
            K = 0;
            k_total_output[t] = K;
            marker = 0;
          }

          else {

            // First Time Washout Stops, Carrying Capacity Is At Max
            if (marker == 0) {
              K = K_static;
              k_total_output[t] = K;
              marker_stop = t;
              Rcpp::Rcout << "Washout stops at " << marker_stop << std::endl;
              marker = 1;
            }

            // Then Carrying Capacity Begins to Decline
            else {
              K = Hill_Function((t - marker_stop), K_static, washout_hill_one, washout_hill_two);
              k_total_output[t] = K;
            }
          }
        }
      }

      if (rainfall_relationship == "linear") {

        if (decline_type == "exponential") {
          double rFsum = 0.0;
          int counter = 1;

          for (int r = (t - tau_with_dt + 1); r <= t; r++) {
            rFsum = rFsum + (((t - tau_with_dt + counter) - t + tau_with_dt) * rF[r]);
            counter = counter + 1;
          }

          double rainfall_average = (2.0 / pow(tau_with_dt, 2)) * rFsum;

          // If Average Rainfall Exceeds Washout Threshold, Set Carrying Capacity to 0
          if (rainfall_average > Washout_Threshold) {
            K = 0;
            k_total_output[t] = K;
            marker = 0;
          }

          else {

            // First Time Washout Stops, Carrying Capacity Is At Max
            if (marker == 0) {
              K = K_static;
              k_total_output[t] = K;
              marker_stop = t;
              Rcpp::Rcout << "Washout stops at " << marker_stop << std::endl;
              marker = 1;
            }

            // Then Carrying Capacity Begins to Decline
            else {
              double calculation = -washout_exp_scaling_factor * (t - marker_stop);
              K = K_static * exp(calculation);
              k_total_output[t] = K;
            }
          }
        }
        else if (decline_type == "hill") {
          double rFsum = 0.0;
          int counter = 1;

          for (int r = (t - tau_with_dt + 1); r <= t; r++) {
            rFsum = rFsum + (((t - tau_with_dt + counter) - t + tau_with_dt) * rF[r]);
            counter = counter + 1;
          }

          double rainfall_average = (2.0 / pow(tau_with_dt, 2)) * rFsum;

          // If Average Rainfall Exceeds Washout Threshold, Set Carrying Capacity to 0
          if (rainfall_average > Washout_Threshold) {
            K = 0;
            k_total_output[t] = K;
            marker = 0;
          }

          else {

            // First Time Washout Stops, Carrying Capacity Is At Max
            if (marker == 0) {
              K = K_static;
              k_total_output[t] = K;
              marker_stop = t;
              Rcpp::Rcout << "Washout stops at " << marker_stop << std::endl;
              marker = 1;
            }

            // Then Carrying Capacity Begins to Decline
            else {
              K = Hill_Function((t - marker_stop), K_static, washout_hill_one, washout_hill_two);
              k_total_output[t] = K;
            }
          }
        }
      }

      if (rainfall_relationship == "exponential") {

        if (decline_type == "exponential") {
          double temp_tau = tau_with_dt;
          double rFsum = 0.0;
          double calc;
          for (int r = 0; r <= t; r++) {
            calc = (-(t - r))/temp_tau;
            rFsum = rFsum + exp(calc) * rF[r];
          }

          double rainfall_average = (1.0 / (temp_tau * (1 - exp(-(t + 1) / temp_tau)))) * rFsum;

          // If Average Rainfall Exceeds Washout Threshold, Set Carrying Capacity to 0
          if (rainfall_average > Washout_Threshold) {
            K = 0;
            k_total_output[t] = K;
            marker = 0;
          }

          else {

            // First Time Washout Stops, Carrying Capacity Is At Max
            if (marker == 0) {
              K = K_static;
              k_total_output[t] = K;
              marker_stop = t;
              Rcpp::Rcout << "Washout stops at " << marker_stop << std::endl;
              marker = 1;
            }

            // Then Carrying Capacity Begins to Decline
            else {
              double calculation = -washout_exp_scaling_factor * (t - marker_stop);
              K = K_static * exp(calculation);
              k_total_output[t] = K;
            }
          }
        }

        else if (decline_type == "hill") {
          double temp_tau = tau_with_dt;
          double rFsum = 0.0;
          double calc;
          for (int r = 0; r <= t; r++) {
            calc = (-(t - r))/temp_tau;
            rFsum = rFsum + exp(calc) * rF[r];
          }

          double rainfall_average = (1.0 / (temp_tau * (1 - exp(-(t + 1) / temp_tau)))) * rFsum;

          // If Average Rainfall Exceeds Washout Threshold, Set Carrying Capacity to 0
          if (rainfall_average > Washout_Threshold) {
            K = 0;
            k_total_output[t] = K;
            marker = 0;
          }

          else {

            // First Time Washout Stops, Carrying Capacity Is At Max
            if (marker == 0) {
              K = K_static;
              k_total_output[t] = K;
              marker_stop = t;
              Rcpp::Rcout << "Washout stops at " << marker_stop << std::endl;
              marker = 1;
            }

            // Then Carrying Capacity Begins to Decline
            else {
              K = Hill_Function((t - marker_stop), K_static, washout_hill_one, washout_hill_two);
              k_total_output[t] = K;
            }
          }
        }
      }
    }

    // Setting the density dependent function regulating larval mortality
    if (mortality_density_function == "power") {
      muE = muE0 * (1 + pow(((E + L) / (K)), dd_pow));
      muL = muL0 * (1 + (lambda * pow(((E + L) / K), dd_pow)));
    }
    else if (mortality_density_function == "exponential") {
      muE = muE0 * exp(((E + L) / K));
      muL = muL0 * exp(lambda * ((E + L) / K));
    }
    else if (mortality_density_function == "linear") {
      muE = muE0 * (1 + ((E + L) / K));
      muL = muL0 * (1 + lambda * ((E + L) / K));
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
                            Rcpp::Named("K_total") = k_total_output,
                            Rcpp::Named("rainfallaverage_Kstatic") = average_rainfall_K_Static);
}
