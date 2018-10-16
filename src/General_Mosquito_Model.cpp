// Specifying all the includes and depends required to run the model
#include "Mosquito_Population_Model.hpp"
#include "Hill_Function.hpp"
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]

//' @export
// [[Rcpp::export]]
Rcpp::List general_mosquito_population_model(int start_time, int end,
                                             Rcpp::NumericVector fitted_parameters, Rcpp::NumericVector static_parameters,
                                             std::vector<double> rainfall,
                                             Rcpp::String mortality_density_function,
                                             Rcpp::String rainfall_relationship,
                                             Rcpp::String rainfall_effect,
                                             Rcpp::String decline_type) {

  // Setting the Start and Endtime
  int t = start_time;
  int end_time = end;

  // Setting the Initial State Variables
  double E = fitted_parameters[24]; double L = fitted_parameters[25];
  double P = fitted_parameters[26]; double M = fitted_parameters[27];
  double Be = 0; double Bl = 0; double Bp = 0; double Bm = 0;

  // Constant Parameters
  double dt = static_parameters[0];   int dd_pow = static_parameters[1]; std::vector<double> rF = rainfall; double offset = fitted_parameters[28];

  // Core Parameters
  double dE = 1/fitted_parameters[0]; double dL = 1/fitted_parameters[1]; double dP = 1/fitted_parameters[2]; double muE0 = fitted_parameters[3];
  double muL0 = fitted_parameters[4]; double muP = fitted_parameters[5]; double muM = fitted_parameters[6];
  double lambda = fitted_parameters[7]; double beta = fitted_parameters[8]; // fitted_parameters[9], [10] and [11] are overdisp, pop_frac and z and these don't directly feature in the model at this stage.

  // Parameters Governing K_Rain(t)
  double tau_rain = fitted_parameters[12];
  double scaling_factor_rainfall = fitted_parameters[13];
  double K_Max_Hill_Rainfall = fitted_parameters[14]; // If a Hill Function is Used, Rather Than the Raw Rainfall
  double Hill_Rainfall_1 = fitted_parameters[15]; // If a Hill Function is Used, Rather Than the Raw Rainfall
  double Hill_Rainfall_2 = fitted_parameters[16]; // If a Hill Function is Used, Rather Than the Raw Rainfall

  // Parameters Governing K_Static(t)
  double tau_static = fitted_parameters[17]; // Number of days rainfall contributing to washout occurrence calculations
  double scaling_factor_static = fitted_parameters[18];
  double K_Max_Static = fitted_parameters[19] * scaling_factor_static; // Max value of K_Static
  double Washout_Threshold = fitted_parameters[20]; // Threshold at which washout starts occurring
  double washout_exp_scaling_factor = fitted_parameters[21]; // Governing exponential decline of K_Static in dry season
  double washout_hill_one = fitted_parameters[22]; // Governing Hill function based decline of K_Static in dry season
  double washout_hill_two = fitted_parameters[23]; // Governing Hill function based decline of K_Static in dry season
  double marker = 0; // Involved in calculating the last timepoint at which washout occurred
  double marker_stop; // Involved in calculating the last timepoint at which washout occurred

  // Initialising Variables Required in the Model
  double K_Static; // permanent component of the environmental carrying capacity
  double K_Rain; // rainfall responsive component of the environmental carrying capacity
  double muE; // mortality rate for Early larvae including density dependence
  double muL; // mortality rate for Late larvae including density dependence
  int tau_with_dt_rain = tau_rain / dt; // the number of timesteps of past rainfall that contribute to K_Rain, taking into account the timestep
  int tau_with_dt_static = tau_static / dt; // the number of timesteps of past rainfall that contribute to K_Static washout calculations, taking into account the timestep
  double rFsum_rain; // rainfall summed over the tau_rain days
  double rFsum_static; // rainfall summed over the tau_static days
  double rFsum_static_prior; // rainfall summed over the previous tau_static days, but for the rainfall before data start
  double mRan; // how many events out of the total pupal events (Bp) to assign as development into mosquitoes

  // Vectors to store the output of the model at each timestep
  Rcpp::NumericVector E_output(end_time - start_time); Rcpp::NumericVector L_output(end_time - start_time);
  Rcpp::NumericVector P_output(end_time - start_time); Rcpp::NumericVector M_output(end_time - start_time);
  Rcpp::NumericVector k_output(end_time - start_time); Rcpp::NumericVector k_rain_output(end);
  Rcpp::NumericVector k_static_output(end); Rcpp::NumericVector k_total_output(end);
  Rcpp::NumericVector average_rainfall_K_Static(end);

  Rcpp::NumericVector prior_data_calculations(start_time + offset - tau_with_dt_static + 1);
  int i = 0;
  int prior_counter = 0;

  // Code to Calculate the Time at Which Washout Last Occurred In the Preceding, Offset Encompassed Rainfall
  for (int p = tau_with_dt_static; p < (offset + start_time); p++) {

    if (rainfall_relationship == "mean") {

      // For Each Timepoint to Start of Data, Calculates the Summed Rainfall
      std::vector<double>::const_iterator first_static_prior = rF.begin() + p - tau_with_dt_static;
      std::vector<double>::const_iterator last_static_prior = rF.begin() + p;
      std::vector<double> rFx_static_prior(first_static_prior, last_static_prior);
      rFsum_static_prior = std::accumulate(rFx_static_prior.begin(), rFx_static_prior.end(), 0.0);

      // Calculating Whether Washout Occurs
      double rainfall_average_static_calc_prior = rFsum_static_prior / tau_with_dt_static;
      if (rainfall_average_static_calc_prior > Washout_Threshold) {
        marker = 0;
      }
      else {
        if (marker == 0) { // First Time Washout Stops, Carrying Capacity Is At Max
          K_Static = K_Max_Static;
          marker_stop = p;
          marker = 1;
        }
        else { // Then Carrying Capacity Begins to Decline
          if (decline_type == "exponential") {
            double calculation = -washout_exp_scaling_factor * (p - marker_stop);
            K_Static = K_Max_Static * exp(calculation);
          }
          else {
            K_Static = Hill_Function((p - marker_stop), K_Max_Static, washout_hill_one, washout_hill_two);
          }
        }
      }
      prior_data_calculations[prior_counter] = K_Static;
      prior_counter = prior_counter + 1;
    }

    else if (rainfall_relationship == "linear") {
      double rFsum_statically_prior = 0.0;
      int counter_static_prior = 1;
      for (int r = (p - tau_with_dt_static); r <= p; r++) {
        rFsum_statically_prior = rFsum_statically_prior + (((p - tau_with_dt_static + counter_static_prior) - p + tau_with_dt_static) * rF[r]);
        counter_static_prior = counter_static_prior + 1;
      }
      // Calculate Whether Washout Occurs
      double rainfall_average_static_calc_prior = ((2.0 / pow(tau_with_dt_static, 2)) * rFsum_statically_prior); // unsure if this is what I want, double check
      if (rainfall_average_static_calc_prior > Washout_Threshold) {
        marker = 0;
      }
      else {
        if (marker == 0) { // First Time Washout Stops, Carrying Capacity Is At Max
          K_Static = K_Max_Static;
          marker_stop = p;
          marker = 1;
        }
        else { // Then Carrying Capacity Begins to Decline
          if (decline_type == "exponential") {
            double calculation = -washout_exp_scaling_factor * (p - marker_stop);
            K_Static = K_Max_Static * exp(calculation);
          }
          else {
            K_Static = Hill_Function((p - marker_stop), K_Max_Static, washout_hill_one, washout_hill_two);
          }
        }
      }
      prior_data_calculations[prior_counter] = K_Static;
      prior_counter = prior_counter + 1;
    }
    else if (rainfall_relationship == "exponential") {

      double temp_tau_static_prior = tau_with_dt_static;
      double rFsum_statically_prior = 0.0;
      double calc;
      for (int r = 0; r <= p; r++) {
        calc = (-(p - r))/temp_tau_static_prior;
        rFsum_statically_prior = rFsum_statically_prior + exp(calc) * rF[r];
      }
      // Calculating Whether Washout Occurs
      double rainfall_average_static_calc_prior = (1.0 / (temp_tau_static_prior * (1 - exp(-(t + offset + 1) / temp_tau_static_prior)))) * rFsum_statically_prior;
      if (rainfall_average_static_calc_prior > Washout_Threshold) {
        marker = 0;
      }
      else {
        if (marker == 0) { // First Time Washout Stops, Carrying Capacity Is At Max
          K_Static = K_Max_Static;
          marker_stop = p;
          marker = 1;
        }
        else { // Then Carrying Capacity Begins to Decline
          if (decline_type == "exponential") {
            double calculation = -washout_exp_scaling_factor * (p - marker_stop);
            K_Static = K_Max_Static * exp(calculation);
          }
          else {
            K_Static = Hill_Function((p - marker_stop), K_Max_Static, washout_hill_one, washout_hill_two);
          }
        }
      }
      prior_data_calculations[prior_counter] = K_Static;
      prior_counter = prior_counter + 1;
    }
  }

  //Rcpp::Rcout << "The value of K_Static at this stage is " << K_Static << std::endl;
  //Rcpp::Rcout << "The value of Marker_Stop at this stage is " << marker_stop << std::endl;

  int within_data_marker = 0;

  // Iterating the model through multiple timepoints
  while (t < end_time) {

    if (t == 0) {
      //Rcpp::Rcout << "The value of K_Static at this next stage is " << K_Static << std::endl;
      //Rcpp::Rcout << "The value of Marker_Stop at this next stage is " << marker_stop << std::endl;
    }

    // Specifying the value of K for all timepoints
    if (rainfall_relationship == "mean") {

      // Calculating K_Rain
      if (rainfall_effect == "raw") {
        // Calculating Raw Rainfall Average
        std::vector<double>::const_iterator first_rain = rF.begin() + t + offset- tau_with_dt_rain;
        std::vector<double>::const_iterator last_rain = rF.begin() + t + offset;
        std::vector<double> rFx_rain(first_rain, last_rain);
        rFsum_rain = std::accumulate(rFx_rain.begin(), rFx_rain.end(), 0.0); // Don't forget to do 0.0!!! Otherwise accumulator will produce an int!!
        K_Rain = scaling_factor_rainfall * (1.0 / tau_with_dt_rain) * rFsum_rain;
      }
      else if (rainfall_effect == "hill") {
        // Calculating Hill Function'd Raw Average Rainfall
        std::vector<double>::const_iterator first_rain = rF.begin() + t + offset- tau_with_dt_rain;
        std::vector<double>::const_iterator last_rain = rF.begin() + t + offset;
        std::vector<double> rFx_rain(first_rain, last_rain);
        rFsum_rain = std::accumulate(rFx_rain.begin(), rFx_rain.end(), 0.0); // Don't forget to do 0.0!!! Otherwise accumulator will produce an int!!
        double mean_rainfall = (1.0 / tau_with_dt_rain) * rFsum_rain;
        double hill_output = Hill_Function(mean_rainfall, K_Max_Hill_Rainfall, Hill_Rainfall_1, Hill_Rainfall_2);
        K_Rain = scaling_factor_rainfall * hill_output;
      }
      // Calculating K_Static
      std::vector<double>::const_iterator first_static = rF.begin() + t + offset - tau_with_dt_static;
      std::vector<double>::const_iterator last_static = rF.begin() + t + offset;
      std::vector<double> rFx_static(first_static, last_static);
      rFsum_static = std::accumulate(rFx_static.begin(), rFx_static.end(), 0.0);
      // Calculating Whether Washout Occurs
      double rainfall_average_static_calc = rFsum_static / tau_with_dt_static;
      average_rainfall_K_Static[t] = rainfall_average_static_calc;
      if (rainfall_average_static_calc > Washout_Threshold) {
        K_Static = 0;
        marker = 0;
        within_data_marker = 1;
      }
      else {
        if (marker == 0) { // First Time Washout Stops, Carrying Capacity Is At Max
          K_Static = K_Max_Static;
          marker_stop = t;
          marker = 1;
        }
        else { // Then Carrying Capacity Begins to Decline
          if (decline_type == "exponential") {
            if (within_data_marker == 0) {
              double calculation = -washout_exp_scaling_factor * (t + offset - marker_stop);
              K_Static = K_Max_Static * exp(calculation);
            }
            else {
              double calculation = -washout_exp_scaling_factor * (t - marker_stop);
              K_Static = K_Max_Static * exp(calculation);
            }
          }
          else {
            if (within_data_marker == 0) {
              K_Static = Hill_Function((t + offset - marker_stop), K_Max_Static, washout_hill_one, washout_hill_two);
            }
            else {
              K_Static = Hill_Function((t - marker_stop), K_Max_Static, washout_hill_one, washout_hill_two);
            }
          }
        }
      }
      k_static_output[t] = K_Static;
      k_rain_output[t] = K_Rain;
      k_total_output[t] = K_Static + K_Rain;
    }

    else if (rainfall_relationship == "linear") {
      // Calculating K_Rain
      if (rainfall_effect == "raw") {
        // Calculating Linearly Weighted Average Rainfall
        double rFsum_rainfall = 0.0;
        int counter_rain = 1;
        for (int r = (t + offset - tau_with_dt_rain + 1); r <= t + offset; r++) {
          rFsum_rainfall = rFsum_rainfall + (((t - tau_with_dt_rain + counter_rain) - t + tau_with_dt_rain) * rF[r]);
          counter_rain = counter_rain + 1;
        }
        K_Rain = scaling_factor_rainfall * (2.0 / pow(tau_with_dt_rain, 2)) * rFsum_rainfall;
      }
      else if (rainfall_effect == "hill") {
        // Calculating Hill Function'd, Linearly Weighted Average Rainfall
        double rFsum_rainfall = 0.0;
        int counter_rain = 1;
        for (int r = (t + offset - tau_with_dt_rain + 1); r <= t + offset; r++) {
          rFsum_rainfall = rFsum_rainfall + (((t - tau_with_dt_rain + counter_rain) - t + tau_with_dt_rain) * rF[r]);
          counter_rain = counter_rain + 1;
        }
        double mean_rainfall = ((2.0 / pow(tau_with_dt_rain, 2)) * rFsum_rainfall);
        double hill_output = Hill_Function(mean_rainfall, K_Max_Hill_Rainfall, Hill_Rainfall_1, Hill_Rainfall_2);
        K_Rain = scaling_factor_rainfall * hill_output;
      }

      // Calculating K_Static
      double rFsum_statically = 0.0;
      int counter_static = 1;
      for (int r = (t + offset- tau_with_dt_static + 1); r <= t + offset; r++) {
        rFsum_statically = rFsum_statically + (((t - tau_with_dt_static + counter_static) - t + tau_with_dt_static) * rF[r]);
        counter_static = counter_static + 1;
      }
      // Calculate Whether Washout Occurs
      double rainfall_average_static_calc = ((2.0 / pow(tau_with_dt_static, 2)) * rFsum_statically); // unsure if this is what I want, double check
      average_rainfall_K_Static[t] = rainfall_average_static_calc;
      if (rainfall_average_static_calc > Washout_Threshold) {
        K_Static = 0;
        marker = 0;
        within_data_marker = 1;
      }
      else {
        if (marker == 0) { // First Time Washout Stops, Carrying Capacity Is At Max
          K_Static = K_Max_Static;
          marker_stop = t;
          marker = 1;
        }
        else { // Then Carrying Capacity Begins to Decline
          if (decline_type == "exponential") {
            if (within_data_marker == 0) {
              double calculation = -washout_exp_scaling_factor * (t + offset - marker_stop);
              K_Static = K_Max_Static * exp(calculation);
            }
            else {
              double calculation = -washout_exp_scaling_factor * (t - marker_stop);
              K_Static = K_Max_Static * exp(calculation);
            }
          }
          else {
            if (within_data_marker == 0) {
              K_Static = Hill_Function((t + offset - marker_stop), K_Max_Static, washout_hill_one, washout_hill_two);
            }
            else {
              K_Static = Hill_Function((t - marker_stop), K_Max_Static, washout_hill_one, washout_hill_two);
            }
          }
        }
      }
      k_static_output[t] = K_Static;
      k_rain_output[t] = K_Rain;
      k_total_output[t] = K_Static + K_Rain;
    }

    else if (rainfall_relationship == "exponential") {
      // Calculating K_Rain
      if (rainfall_effect == "raw") {
        // Calculating Exponentially Weighted Rainfall Average
        double temp_tau_rain = tau_with_dt_rain;
        double rFsum_rainfall = 0.0;
        double calc;
        for (int r = 0; r <= t + offset; r++) { // need to decide whether to start at beginning of rain or t = 0 ??
          calc = (-(t + offset - r))/temp_tau_rain;
          rFsum_rainfall = rFsum_rainfall + exp(calc) * rF[r];
        }
        K_Rain = scaling_factor_rainfall * (1.0 / (temp_tau_rain * (1 - exp(-(t + offset + 1) / temp_tau_rain)))) * rFsum_rainfall;
      }
      else if (rainfall_effect == "hill") {
        // Calculating Hill Function'd Exponentially Weighted Rainfall Average
        double temp_tau_rain = tau_with_dt_rain;
        double rFsum_rainfall = 0.0;
        double calc;
        for (int r = 0; r <= t + offset; r++) {
          calc = (-(t + offset - r))/temp_tau_rain;
          rFsum_rainfall = rFsum_rainfall + exp(calc) * rF[r];
        }
        double mean_rainfall = (1.0 / (temp_tau_rain * (1 - exp(-(t + offset + 1) / temp_tau_rain)))) * rFsum_rainfall;
        double hill_output = Hill_Function(mean_rainfall, K_Max_Hill_Rainfall, Hill_Rainfall_1, Hill_Rainfall_2);
        K_Rain = scaling_factor_rainfall * hill_output;
      }

      // Calculating K_Static
      double temp_tau_static = tau_with_dt_static;
      double rFsum_statically = 0.0;
      double calc;
      for (int r = 0; r <= t + offset; r++) {
        calc = (-(t + offset - r))/temp_tau_static;
        rFsum_statically = rFsum_statically + exp(calc) * rF[r];
      }
      // Calculating Whether Washout Occurs
      double rainfall_average_static_calc = (1.0 / (temp_tau_static * (1 - exp(-(t + offset + 1) / temp_tau_static)))) * rFsum_statically;
      average_rainfall_K_Static[t] = rainfall_average_static_calc;
      if (rainfall_average_static_calc > Washout_Threshold) {
        K_Static = 0;
        marker = 0;
        within_data_marker = 1;
      }
      else {
        if (marker == 0) { // First Time Washout Stops, Carrying Capacity Is At Max
          K_Static = K_Max_Static;
          marker_stop = t;
          marker = 1;
        }
        else { // Then Carrying Capacity Begins to Decline
          if (decline_type == "exponential") {
            if (within_data_marker == 0) {
              double calculation = -washout_exp_scaling_factor * (t + offset - marker_stop);
              K_Static = K_Max_Static * exp(calculation);
            }
            else {
              double calculation = -washout_exp_scaling_factor * (t - marker_stop);
              K_Static = K_Max_Static * exp(calculation);
            }
          }
          else {
            if (within_data_marker == 0) {
              K_Static = Hill_Function((t + offset - marker_stop), K_Max_Static, washout_hill_one, washout_hill_two);
            }
            else {
              K_Static = Hill_Function((t - marker_stop), K_Max_Static, washout_hill_one, washout_hill_two);
            }
          }
        }
      }
      k_static_output[t] = K_Static;
      k_rain_output[t] = K_Rain;
      k_total_output[t] = K_Static + K_Rain;
    }

    // Setting the density dependent function regulating larval mortality
    if (mortality_density_function == "power") {
      muE = muE0 * (1 + pow(((E + L) / (K_Rain + K_Static)), dd_pow));
      muL = muL0 * (1 + (lambda * pow(((E + L) / (K_Rain + K_Static)), dd_pow)));
    }
    else if (mortality_density_function == "exponential") {
      muE = muE0 * exp(((E + L) / (K_Rain + K_Static)));
      muL = muL0 * exp(lambda * ((E + L) / (K_Rain + K_Static)));
    }
    else if (mortality_density_function == "linear") {
      muE = muE0 * (1 + ((E + L) / (K_Rain + K_Static)));
      muL = muL0 * (1 + lambda * ((E + L) / (K_Rain + K_Static)));
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

    // Adding the Current State Variables to the Output
    E_output[i] = E; L_output[i] = L; P_output[i] = P; M_output[i] = M;

    // Stepping Forward Another Timestep
    t = t + 1;
    i = i + 1; // can't remember why I have this. Make sure to check
  }

  return Rcpp::List::create(Rcpp::Named("E_Output") = E_output,
                            Rcpp::Named("L_Output") = L_output,
                            Rcpp::Named("P_Output") = P_output,
                            Rcpp::Named("M_Output") = M_output,
                            Rcpp::Named("K_Rain") = k_rain_output,
                            Rcpp::Named("K_Static") = k_static_output,
                            Rcpp::Named("K_Total") = k_total_output,
                            Rcpp::Named("rainfallaverage_Kstatic") = average_rainfall_K_Static,
                            Rcpp::Named("prior K") = prior_data_calculations);
}
