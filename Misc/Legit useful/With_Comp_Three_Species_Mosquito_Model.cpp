// Specifying all the includes and depends required to run the model
#include "Mosquito_Population_Model.hpp"
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]

//' @export
// [[Rcpp::export]]
Rcpp::List three_species_mosquito_population_model_with_comp_k_part(int start_time, int end,
                                                                  Rcpp::NumericVector fitted_parameters_species_one, Rcpp::NumericVector fitted_parameters_species_two,
                                                                  Rcpp::NumericVector fitted_parameters_species_three, Rcpp::NumericVector static_parameters,
                                                                  std::vector<double> rainfall,
                                                                  Rcpp::String mortality_density_function_species_one, Rcpp::String mortality_density_function_species_two, Rcpp::String mortality_density_function_species_three,
                                                                  Rcpp::String rainfall_relationship_species_one, Rcpp::String rainfall_relationship_species_two, Rcpp::String rainfall_relationship_species_three,
                                                                  double competition_parameter_two_on_one, double competition_parameter_three_on_one,
                                                                  double competition_parameter_one_on_two, double competition_parameter_three_on_two,
                                                                  double competition_parameter_one_on_three, double competition_parameter_two_on_three,
                                                                  double species_one_static_K, double species_two_static_K, double species_three_static_K,
                                                                  double threshold_dd_1, double threshold_dd_2, double threshold_dd_3, // currently dd thresholds same for the 2 species exerting competitive influence on 3rd. Consider changing so it can be different.
                                                                  double max_K_1, double max_K_2, double max_K_3) {

  // Setting the Start and Endtime
  int t = start_time;
  int end_time = end;

  // Setting Various Fitted Model Parameters

  // Species One
  double dE_1 = 1/fitted_parameters_species_one[0]; double  dL_1 = 1/fitted_parameters_species_one[1]; double dP_1 = 1/fitted_parameters_species_one[2];
  double muE0_1 = fitted_parameters_species_one[3]; double muL0_1 = fitted_parameters_species_one[4]; double muP_1 = fitted_parameters_species_one[5];
  double muM_1 = fitted_parameters_species_one[6]; double lambda_1 = fitted_parameters_species_one[7]; double tau_1 = fitted_parameters_species_one[8];
  double beta_1 = fitted_parameters_species_one[9];

  // Species Two
  double dE_2 = 1/fitted_parameters_species_two[0]; double dL_2 = 1/fitted_parameters_species_two[1]; double dP_2 = 1/fitted_parameters_species_two[2];
  double muE0_2 = fitted_parameters_species_two[3]; double muL0_2 = fitted_parameters_species_two[4]; double muP_2 = fitted_parameters_species_two[5];
  double muM_2 = fitted_parameters_species_two[6]; double lambda_2 = fitted_parameters_species_two[7]; double tau_2 = fitted_parameters_species_two[8];
  double beta_2 = fitted_parameters_species_two[9];

  // Species Three
  double dE_3 = 1/fitted_parameters_species_three[0]; double dL_3 = 1/fitted_parameters_species_three[1]; double dP_3 = 1/fitted_parameters_species_three[2];
  double muE0_3 = fitted_parameters_species_three[3]; double muL0_3 = fitted_parameters_species_three[4]; double muP_3 = fitted_parameters_species_three[5];
  double muM_3 = fitted_parameters_species_three[6]; double lambda_3 = fitted_parameters_species_three[7]; double tau_3 = fitted_parameters_species_three[8];
  double beta_3 = fitted_parameters_species_three[9];

  // Setting the Static Parameters
  double dt = static_parameters[0];
  int dd_pow = static_parameters[1];

  // Setting the Initial State Variables
  std::vector<double> rF = rainfall;

  // Species One
  double E_1 = fitted_parameters_species_one[12]; double L_1 = fitted_parameters_species_one[13];
  double P_1 = fitted_parameters_species_one[14]; double M_1 = fitted_parameters_species_one[15];
  double Be_1 = 0; double Bl_1 = 0; double Bp_1 = 0; double Bm_1 = 0;

  // Species Two
  double E_2 = fitted_parameters_species_two[12]; double L_2 = fitted_parameters_species_two[13];
  double P_2 = fitted_parameters_species_two[14]; double M_2 = fitted_parameters_species_two[15];
  double Be_2 = 0; double Bl_2 = 0; double Bp_2 = 0; double Bm_2 = 0;

  // Species Three
  double E_3 = fitted_parameters_species_three[12]; double L_3 = fitted_parameters_species_three[13];
  double P_3 = fitted_parameters_species_three[14]; double M_3 = fitted_parameters_species_three[15];
  double Be_3 = 0; double Bl_3 = 0; double Bp_3 = 0; double Bm_3 = 0;

  // Declaring Additional Parameters

  // Species One
  double K_1_rainfall; // rainfall driven component of the ecological carrying capacity of the environment
  double K_1_static = species_one_static_K; // rainfall invaraiant component of the ecological carrying capacity of the environment
  double muE_1; // mortality rate for Early larvae including density dependence
  double muL_1; // mortality rate for Late larvae including density dependence
  int tau_with_dt_1 = tau_1 / dt; // the number of timesteps of past rainfall that contribute to K, taking into account the timestep
  double rFsum_1; // rainfall summed over the tau days
  double rFsum_fudge_1; // rainfall summed over the tau days
  double mRan_1; // how many events out of the total pupal events (Bp) to assign as development into mosquitoes

  // Species Two
  double K_2_rainfall; // rainfall driven component of the ecological carrying capacity of the environment
  double K_2_static = species_two_static_K; // rainfall invaraiant component of the ecological carrying capacity of the environment
  double muE_2; // mortality rate for Early larvae including density dependence
  double muL_2; // mortality rate for Late larvae including density dependence
  int tau_with_dt_2 = tau_2 / dt; // the number of timesteps of past rainfall that contribute to K, taking into account the timestep
  double rFsum_2; // rainfall summed over the tau days
  double rFsum_fudge_2; // rainfall summed over the tau days
  double mRan_2; // how many events out of the total pupal events (Bp) to assign as development into mosquitoes

  // Species Three
  double K_3_rainfall; // rainfall driven component of the ecological carrying capacity of the environment
  double K_3_static = species_three_static_K; // rainfall invaraiant component of the ecological carrying capacity of the environment
  double muE_3; // mortality rate for Early larvae including density dependence
  double muL_3; // mortality rate for Late larvae including density dependence
  int tau_with_dt_3 = tau_3 / dt; // the number of timesteps of past rainfall that contribute to K, taking into account the timestep
  double rFsum_3; // rainfall summed over the tau days
  double rFsum_fudge_3; // rainfall summed over the tau days
  double mRan_3; // how many events out of the total pupal events (Bp) to assign as development into mosquitoes

  // Vectors to store the output of the model at each timestep

  // Species One
  Rcpp::NumericVector E_output_1(end_time - start_time); // why can't I use vector<double> without stating size?
  Rcpp::NumericVector L_output_1(end_time - start_time); // also when I was initialising these without the size, it was crashing automatically
  Rcpp::NumericVector P_output_1(end_time - start_time);
  Rcpp::NumericVector M_output_1(end_time - start_time);
  Rcpp::NumericVector k_output_1(end);


  // Species Two
  Rcpp::NumericVector E_output_2(end_time - start_time);
  Rcpp::NumericVector L_output_2(end_time - start_time);
  Rcpp::NumericVector P_output_2(end_time - start_time);
  Rcpp::NumericVector M_output_2(end_time - start_time);
  Rcpp::NumericVector k_output_2(end);

  // Species Three
  Rcpp::NumericVector E_output_3(end_time - start_time);
  Rcpp::NumericVector L_output_3(end_time - start_time);
  Rcpp::NumericVector P_output_3(end_time - start_time);
  Rcpp::NumericVector M_output_3(end_time - start_time);
  Rcpp::NumericVector k_output_3(end);

  int i = 0;

  // Iterating the model through multiple timepoints
  while (t < end_time) {

    for (int j = 0; j < 3; j++) {

      if (j == 0) {

        // Specifying the value of K for all timepoints
        if (rainfall_relationship_species_one == "mean") {

          if (t <= tau_with_dt_1) {

            std::vector<double>::const_iterator first_fudge = rF.begin();
            std::vector<double>::const_iterator last_fudge = rF.begin() + t;
            std::vector<double> rFx_fudge(first_fudge, last_fudge);
            rFsum_fudge_1 = std::accumulate(rFx_fudge.begin(), rFx_fudge.end(), 0.0); // Don't forget to do 0.0!!! Otherwise accumulator will produce an int!!

            if (t == 0) {
              K_1_rainfall = rainfall[t];
            }
            else {
              K_1_rainfall = ((1.0 / t) * rFsum_fudge_1);
            }

            k_output_1[t] = K_1_rainfall;
          }

          else {

            std::vector<double>::const_iterator first = rF.begin() + t - tau_with_dt_1;
            std::vector<double>::const_iterator last = rF.begin() + t;
            std::vector<double> rFx(first, last);
            rFsum_1 = std::accumulate(rFx.begin(), rFx.end(), 0.0); // Don't forget to do 0.0!!! Otherwise accumulator will produce an int!!
            K_1_rainfall = (1.0 + ((1.0 / tau_with_dt_1) * rFsum_1));

            if (K_1_rainfall > max_K_1) {
              K_1_rainfall = max_K_1;
            }

            k_output_1[t] = K_1_rainfall;

          }
        }


        // Setting the density dependent function regulating larval mortality
        if (mortality_density_function_species_two == "linear") {

          if (M_2 < threshold_dd_1 & M_3 < threshold_dd_1) {

            muE_1 = muE0_1 * (1 + ((E_1 + L_1) / (K_1_rainfall + K_1_static)));
            muL_1 = muL0_1 * (1 + lambda_1 * (((E_1 + L_1) / (K_1_rainfall + K_1_static))));

          }

          else if (M_2 < threshold_dd_1 & M_3 >= threshold_dd_1) {

            muE_1 = muE0_1 * (1 + ((E_1 + L_1 + E_3 * competition_parameter_three_on_one + L_3 * competition_parameter_three_on_one) / (K_1_rainfall + K_1_static)));
            muL_1 = muL0_1 * (1 + lambda_1 * (((E_1 + L_1 + E_3 * competition_parameter_three_on_one + L_3 * competition_parameter_three_on_one) / (K_1_rainfall + K_1_static))));

          }

          else if (M_2 >= threshold_dd_1 & M_3 < threshold_dd_1) {

            muE_1 = muE0_1 * (1 + ((E_1 + L_1 + E_2 * competition_parameter_two_on_one + L_2 * competition_parameter_two_on_one) / (K_1_rainfall + K_1_static)));
            muL_1 = muL0_1 * (1 + lambda_1 * (((E_1 + L_1 + E_2 * competition_parameter_two_on_one + L_2 * competition_parameter_two_on_one) / (K_1_rainfall + K_1_static))));

          }

          else if (M_2 >= threshold_dd_1 & M_3 >= threshold_dd_1) {

            muE_1 = muE0_1 * (1 + ((E_1 + L_1 + E_3 * competition_parameter_three_on_one + L_3 * competition_parameter_three_on_one + E_2 * competition_parameter_two_on_one + L_2 * competition_parameter_two_on_one) / (K_1_rainfall + K_1_static)));
            muL_1 = muL0_1 * (1 + lambda_1 * (((E_1 + L_1 + E_3 * competition_parameter_three_on_one + L_3 * competition_parameter_three_on_one + E_2 * competition_parameter_two_on_one + L_2 * competition_parameter_two_on_one) / (K_1_rainfall + K_1_static))));

          }

        }

        // Specifying the total number of events occurring for each compartment at each timepoint

        // Be
        if ((dE_1 + muE_1)*dt < 1) {
          Be_1 = R::rbinom(E_1, (dE_1 + muE_1)*dt);
        }
        else {
          Be_1 = R::rbinom(E_1, 1);
        }

        // Bl
        if ((dL_1 + muL_1)*dt < 1) {
          Bl_1 = R::rbinom(L_1, (dL_1 + muL_1)*dt);
        }
        else {
          Bl_1 = R::rbinom(L_1, 1);
        }

        // Bp
        if ((dP_1 + muP_1)*dt < 1) {
          Bp_1 = R::rbinom(P_1, (dP_1 + muP_1)*dt);
        }
        else {
          Bp_1 = R::rbinom(P_1, 1);
        }

        // Bm (dealing with instances of 0 mosquitoes) - Mosquito Deaths
        if (M_1 >= 1) { // Aaron has >= 1, surely should be > 1??
          Bm_1 = R::rbinom(M_1, muM_1 * dt);
        }
        else {
          Bm_1 = 0;
        }

        // Updating the State Variables

        // E
        E_1 = round(E_1 - Be_1 + M_1 * beta_1 * dt);
        if (E_1 <= 0) {
          E_1 = 0;
        }

        // L
        if (Be_1 >= 1) {
          L_1 = round(L_1 - Bl_1 + R::rbinom(Be_1, (dE_1 / (muE_1 + dE_1))));
        }
        else {
          L_1 = round(L_1 - Bl_1);
        }

        // P
        if (Bl_1 >= 1) {
          P_1 = round(P_1 - Bp_1 + R::rbinom(Bl_1, (dL_1 / (muL_1 + dL_1))));
        }
        else {
          P_1 = round(P_1 - Bp_1);
        }

        // M
        if (Bp_1 > 1) {
          mRan_1 = R::rbinom(Bp_1, (dP_1 / (muP_1 + dP_1)));
        }
        else {
          mRan_1 = 0;
        }
        M_1 = round(M_1 + (0.5 * mRan_1) - Bm_1);
        if (M_1 < 1) {
          M_1 = 1;
        }

        E_output_1[i] = E_1;
        L_output_1[i] = L_1;
        P_output_1[i] = P_1;
        M_output_1[i] = M_1;

      }

      if (j == 1) {

        // Specifying the value of K for all timepoints
        if (rainfall_relationship_species_two == "mean") {

          if (t <= tau_with_dt_2) {

            std::vector<double>::const_iterator first_fudge = rF.begin();
            std::vector<double>::const_iterator last_fudge = rF.begin() + t;
            std::vector<double> rFx_fudge(first_fudge, last_fudge);
            rFsum_fudge_2 = std::accumulate(rFx_fudge.begin(), rFx_fudge.end(), 0.0); // Don't forget to do 0.0!!! Otherwise accumulator will produce an int!!

            if (t == 0) {
              K_2_rainfall = rainfall[t];
            }
            else {
              K_2_rainfall = ((1.0 / t) * rFsum_fudge_2);
            }

            k_output_2[t] = K_2_rainfall;
          }

          else {

            std::vector<double>::const_iterator first = rF.begin() + t - tau_with_dt_2;
            std::vector<double>::const_iterator last = rF.begin() + t;
            std::vector<double> rFx(first, last);
            rFsum_2 = std::accumulate(rFx.begin(), rFx.end(), 0.0); // Don't forget to do 0.0!!! Otherwise accumulator will produce an int!!
            K_2_rainfall = (1.0 + ((1.0 / tau_with_dt_2) * rFsum_2));

            if (K_2_rainfall > max_K_2) {
              K_2_rainfall = max_K_2;
            }

            k_output_2[t] = K_2_rainfall;

          }
        }


        // Setting the density dependent function regulating larval mortality
        if (mortality_density_function_species_two == "linear") {

          if (M_1 < threshold_dd_2 & M_3 < threshold_dd_2) {

            muE_2 = muE0_2 * (1 + ((E_2 + L_2) / (K_2_rainfall + K_2_static)));
            muL_2 = muL0_2 * (1 + lambda_2 * (((E_2 + L_2) / (K_2_rainfall + K_2_static))));

          }

          else if (M_1 < threshold_dd_2 & M_3 >= threshold_dd_2) {

            muE_2 = muE0_2 * (1 + ((E_2 + L_2 + E_3 * competition_parameter_three_on_two + L_3 * competition_parameter_three_on_two) / (K_2_rainfall + K_2_static)));
            muL_2 = muL0_2 * (1 + lambda_2 * (((E_2 + L_2 + E_3 * competition_parameter_three_on_two + L_3 * competition_parameter_three_on_two) / (K_2_rainfall + K_2_static))));

          }

          else if (M_1 >= threshold_dd_2 & M_3 < threshold_dd_2) {

            muE_2 = muE0_2 * (1 + ((E_2 + L_2 + E_1 * competition_parameter_one_on_two + L_1 * competition_parameter_one_on_two) / (K_2_rainfall + K_2_static)));
            muL_2 = muL0_2 * (1 + lambda_2 * (((E_2 + L_2 + E_1 * competition_parameter_one_on_two + L_1 * competition_parameter_one_on_two) / (K_2_rainfall + K_2_static))));

          }

          else if (M_1 >= threshold_dd_2 & M_3 >= threshold_dd_2) {

            muE_2 = muE0_2 * (1 + ((E_2 + L_2 + E_1 * competition_parameter_one_on_two + L_1 * competition_parameter_one_on_two + E_3 * competition_parameter_three_on_two + L_3 * competition_parameter_three_on_two) / (K_2_rainfall + K_2_static)));
            muL_2 = muL0_2 * (1 + lambda_2 * (((E_2 + L_2 + E_1 * competition_parameter_one_on_two + L_1 * competition_parameter_one_on_two + E_3 * competition_parameter_three_on_two + L_3 * competition_parameter_three_on_two) / (K_2_rainfall + K_2_static))));

          }

        }


        // Specifying the total number of events occurring for each compartment at each timepoint

        // Be
        if ((dE_2 + muE_2)*dt < 1) {
          Be_2 = R::rbinom(E_2, (dE_2 + muE_2)*dt);
        }
        else {
          Be_2 = R::rbinom(E_2, 1);
        }

        // Bl
        if ((dL_2 + muL_2)*dt < 1) {
          Bl_2 = R::rbinom(L_2, (dL_2 + muL_2)*dt);
        }
        else {
          Bl_2 = R::rbinom(L_2, 1);
        }

        // Bp
        if ((dP_2 + muP_2)*dt < 1) {
          Bp_2 = R::rbinom(P_2, (dP_2 + muP_2)*dt);
        }
        else {
          Bp_2 = R::rbinom(P_2, 1);
        }

        // Bm (dealing with instances of 0 mosquitoes) - Mosquito Deaths
        if (M_2 >= 1) { // Aaron has >= 1, surely should be > 1??
          Bm_2 = R::rbinom(M_2, muM_2 * dt);
        }
        else {
          Bm_2 = 0;
        }

        // Updating the State Variables

        // E
        E_2 = round(E_2 - Be_2 + M_2 * beta_2 * dt);
        if (E_2 <= 0) {
          E_2 = 0;
        }

        // L
        if (Be_2 >= 1) {
          L_2 = round(L_2 - Bl_2 + R::rbinom(Be_2, (dE_2 / (muE_2 + dE_2))));
        }
        else {
          L_2 = round(L_2 - Bl_2);
        }

        // P
        if (Bl_2 >= 1) {
          P_2 = round(P_2 - Bp_2 + R::rbinom(Bl_2, (dL_2 / (muL_2 + dL_2))));
        }
        else {
          P_2 = round(P_2 - Bp_2);
        }

        // M
        if (Bp_2 > 1) {
          mRan_2 = R::rbinom(Bp_2, (dP_2 / (muP_2 + dP_2)));
        }
        else {
          mRan_2 = 0;
        }
        M_2 = round(M_2 + (0.5 * mRan_2) - Bm_2);
        if (M_2 < 1) {
          M_2 = 1;
        }

        E_output_2[i] = E_2;
        L_output_2[i] = L_2;
        P_output_2[i] = P_2;
        M_output_2[i] = M_2;

      }

      if (j == 2) {

        // Specifying the value of K for all timepoints
        if (rainfall_relationship_species_three == "mean") {

          if (t <= tau_with_dt_3) {

            std::vector<double>::const_iterator first_fudge = rF.begin();
            std::vector<double>::const_iterator last_fudge = rF.begin() + t;
            std::vector<double> rFx_fudge(first_fudge, last_fudge);
            rFsum_fudge_3 = std::accumulate(rFx_fudge.begin(), rFx_fudge.end(), 0.0); // Don't forget to do 0.0!!! Otherwise accumulator will produce an int!!

            if (t == 0) {
              K_3_rainfall = rainfall[t];
            }
            else {
              K_3_rainfall = ((1.0 / t) * rFsum_fudge_3);
            }

            k_output_3[t] = K_3_rainfall;
          }

          else {

            std::vector<double>::const_iterator first = rF.begin() + t - tau_with_dt_3;
            std::vector<double>::const_iterator last = rF.begin() + t;
            std::vector<double> rFx(first, last);
            rFsum_3 = std::accumulate(rFx.begin(), rFx.end(), 0.0); // Don't forget to do 0.0!!! Otherwise accumulator will produce an int!!
            K_3_rainfall = (1.0 + ((1.0 / tau_with_dt_3) * rFsum_3));

            if (K_3_rainfall > max_K_3) {
              K_3_rainfall = max_K_3;
            }

            k_output_3[t] = K_3_rainfall;

          }
        }


        // Setting the density dependent function regulating larval mortality
        if (mortality_density_function_species_three == "linear") {

          if (M_1 < threshold_dd_3 & M_2 < threshold_dd_3) {

            muE_3 = muE0_3 * (1 + ((E_3 + L_3) / (K_3_rainfall + K_3_static)));
            muL_3 = muL0_3 * (1 + lambda_3 * (((E_3 + L_3) / (K_3_rainfall + K_3_static))));

          }

          else if (M_1 < threshold_dd_3 & M_2 >= threshold_dd_3) {

            muE_3 = muE0_3 * (1 + ((E_3 + L_3 + E_2 * competition_parameter_two_on_three + L_2 * competition_parameter_two_on_three) / (K_3_rainfall + K_3_static)));
            muL_3 = muL0_3 * (1 + lambda_3 * (((E_3 + L_3 + E_2 * competition_parameter_two_on_three + L_2 * competition_parameter_two_on_three) / (K_3_rainfall + K_3_static))));

          }

          else if (M_1 >= threshold_dd_3 & M_2 < threshold_dd_3) {

            muE_3 = muE0_3 * (1 + ((E_3 + L_3 + E_1 * competition_parameter_one_on_three + L_1 * competition_parameter_one_on_three) / (K_3_rainfall + K_3_static)));
            muL_3 = muL0_3 * (1 + lambda_3 * (((E_3 + L_3 + E_1 * competition_parameter_one_on_three + L_1 * competition_parameter_one_on_three) / (K_3_rainfall + K_3_static))));

          }

          else if (M_1 >= threshold_dd_3 & M_2 >= threshold_dd_3) {

            muE_3 = muE0_3 * (1 + ((E_3 + L_3 + E_1 * competition_parameter_one_on_three + L_1 * competition_parameter_one_on_three + E_2 * competition_parameter_two_on_three + L_2 * competition_parameter_two_on_three) / (K_3_rainfall + K_3_static)));
            muL_3 = muL0_3 * (1 + lambda_3 * (((E_3 + L_3 + E_1 * competition_parameter_one_on_three + L_1 * competition_parameter_one_on_three + E_2 * competition_parameter_two_on_three + L_2 * competition_parameter_two_on_three) / (K_3_rainfall + K_3_static))));

          }

        }

        // Specifying the total number of events occurring for each compartment at each timepoint

        // Be
        if ((dE_3 + muE_3)*dt < 1) {
          Be_3 = R::rbinom(E_3, (dE_3 + muE_3)*dt);
        }
        else {
          Be_3 = R::rbinom(E_3, 1);
        }

        // Bl
        if ((dL_3 + muL_3)*dt < 1) {
          Bl_3 = R::rbinom(L_3, (dL_3 + muL_3)*dt);
        }
        else {
          Bl_3 = R::rbinom(L_3, 1);
        }

        // Bp
        if ((dP_3 + muP_3)*dt < 1) {
          Bp_3 = R::rbinom(P_3, (dP_3 + muP_3)*dt);
        }
        else {
          Bp_3 = R::rbinom(P_3, 1);
        }

        // Bm (dealing with instances of 0 mosquitoes) - Mosquito Deaths
        if (M_3 >= 1) { // Aaron has >= 1, surely should be > 1??
          Bm_3 = R::rbinom(M_3, muM_3 * dt);
        }
        else {
          Bm_3 = 0;
        }

        // Updating the State Variables

        // E
        E_3 = round(E_3 - Be_3 + M_3 * beta_3 * dt);
        if (E_3 <= 0) {
          E_3 = 0;
        }

        // L
        if (Be_3 >= 1) {
          L_3 = round(L_3 - Bl_3 + R::rbinom(Be_3, (dE_3 / (muE_3 + dE_3))));
        }
        else {
          L_3 = round(L_3 - Bl_3);
        }

        // P
        if (Bl_3 >= 1) {
          P_3 = round(P_3 - Bp_3 + R::rbinom(Bl_3, (dL_3 / (muL_3 + dL_3))));
        }
        else {
          P_3 = round(P_3 - Bp_3);
        }

        // M
        if (Bp_3 > 1) {
          mRan_3 = R::rbinom(Bp_3, (dP_3 / (muP_3 + dP_3)));
        }
        else {
          mRan_3 = 0;
        }
        M_3 = round(M_3 + (0.5 * mRan_3) - Bm_3);
        if (M_3 < 1) {
          M_3 = 1;
        }

        E_output_3[i] = E_3;
        L_output_3[i] = L_3;
        P_output_3[i] = P_3;
        M_output_3[i] = M_3;

      }

    }

    // Stepping forward another timestep
    t = t + 1;
    i = i + 1;

  }

  return Rcpp::List::create(Rcpp::Named("Species_1_E_Output") = E_output_1,
                            Rcpp::Named("Species_1_L_Output") = L_output_1,
                            Rcpp::Named("Species_1_P_Output") = P_output_1,
                            Rcpp::Named("Species_1_M_Output") = M_output_1,
                            Rcpp::Named("Species_1_K") = k_output_1,
                            Rcpp::Named("Species_2_E_Output") = E_output_2,
                            Rcpp::Named("Species_2_L_Output") = L_output_2,
                            Rcpp::Named("Species_2_P_Output") = P_output_2,
                            Rcpp::Named("Species_2_M_Output") = M_output_2,
                            Rcpp::Named("Species_2_K") = k_output_2,
                            Rcpp::Named("Species_3_E_Output") = E_output_3,
                            Rcpp::Named("Species_3_L_Output") = L_output_3,
                            Rcpp::Named("Species_3_P_Output") = P_output_3,
                            Rcpp::Named("Species_3_M_Output") = M_output_3,
                            Rcpp::Named("Species_3_K") = k_output_3);
}
