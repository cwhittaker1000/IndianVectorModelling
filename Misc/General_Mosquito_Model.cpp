/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
///                                                                                                   ///
///  Stochastic Model of Vector Population Dynamics                                                   ///
///                                                                                                   ///
///  Charlie Whittaker                                                                                ///
///  Imperial College London                                                                          ///
///  charles.whittaker16@imperial.ac.uk                                                               ///
///                                                                                                   ///
///  The below model represents a stochastic, compartmental model of                                  ///
///  of mosquito population dynamics based on that developed by Michael                               ///
///  White and detailed in White et al., Parasites & Vectors, 2011 and                                ///
///  available here: https://parasitesandvectors.biomedcentral.com/articles/10.1186/1756-3305-4-153   ///
///                                                                                                   ///
///  It extends it in a number of key ways. Specifically it introduces representations of different   ///
///  hydrological phenomena in order to account for the diversity of different population dynamics    ///
///  observed for mosquitoes in Indian settings. These include:                                       ///
///                                                                                                   ///
///       1) A Static, Permanent Component to the Carrying Capacity- meant to represent               ///
///          standing bodies of water or similar aquatic habitats that are available and              ///
///          suitable for breeding year round, in a manner indepenent of incipient rainfall.          ///
///       2) Washout and Dryout- meant to represent the streams and rivers where high levels          ///
///          of rainfall (and associated river flow) might render them unsuitable for breeding        ///
///          (Washout) but otherwise they are suitable. Additionally can incorporate reductions       ///
///          in the suitability of these environments in the absence of rain (Dryout).                ///
///                                                                                                   ///
///  More information on this model formulation available here:                                       ///
///                                                                                                   ///
///  The code below (specifying the model run function general_mosquito_population_model) is          ///
///  structured as follows:                                                                           ///
///                                                                                                   ///
///       1) PARAMETER SETTING & INITIALISATION                                                       ///
///       2) (optional) PERFORM COMPUTATIONALLY EXPENSIVE EXPONENTIAL CARRYING CAPACITIES OUTSIDE     ///
///                     THE MAIN LOOP SO THAT IT ONLY NEEDS TO BE DONE ONCE. (ONLY REQUIRED IF THE    ///
///                     CARRYING CAPACITY RELATIONSHIP FOR K_RAIN IS SET TO BE EXPONENTIAL- CAN       ///
///                     BE DONE INSIDE THE MODEL OR OUTSIDE, IN THE PARTICLE FILTER).                 ///
///       3) (optional) CALCULATE CARRYING CAPACITY FOR TIMEPOINT PRECEDING THE TIMEPOINTS FOR WHICH  ///
///                     DATA IS AVAILABLE. DONE IN ORDER TO TRACK HOW RECENTLY PRIOR TO THE DATA      ///
///                     WASHOUT WAS LIKELY TO HAVE OCCURRED. (ONLY REQUIRED IF WASHOUT IS OPERATING   ///
///                     IN THE MODEL)                                                                 ///
///       4) CALCULATING THE VALUE OF THE CARRYING CAPACITY (K) FOR THE CURRENT TIMEPOINT             ///
///       5) MODEL RUNNING BETWEEN CURRENT TIMEPOINT AND THE NEXT                                     ///
///                                                                                                   ///
///                                                                                                   ///
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

// Specifying the Includes and Depends Required
#include "Mosquito_Population_Model.hpp"
#include "Hill_Function.hpp"
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]


//' @export
// [[Rcpp::export]]
Rcpp::List general_mosquito_population_model(int start_time, int end, // Start and Endtimes
                                             Rcpp::NumericVector fitted_parameters, Rcpp::NumericVector static_parameters, // Model Parameters
                                             std::vector<double> rainfall, // Rainfall
                                             Rcpp::String mortality_density_function,  // Vector Mortality Density Dependence Formulation
                                             Rcpp::String rainfall_relationship, // Weighting of Recent and Past Rainfall (Mean, Linear, Exponential)
                                             Rcpp::String rainfall_effect, // The Effect of Rainfall on Carrying Capacity (Raw or Hill)
                                             Rcpp::String decline_type, // How the Carrying Capacity Declines During Absence of Rainfall (Exponential or Hill)
                                             Rcpp::String calc_inside_mosquito_model, // If rainfall_relationship = Exponential, Whether to Calculate Inside Model or Not (Outside = Faster)
                                             std::vector<double> input_Exponential_Weighting_Factors_Static, // Related to calc_inside_mosquito_model
                                             std::vector<double> input_Exponential_Weighting_Factors_Rainfall, // Related to calc_inside_mosquito_model
                                             std::vector<double> input_Exponential_Normalisation_Factors_Static, // Related to calc_inside_mosquito_model
                                             std::vector<double> input_Exponential_Normalisation_Factors_Rainfall, // Related to calc_inside_mosquito_model
                                             int full_output) { // Specifies whether to output lots of things or just model compartment outputs


  ///////////////////////////////////////////////////////////////////
  //                                                               //
  //  1.  Parameter Setting and Initialisation                     //
  //                                                               //
  ///////////////////////////////////////////////////////////////////


  // Start and Endtime for Model Running
  int t = start_time;
  int end_time = end;

  // Initial State Variables
  double E = fitted_parameters[25]; double L = fitted_parameters[26];  double P = fitted_parameters[27]; double M = fitted_parameters[28];
  double Be = 0; double Bl = 0; double Bp = 0; double Bm = 0;

  // Constant Parameters
  double dt = static_parameters[0];   int dd_pow = static_parameters[1]; double offset = static_parameters[2]; std::vector<double> rF = rainfall;

  // Core Parameters (Present Across All Three Model Modalities )
  double dE = 1/fitted_parameters[0]; double dL = 1/fitted_parameters[1]; double dP = 1/fitted_parameters[2]; double muE0 = fitted_parameters[3];
  double muL0 = fitted_parameters[4]; double muP = fitted_parameters[5]; double muM = fitted_parameters[6];
  double lambda = fitted_parameters[7]; double beta = fitted_parameters[8]; // fitted_parameters[9], [10] and [11] are overdisp, pop_frac and z and these don't directly feature in the model at this stage.

  // Parameters Governing K_Rain(t) (Present in A.annularis & A.culicifacies Modalities Only)
  double tau_rain = fitted_parameters[12]; // Involved in the weighting of past and present rainfall, see White et al for full description
  double scaling_factor_rainfall = fitted_parameters[13]; // Scaling factor using to adjust the absolute magnitude of the carrying capacity calculated
  double K_Max_Hill_Rainfall = fitted_parameters[14]; // Parameters for the Hill Function that can be used to calculate the carrying capacity instead of raw rainfall total
  double Hill_Rainfall_1 = fitted_parameters[15]; // Hill Function Parameter
  double Hill_Rainfall_2 = fitted_parameters[16]; // Hill Function Parameter

  // Parameters Governing K_Static(t) (Present in A.annularis & A.fluviatilis Modalities Only)
  double tau_static = fitted_parameters[17]; // Involved in the weighting of past and present rainfall, see White et al for full description
  double scaling_factor_static = fitted_parameters[18]; // Scaling factor using to adjust the absolute magnitude of the carrying capacity calculated
  double K_Max_Static = fitted_parameters[19] * scaling_factor_static; // Max value of K_Static
  double Washout_Occurrence = fitted_parameters[20]; // Whether Washout is Occurring or Not
  double Washout_Threshold = fitted_parameters[21]; // Threshold of incipient rainfall at which washout starts occurring
  double washout_exp_scaling_factor = fitted_parameters[22]; // Parameters for the Hill Function that can be used to calculate the carrying capacity instead of raw rainfall total
  double washout_hill_one = fitted_parameters[23]; // Hill Function Parameter
  double washout_hill_two = fitted_parameters[24]; // Hill Function Parameter
  double marker = 0; // Dummy variable involved in tracking when washout last occured. 1 if Washout just occured, 0 otherwise (check).
  double marker_stop; // Dummy variable involed in tracking when washout last occured.

  // Initialising Variables Required in the Model
  double K_Static; // permanent component of the environmental carrying capacity - can either be constant or temporally variable
  double K_Rain; // rainfall responsive component of the environmental carrying capacity - temporally variable
  double muE; // mortality rate for Early larvae including density dependence
  double muL; // mortality rate for Late larvae including density dependence
  int tau_with_dt_rain = tau_rain / dt; // the number of timesteps of past rainfall that contribute to K_Rain, taking into account the timestep
  int tau_with_dt_static = tau_static / dt; // the number of timesteps of past rainfall that contribute to K_Static washout calculations, taking into account the timestep
  double rFsum_rain; // rainfall summed over the tau_rain days
  double rFsum_static; // rainfall summed over the tau_static days
  double rFsum_static_prior; // rainfall summed over the previous tau_static days, but for the rainfall timepoints before data start
  double mRan; // how many events out of the total pupal events (Bp) to assign as development into mosquitoes

  // Vectors to store the output of the model at each timestep
  Rcpp::NumericVector E_output(end_time - start_time); Rcpp::NumericVector L_output(end_time - start_time);
  Rcpp::NumericVector P_output(end_time - start_time); Rcpp::NumericVector M_output(end_time - start_time);
  Rcpp::NumericVector k_output(end_time - start_time); Rcpp::NumericVector k_rain_output(end);
  Rcpp::NumericVector k_static_output(end); Rcpp::NumericVector k_total_output(end);
  Rcpp::NumericVector average_rainfall_K_Static(end);

  // Checking some of the exponential calculations to make sure they're right
  Rcpp::NumericVector prior_K_Static_values(start_time + offset - tau_with_dt_static + 1);
  int i = 0;
  int prior_counter = 0;

  // Rcpp::Rcout << "Step 1: Parameter Initialisation COMPLETE" << std::endl;

  ///////////////////////////////////////////////////////////////////
  //                                                               //
  //  2.  (optional) Perform Exponential Carrying Capacity         //
  //                 Equations Outside Loop                        //
  //                                                               //
  ///////////////////////////////////////////////////////////////////

  // Using the exponential relationship to calculate the contribution of rainfall at different timepoints is very computationally expensive.
  // Particularly when implemented as part of the Particle Filter (requiring running many model simulation simultaneously) and in this
  // model context (as we consider time prior to when we have data, increasing the number of timepoints the calculations have to be done for).
  //
  // This code attempts to ease this computational burden by calculating the quantities that comprise the exponential weighting function,
  // specifically the normalisation factor and the weighting factor, outside of the main loop. This avoids having to calculate the same values
  // multiple times, reducing the computational intensity of running the model with this relationship.
  //
  // This section produces two outputs:
  //
  //        1. Exponential_Weighting_Factor - A vector of weighting factors where the first value in the vector is the weighting factor
  //                                          for the nearest day to the current day (i.e. the current day). The second value is the
  //                                          weighting factor that should be used to weight the rainfall that occurred the day before
  //                                          the current day etc.
  //        2. Exponential_Normalisation_Factor - A vector of normalisation factors where the first value in the vector is the normalisation
  //                                              factor to be used at Day 1 (where we only consider 1 day's worth of rainfall), the second
  //                                              value the factor to be used at Day 2 (where we consider 2 day's worth of rainfall)etc.
  //
  // Remember:              K  =  1/(tau * (1 - exp(-t/tau))) * integral (exp(-(t-t')/tau)) * raint' dt'
  //                              --- Normalisation Factor---            --Weighting Factor--
  //
  // Note: tau_static and tau_rain can be different so this has to be done twice, for both components of the carrying capacity.
  //
  // See White et al., for more details

  // Initialising temporary/dummy variables and vectors to store calculations within the loop
  double temp_static_weighting_factor;
  double temp_rain_weighting_factor;
  double temp_static_normalisation_factor;
  double temp_rain_normalisation_factor;
  std::vector<double> temp_vector_Exponential_Weighting_Factors_Static;  // These are not superfluous- produces different results if I'm adding to the
  std::vector<double> temp_vector_Exponential_Weighting_Factors_Rainfall; // vector directly as opposed to creating a temporary one and then copying the
  std::vector<double> temp_vector_Exponential_Normalisation_Factors_Static; // results over after the temporary one has been filled.
  std::vector<double> temp_vector_Exponential_Normalisation_Factors_Rainfall;  // ASK RICH WHY THIS IS!!!!
  std::vector<double> Exponential_Weighting_Factors_Static;
  std::vector<double> Exponential_Weighting_Factors_Rainfall;
  std::vector<double> Exponential_Normalisation_Factors_Static;
  std::vector<double> Exponential_Normalisation_Factors_Rainfall;

  // Only gets run if calc_inside_mosquito_model is set to Yes. If no, these calculations get done in Particle Filter, so that it can be done once for all Particles,
  // not individually for each Particle.
  if (rainfall_relationship == "exponential" & calc_inside_mosquito_model == "Yes") {

    // Timeperiod over which the factors must be calculated
    int overall_time_length = offset + end_time;

    // Generating the Vector of Exponential Weighting Factors
    for (int i = 0; i <= overall_time_length; i++) {

      temp_static_weighting_factor = -((offset + end_time) - i)/tau_with_dt_static;;
      temp_rain_weighting_factor = -((offset + end_time) - i)/tau_with_dt_rain;

      temp_vector_Exponential_Weighting_Factors_Static.insert(temp_vector_Exponential_Weighting_Factors_Static.begin(), exp(temp_static_weighting_factor));
      temp_vector_Exponential_Weighting_Factors_Rainfall.insert(temp_vector_Exponential_Weighting_Factors_Rainfall.begin(), exp(temp_rain_weighting_factor));

    }

    Exponential_Weighting_Factors_Static = temp_vector_Exponential_Weighting_Factors_Static;
    Exponential_Weighting_Factors_Rainfall = temp_vector_Exponential_Weighting_Factors_Rainfall;

    // Generating the Vector of Exponential Normalisation Factors
    for (int x = 0; x <= overall_time_length; x++) {

      temp_static_normalisation_factor = (1.0 / (tau_with_dt_static * (1 - exp(- (x + 1) / tau_with_dt_static))));
      temp_rain_normalisation_factor = (1.0 / (tau_with_dt_rain * (1 - exp(-(x + 1) / tau_with_dt_rain))));

      temp_vector_Exponential_Normalisation_Factors_Static.emplace_back(temp_static_normalisation_factor);
      temp_vector_Exponential_Normalisation_Factors_Rainfall.emplace_back(temp_rain_normalisation_factor);

    }

    Exponential_Normalisation_Factors_Static = temp_vector_Exponential_Normalisation_Factors_Static;
    Exponential_Normalisation_Factors_Rainfall = temp_vector_Exponential_Normalisation_Factors_Rainfall;

  }

  else if (rainfall_relationship == "exponential" & calc_inside_mosquito_model == "No"){

    Exponential_Weighting_Factors_Static = input_Exponential_Weighting_Factors_Static;
    Exponential_Weighting_Factors_Rainfall = input_Exponential_Weighting_Factors_Rainfall;
    Exponential_Normalisation_Factors_Static = input_Exponential_Normalisation_Factors_Static;
    Exponential_Normalisation_Factors_Rainfall = input_Exponential_Normalisation_Factors_Rainfall;

  }

  // Rcpp::Rcout << "Step 2: Exponential Weightings COMPLETE" << std::endl;

  ///////////////////////////////////////////////////////////////////
  //                                                               //
  //  3.  (optional) Calculate Weighted Rainfall for Timepoints    //
  //                 Preceding Data Timepoints. Done to Assess     //
  //                 Whether and When Washout Occurred.            //
  //                                                               //
  ///////////////////////////////////////////////////////////////////

  // Incorporating Washout (and Dryout) into the mosquito population model allows us to model a wider array of hydrological
  // breeding habitats. However, it poses an issue. Specifically, that we need to know where within the Washout/Dryout cycle
  // we're are situated when our data begins. To get round this, the model incorporates rainfall prior to when our data begin
  // (encompassed by the "Offset" term). During this "Offset" period, the weighted/normalised rainfall is calculated, and
  // whether Washout/Dryout is occurring is also determined. This allows us to estimate the carrying capacity at the point where
  // our data begins more accurately, and allowing the model to better capture the observed dynamics.
  //
  // This section of the code therefore takes rainfall data for the period spanning tau_with_dt_static --> offset (both of which
  // are timepoints) and calculates how the carrying capacity should have behaved prior to our observed data given a particular
  // value of WashoutThreshold and washout_exp_scaling_factor (which are two parameters we're estimating).
  //
  // It produces three outputs:
  //
  //        1. marker - A dummy variable that tracks whether washout has just occured. If Washout is occurring, it gets set to 0.
  //                    Else, the first time washout stops, it gets set to 1.
  //        2. marker_stop - A variable that tracks the last timepoint when washout occurred. Required to calculate the value of
  //                         the carrying capacity as it declines due to dryout.
  //        3. K_Static - The value of K_Static at the timepoint immediately preceding the time we have data starting from;
  //                      taking into account the cycles of Washout and Dryout that have preceded our data.
  //

  // Calculates Washout occurrence, Time Since Washout Occurrence and Carrying Capacity value for all timepoints prior to data
  // beginning (note we skip the first tau_with_dt_static timepoints for convenience)
  if (Washout_Occurrence == 1) {

    for (int p = tau_with_dt_static; p < (offset + start_time); p++) {


      // Rainfall Weighting is the Mean of Rainfall Fallen in the Past tau_with_dt_static Days
      if (rainfall_relationship == "mean") {

        // For Each Timepoint, Calculates the Summed Rainfall Over the Past tau_with_dt_static Days
        std::vector<double>::const_iterator first_point_offset_period = rF.begin() + p - tau_with_dt_static;
        std::vector<double>::const_iterator last_point_offset_period = rF.begin() + p;
        std::vector<double> rFx_static_offset_period(first_point_offset_period, last_point_offset_period);
        rFsum_static_prior = std::accumulate(rFx_static_offset_period.begin(), rFx_static_offset_period.end(), 0.0);
        double rainfall_average_static_calc_prior = (1/tau_with_dt_static) * rFsum_static_prior;

        // Calculating Whether Washout Occurs
        if (rainfall_average_static_calc_prior > Washout_Threshold) { // Washout occurs
          marker = 0; // marker set to 0 the first time Washout occurs. Will continue to be set to 0 as long as washout keeps occurring.
        }
        else {
          if (marker == 0) { // First Time Washout Stops
            K_Static = K_Max_Static; //Carrying Capacity Is At Max
            marker_stop = p; // marker_stop set to the timepoint at which Washout stops.
            marker = 1; // marker set to 1 the first timepoint at which Washout stops.
          }
          else { // Following Washout cessation, the Carrying Capacity begins to decline
            if (decline_type == "exponential") {
              double calculation = -washout_exp_scaling_factor * (p - marker_stop);  // p - marker_stop repsents how much time has elapsed since washout ceased
              K_Static = K_Max_Static * exp(calculation);
            }
            else {
              K_Static = Hill_Function((p - marker_stop), K_Max_Static, washout_hill_one, washout_hill_two);
            }
          }
        }

        // Storing these Results in a Vector to Assess Whether they're Correct Via the Output
        prior_K_Static_values[prior_counter] = K_Static;
        prior_counter = prior_counter + 1;
      }

      // Rainfall Weighting is the Linearly Weighted Mean of Rainfall Fallen in the Past tau_with_dt_static Days
      else if (rainfall_relationship == "linear") {

        // For Each Timepoint, Calculates the Linearly Weighted Mean Rainfall Over the Past tau_with_dt_static Days
        double rFsum_statically_prior = 0.0;
        int counter_static_prior = 1;
        for (int r = (p - tau_with_dt_static); r <= p; r++) {
          rFsum_statically_prior = rFsum_statically_prior + (((p - tau_with_dt_static + counter_static_prior) - p + tau_with_dt_static) * rF[r]);
          counter_static_prior = counter_static_prior + 1;
        }
        double rainfall_average_static_calc_prior = ((2.0 / pow(tau_with_dt_static, 2)) * rFsum_statically_prior); // unsure if this is what I want, double check

        // Calculating Whether Washout Occurs
        if (rainfall_average_static_calc_prior > Washout_Threshold) { // Washout occurs
          marker = 0; // marker set to 0 the first time Washout occurs. Will continue to be set to 0 as long as washout keeps occurring.
        }
        else {
          if (marker == 0) { // First Time Washout Stops
            K_Static = K_Max_Static; //Carrying Capacity Is At Max
            marker_stop = p; // marker_stop set to the timepoint at which Washout stops.
            marker = 1; // marker set to 1 the first timepoint at which Washout stops.
          }
          else { // Then Carrying Capacity Begins to Decline
            if (decline_type == "exponential") {
              double calculation = -washout_exp_scaling_factor * (p - marker_stop);  // p - marker_stop repsents how much time has elapsed since washout ceased
              K_Static = K_Max_Static * exp(calculation);
            }
            else {
              K_Static = Hill_Function((p - marker_stop), K_Max_Static, washout_hill_one, washout_hill_two);
            }
          }
        }
        prior_K_Static_values[prior_counter] = K_Static;
        prior_counter = prior_counter + 1;
      }
      else if (rainfall_relationship == "exponential") {

        // For Each Timepoint, Calculates the Exponentially Weighted Mean Rainfall Over the Past tau_with_dt_static Days
        double temp_tau_static_prior = tau_with_dt_static;
        double rFsum_statically_prior = 0.0;
        double calc;
        for (int r = 0; r <= p; r++) {
          calc = Exponential_Weighting_Factors_Static[r];
          rFsum_statically_prior = rFsum_statically_prior + calc * rF[p - r];
        }
        double rainfall_average_static_calc_prior = Exponential_Normalisation_Factors_Static[p] * rFsum_statically_prior;

        // Calculating Whether Washout Occurs
        if (rainfall_average_static_calc_prior > Washout_Threshold) { // Washout occurs
          marker = 0; // marker set to 0 the first time Washout occurs. Will continue to be set to 0 as long as washout keeps occurring.
        }
        else {
          if (marker == 0) { // First Time Washout Stops
            K_Static = K_Max_Static; //Carrying Capacity Is At Max
            marker_stop = p; // marker_stop set to the timepoint at which Washout stops.
            marker = 1; // marker set to 1 the first timepoint at which Washout stops.
          }
          else { // Then Carrying Capacity Begins to Decline
            if (decline_type == "exponential") {
              double calculation = -washout_exp_scaling_factor * (p - marker_stop);  // p - marker_stop repsents how much time has elapsed since washout ceased
              K_Static = K_Max_Static * exp(calculation);
            }
            else {
              K_Static = Hill_Function((p - marker_stop), K_Max_Static, washout_hill_one, washout_hill_two);
            }
          }
        }
        prior_K_Static_values[prior_counter] = K_Static;
        prior_counter = prior_counter + 1;
      }
    }
    //Rcpp::Rcout << "Test that Washout Occurrence = 1 Returns this Message" << std::endl;
  }

  int within_data_marker = 0; // dummy variable to span the offset period and the time-period for the data re washout. See Lines 419 - 428 for clarification
  // Rcpp::Rcout << "Step 3: Offset Timeperiod Running COMPLETE" << std::endl;

  ///////////////////////////////////////////////////////////////////
  //                                                               //
  //                      RUNNING THE MODEL                        //
  //                                                               //
  ///////////////////////////////////////////////////////////////////

  // Iterating the model through timepoints spanning the period for which we have data available
  while (t < end_time) {

    ///////////////////////////////////////////////////////////////////
    //                                                               //
    //  4.  Calculating the Value of the Carrying Capacity (K)       //
    //                                                               //
    ///////////////////////////////////////////////////////////////////

    // This section of code calculates the value of K_Rain and K_static for the timepoints which we have data available for.
    // It weights previous and past rainfall in one of three ways based on the relationship specified: Mean, Linear or Exponential.
    // This weighted rainfall is calculated and is then used either directly to calculate the carrying capacity (rainfall_effect = "raw")
    // or passed through a Hill Function (rainfall_effect = "Hill") to calculate it.

    // Mean Weighting Relationship

    if (rainfall_relationship == "mean") {

      // Calculating K_Rain
      std::vector<double>::const_iterator first_rain = rF.begin() + t + offset- tau_with_dt_rain;
      std::vector<double>::const_iterator last_rain = rF.begin() + t + offset;
      std::vector<double> rFx_rain(first_rain, last_rain);
      rFsum_rain = std::accumulate(rFx_rain.begin(), rFx_rain.end(), 0.0); // Don't forget to do 0.0!!! Otherwise accumulator will produce an int!!

      if (rainfall_effect == "raw") {
        K_Rain = scaling_factor_rainfall * (1.0 / tau_with_dt_rain) * rFsum_rain;

      }
      else if (rainfall_effect == "hill") {
        double mean_rainfall = (1.0 / tau_with_dt_rain) * rFsum_rain;
        double hill_output = Hill_Function(mean_rainfall, K_Max_Hill_Rainfall, Hill_Rainfall_1, Hill_Rainfall_2);
        K_Rain = scaling_factor_rainfall * hill_output;
      }

      // Calculating K_Static
      if (Washout_Occurrence == 0) {
        K_Static = K_Max_Static;
      }

      else if (Washout_Occurrence == 1) {

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
              if (within_data_marker == 0) { // if still 0, means the decline needs to carry on from the preceding offset period (marker will = 1 and within_data_marker will = 0 so this is run)
                double calculation = -washout_exp_scaling_factor * (t + offset - marker_stop);
                K_Static = K_Max_Static * exp(calculation);
              }
              else { // if within_data_marker = 1, it means a washout has occurred since the offset period, and hence can just use the marker_stop demarcating the most recent cessation of Washout
                double calculation = -washout_exp_scaling_factor * (t - marker_stop);
                K_Static = K_Max_Static * exp(calculation);
              }
            }
            else { // for when the decline type is Hill shaped
              if (within_data_marker == 0) {
                K_Static = Hill_Function((t + offset - marker_stop), K_Max_Static, washout_hill_one, washout_hill_two);
              }
              else {
                K_Static = Hill_Function((t - marker_stop), K_Max_Static, washout_hill_one, washout_hill_two);
              }
            }
          }
        }
      }
      k_static_output[t] = K_Static;
      k_rain_output[t] = K_Rain;
      k_total_output[t] = K_Static + K_Rain;
    }

    // Linear Weighting Relationship

    else if (rainfall_relationship == "linear") {

      // Calculating K_Rain
      double rFsum_rainfall = 0.0;
      int counter_rain = 1;
      for (int r = (t + offset - tau_with_dt_rain + 1); r <= t + offset; r++) {
        rFsum_rainfall = rFsum_rainfall + (((t - tau_with_dt_rain + counter_rain) - t + tau_with_dt_rain) * rF[r]);
        counter_rain = counter_rain + 1;
      }

      if (rainfall_effect == "raw") {
        K_Rain = scaling_factor_rainfall * (2.0 / pow(tau_with_dt_rain, 2)) * rFsum_rainfall;
      }
      else if (rainfall_effect == "hill") {
        double mean_rainfall = ((2.0 / pow(tau_with_dt_rain, 2)) * rFsum_rainfall);
        double hill_output = Hill_Function(mean_rainfall, K_Max_Hill_Rainfall, Hill_Rainfall_1, Hill_Rainfall_2);
        K_Rain = scaling_factor_rainfall * hill_output;
      }

      // Calculating K_Static
      if (Washout_Occurrence == 0) {
        K_Static = K_Max_Static;
      }

      else if (Washout_Occurrence == 1) {

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
              if (within_data_marker == 0) { // if still 0, means the decline needs to carry on from the preceding offset period (marker will = 1 and within_data_marker will = 0 so this is run)
                double calculation = -washout_exp_scaling_factor * (t + offset - marker_stop);
                K_Static = K_Max_Static * exp(calculation);
              }
              else { // if within_data_marker = 1, it means a washout has occurred since the offset period, and hence can just use the marker_stop demarcating the most recent cessation of Washout
                double calculation = -washout_exp_scaling_factor * (t - marker_stop);
                K_Static = K_Max_Static * exp(calculation);
              }
            }
            else { // for when the decline type is Hill shaped
              if (within_data_marker == 0) {
                K_Static = Hill_Function((t + offset - marker_stop), K_Max_Static, washout_hill_one, washout_hill_two);
              }
              else {
                K_Static = Hill_Function((t - marker_stop), K_Max_Static, washout_hill_one, washout_hill_two);
              }
            }
          }
        }
      }
      k_static_output[t] = K_Static;
      k_rain_output[t] = K_Rain;
      k_total_output[t] = K_Static + K_Rain;
    }

    // Exponential  Weighting Relationship

    else if (rainfall_relationship == "exponential") {

      // Calculating K_Rain
      double temp_tau_rain = tau_with_dt_rain;
      double rFsum_rainfall = 0.0;
      double calc;
      for (int r = 0; r <= t + offset; r++) { // need to decide whether to start at beginning of rain or t = 0 ??
        calc = Exponential_Weighting_Factors_Rainfall[r];
        rFsum_rainfall = rFsum_rainfall + calc * rF[(t + offset) - r];
      }

      if (rainfall_effect == "raw") {
        K_Rain = scaling_factor_rainfall * Exponential_Normalisation_Factors_Rainfall[(t + offset)] * rFsum_rainfall;
      }
      else if (rainfall_effect == "hill") {
        double mean_rainfall = Exponential_Normalisation_Factors_Rainfall[(t + offset)] * rFsum_rainfall;
        double hill_output = Hill_Function(mean_rainfall, K_Max_Hill_Rainfall, Hill_Rainfall_1, Hill_Rainfall_2);
        K_Rain = scaling_factor_rainfall * hill_output;
      }

      // Calculating K_Static
      if (Washout_Occurrence == 0) {
        K_Static = K_Max_Static;
      }

      else if (Washout_Occurrence == 1) {

        double temp_tau_static = tau_with_dt_static;
        double rFsum_statically = 0.0;
        for (int r = 0; r <= t + offset; r++) { // need to decide whether to start at beginning of rain or t = 0 ??
          calc = Exponential_Weighting_Factors_Static[r];
          rFsum_statically = rFsum_statically + calc * rF[(t + offset) - r];
        }

        // Calculating Whether Washout Occurs
        double rainfall_average_static_calc = Exponential_Normalisation_Factors_Static[(t + offset)] * rFsum_statically;
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
            else { // If Carrying Capacity Decline Follows a Hill Function
              if (within_data_marker == 0) {
                K_Static = Hill_Function((t + offset - marker_stop), K_Max_Static, washout_hill_one, washout_hill_two);
              }
              else {
                K_Static = Hill_Function((t - marker_stop), K_Max_Static, washout_hill_one, washout_hill_two);
              }
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

    ///////////////////////////////////////////////////////////////////
    //                                                               //
    //  5.  Model Running - Stepping the Model Forward One Timepoint //
    //                                                               //
    ///////////////////////////////////////////////////////////////////


    // Actual Model Running Now That the Values of K, muE and muL Have Been Calculated.
    //      Specifying the total number of events occurring for each compartment at each timepoint

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

  // Rcpp::Rcout << "Step 4: Model Running COMPLETE" << std::endl;

  // Returns Model Outputs and Other Outputs Required to Check Everything's Functional
  if(full_output == 0) {
    return Rcpp::List::create(Rcpp::Named("E_Output") = E_output,
                              Rcpp::Named("L_Output") = L_output,
                              Rcpp::Named("P_Output") = P_output,
                              Rcpp::Named("M_Output") = M_output);
  }
  else if(full_output == 1) {
    return Rcpp::List::create(Rcpp::Named("E_Output") = E_output,
                              Rcpp::Named("L_Output") = L_output,
                              Rcpp::Named("P_Output") = P_output,
                              Rcpp::Named("M_Output") = M_output,
                              Rcpp::Named("K_Rain") = k_rain_output,
                              Rcpp::Named("K_Static") = k_static_output,
                              Rcpp::Named("K_Total") = k_total_output,
                              Rcpp::Named("rainfallaverage_Kstatic") = average_rainfall_K_Static,
                              Rcpp::Named("prior K") = prior_K_Static_values,
                              Rcpp::Named("new_exp_static") = Exponential_Weighting_Factors_Static,
                              Rcpp::Named("new_exp_rainfall") = Exponential_Weighting_Factors_Rainfall,
                              Rcpp::Named("exponential_normaliser") = Exponential_Normalisation_Factors_Rainfall,
                              Rcpp::Named("within_loop") = temp_vector_Exponential_Weighting_Factors_Rainfall,
                              Rcpp::Named("outside_loop")= Exponential_Weighting_Factors_Rainfall);
  }
}

