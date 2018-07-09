#include "Mosquito_Population_Model.hpp"
#include "Negative_Binomial.hpp"
#include "Particle_Weight_Normalisation.hpp"
#include "weighted_sampling_with_replacement.hpp"
#include "Initial_State_Sample.hpp"


// [[Rcpp::export]]
double min_output_particle_filter(int N, std::vector <double> rainfall, std::vector <int> obsData, int number_of_datapoints, int data_timeframe, Rcpp::NumericVector fitted_parameters, Rcpp::NumericVector static_parameters, Rcpp::String density_function) {

  // Specifying the INITIAL VALUES E, L, P and M and generating the PARTICLES
  Rcpp::NumericVector initial_states = initial_state_sample(fitted_parameters, static_parameters);
  Rcpp::NumericMatrix particles(N, 4);
  for (int x = 0; x < N; x++) {
    particles(x, 0) = initial_states[0];
    particles(x, 1) = initial_states[1];
    particles(x, 2) = initial_states[2];
    particles(x, 3) = initial_states[3];
  }

  // Specifying the TIMEPOINTS to run the model between
  //    - RETURNS A VECTOR LENGTH #TIMEPOINTS + 1, DON'T FORGET THE +1 - E.G. 20 DATAPOINTS = 21 TIMEPOINTS TO RUN MODEL BETWEEN
  int total_length_in_steps = (data_timeframe * number_of_datapoints) / static_parameters[0]; // total number of days the model is run for
  int length_of_one_step = data_timeframe / static_parameters[0]; // static_parameters[0] is dt.
  Rcpp::NumericVector timepoints(number_of_datapoints + 1)                   ;
  for(int p = 0; p < number_of_datapoints + 1; p++) {
    if (p == 0) {
      timepoints[p] = 0;
    }
    else {
      timepoints[p] = (length_of_one_step/2) + length_of_one_step * (p - 1);
    }
  }

  // Pre-allocating storage for the MODEL OUTPUTS
  Rcpp::NumericMatrix loglikelihood(number_of_datapoints, N);
  double final_log_likelihood = 0;

  // Running the PARTICLE FILTER
  for (int i = 0; i < number_of_datapoints; i++) {  // need to check that it's -2 and not -1. Think I'm right

    // Pre-allocating storage for TEMPORARY OUTPUTS
    std::vector <double> particle_weights(N);
    Rcpp::NumericMatrix final_output_particles_before_weighting(N, 4);

    // Iterating over the NUMBER OF PARTICLES we're using
    for (int j = 0; j < N; j++) {

      // Initialising E0, L0, P0 and M0
      fitted_parameters[12] = particles(j, 0); // E  // need to make sure this is doing what I want it to!
      fitted_parameters[13] = particles(j, 1); // L
      fitted_parameters[14] = particles(j, 2); // P
      fitted_parameters[15] = particles(j, 3); // M

      // RUNS THE MODEL for a SINGLE PARTICLE between the TWO TIMEPOINTS
      Rcpp::List particle_output = mosquito_population_model(timepoints[i], timepoints[i + 1], fitted_parameters, static_parameters, rainfall, density_function);

      // Extracts the FULL OUTPUT for each compartment for that particle
      Rcpp::NumericVector E_output = particle_output["E_Output"];
      Rcpp::NumericVector L_output = particle_output["L_Output"];
      Rcpp::NumericVector P_output = particle_output["P_Output"];
      Rcpp::NumericVector M_output = particle_output["M_Output"];

      // Extracts the FINAL OUTPUT VALUE for each compartment for that particle
      final_output_particles_before_weighting(j, 0) = E_output[(timepoints[i + 1] - timepoints[i]) - 1];
      final_output_particles_before_weighting(j, 1) = L_output[(timepoints[i + 1] - timepoints[i]) - 1];
      final_output_particles_before_weighting(j, 2) = P_output[(timepoints[i + 1] - timepoints[i]) - 1];
      final_output_particles_before_weighting(j, 3) = M_output[(timepoints[i + 1] - timepoints[i]) - 1];

      // Calculate the PARTICLE WEIGHTS // MAKE SURE TO ENSURE POP_FRAC AND OVERDISP ARE AT THE CORRECT POSITIONS IN THE VECTOR TO BE INCLUDED HERE
      double current_M = final_output_particles_before_weighting(j, 3);
      particle_weights[j] = Negative_Binomial(obsData[i], current_M, fitted_parameters[10], fitted_parameters[11]); // double check I've got the overdispersion and population fraction the right way around!!
    }

    // Adding the unnormalised particle weights to the loglikelihood vector
    for (int k = 0; k < N; k++) {
      loglikelihood(i, k) = particle_weights[k];
    }

    // Normalising particle weights
    std::vector <double> normalised_particle_weights = Particle_Weight_Normalisation(particle_weights);

    // Create weighted distribution for resampling, resampling and then generating the required outputs
    std::vector <int> rows = weighted_sampling(N, normalised_particle_weights, 1, N);

    for (int d = 0; d < N; d++) {
      particles(d, 0) = final_output_particles_before_weighting(rows[d], 0);
      particles(d, 1) = final_output_particles_before_weighting(rows[d], 1);
      particles(d, 2) = final_output_particles_before_weighting(rows[d], 2);
      particles(d, 3) = final_output_particles_before_weighting(rows[d], 3);
    }
  }

  for (int l = 0; l < number_of_datapoints; l++) {
    final_log_likelihood = final_log_likelihood + Rcpp::mean(loglikelihood(l,Rcpp::_));
  }

  return(final_log_likelihood);
}
