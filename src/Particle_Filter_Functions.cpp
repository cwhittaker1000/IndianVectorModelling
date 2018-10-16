// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]

#include "Initial_State_Sampler.hpp"
#include "Mosquito_Population_Model.hpp"
#include "Particle_Filter_Functions.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include "Hill_Function.hpp"

//double Hill_Function(double rainfall, double K_max, double a, double b);


//' @export
// [[Rcpp::export]]
double min_output_particle_filter(int N, std::vector <double> rainfall, std::vector <int> obsData,
                                  Rcpp::NumericVector fitted_parameters, Rcpp::NumericVector static_parameters,
                                  Rcpp::String mortality_density_function, Rcpp::String rainfall_relationship,
                                  Rcpp::String rainfall_effect, Rcpp::String decline_type,
                                  double sampling_point, Rcpp::StringVector offset_month_vector, Rcpp::StringVector sampling_month_vector, Rcpp::String likelihood_choice) {

  // Initialising Number of Datapoints and the Initial States Vector
  int number_of_datapoints = obsData.size();
  std::vector <int> initial_states = test_initial_state_sample(fitted_parameters, static_parameters, mortality_density_function);

  // Filling the particles vector with the initial states
  boost::numeric::ublas::matrix <int> particles(N, 4);
  for (int x = 0; x < N; x++) {
    particles(x, 0) = initial_states[0];
    particles(x, 1) = initial_states[1];
    particles(x, 2) = initial_states[2];
    particles(x, 3) = initial_states[3];
  }

  // Calculating the Size/Value of the Offset Term and Inputting Its Value Into the Fitted Parameters Vector
  int offset = 0;
  for(int o = 0; o <  offset_month_vector.size(); o++) {
    if (sampling_month_vector[o] == "January" | sampling_month_vector[o] == "March" | sampling_month_vector[o] == "May" | sampling_month_vector[o] == "July" | sampling_month_vector[o] == "August" | sampling_month_vector[o] == "October" | sampling_month_vector[o] == "December") {
      offset = offset + (31 / static_parameters[0]);
    }

    else if(sampling_month_vector[o] == "April" | sampling_month_vector[o] == "June" | sampling_month_vector[o] == "September" | sampling_month_vector[o] == "November") {
      offset = offset + (30 / static_parameters[0]);
    }

    else if(sampling_month_vector[o] == "February") {
      offset = offset + (28 / static_parameters[0]);
    }

  }
  Rcpp::Rcout << "The offset is " << offset << std::endl;
  fitted_parameters[28] = offset;

  // Specifying the TIMEPOINTS to run the model between
  //    - RETURNS A VECTOR LENGTH #TIMEPOINTS + 1, DON'T FORGET THE +1 - E.G. 20 DATAPOINTS = 21 TIMEPOINTS TO RUN MODEL BETWEEN
  std::vector <int> timepoints(number_of_datapoints + 1);
  for(int p = 0; p < number_of_datapoints + 1; p++) {

    if (p == 0) {
      timepoints[p] = 0;
    }

    else if (p == 1) {

      if (sampling_month_vector[p - 1] == "January" | sampling_month_vector[p - 1] == "March" | sampling_month_vector[p - 1] == "May" | sampling_month_vector[p - 1] == "July" | sampling_month_vector[p - 1] == "August" | sampling_month_vector[p - 1] == "October" | sampling_month_vector[p - 1] == "December") {
        timepoints[p] =  round(sampling_point * (31 / static_parameters[0]));
      }

      else if(sampling_month_vector[p - 1] == "April" | sampling_month_vector[p - 1] == "June" | sampling_month_vector[p - 1] == "September" | sampling_month_vector[p - 1] == "November") {
        timepoints[p] =  round(sampling_point * (30 / static_parameters[0]));
      }

      else if(sampling_month_vector[p - 1] == "February") {
        timepoints[p] =  round(sampling_point * (28 / static_parameters[0]));
      }

    }

    else if (p > 1) {

      if (sampling_month_vector[p - 1] == "January" | sampling_month_vector[p - 1] == "March" | sampling_month_vector[p - 1] == "May" | sampling_month_vector[p - 1] == "July" | sampling_month_vector[p - 1] == "August" | sampling_month_vector[p - 1] == "October" | sampling_month_vector[p - 1] == "December") {
        timepoints[p] = timepoints[p - 1] + 31 / static_parameters[0];
      }

      else if(sampling_month_vector[p - 1] == "April" | sampling_month_vector[p - 1] == "June" | sampling_month_vector[p - 1] == "September" | sampling_month_vector[p - 1] == "November") {
        timepoints[p] = timepoints[p - 1] + 30 / static_parameters[0];
      }

      else if(sampling_month_vector[p - 1] == "February") {
        timepoints[p] = timepoints[p - 1] + 28 / static_parameters[0];
      }
    }
  }
  int total_length_in_steps = timepoints[number_of_datapoints]; // total number of days the model is run for

  // Pre-allocating storage for the MODEL OUTPUTS
  Rcpp::NumericMatrix loglikelihood(number_of_datapoints, N);
  double final_log_likelihood = 0;

  // Running the PARTICLE FILTER
  for (int i = 0; i < number_of_datapoints; i++) {  // need to check that it's -2 and not -1. Think I'm right

    // Pre-allocating storage for TEMPORARY OUTPUTS
    std::vector <double> particle_weights(N);
    boost::numeric::ublas::matrix <double> final_output_particles_before_weighting(N, 4);

    // Iterating over the NUMBER OF PARTICLES we're using
    for (int j = 0; j < N; j++) {

      // Initialising E0, L0, P0 and M0
      fitted_parameters[24] = particles(j, 0); // E  // need to make sure this is doing what I want it to!
      fitted_parameters[25] = particles(j, 1); // L
      fitted_parameters[26] = particles(j, 2); // P
      fitted_parameters[27] = particles(j, 3); // M

      // RUNS THE MODEL for a SINGLE PARTICLE between the TWO TIMEPOINTS
      Rcpp::List particle_output = general_mosquito_population_model(timepoints[i], timepoints[i + 1], fitted_parameters, static_parameters, rainfall, mortality_density_function, rainfall_relationship, rainfall_effect, decline_type);

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

      if (likelihood_choice == "poisson") {
        particle_weights[j] = Poisson(obsData[i], current_M, fitted_parameters[10]);
      }
      else if (likelihood_choice == "negative binomial") {
        particle_weights[j] = Negative_Binomial(obsData[i], current_M, fitted_parameters[9], fitted_parameters[10]); // double check I've got the overdispersion and population fraction the right way around!!
      }

    }

    // Adding the unnormalised particle weights to the loglikelihood vector
    for (int k = 0; k < N; k++) {
      loglikelihood(i, k) = particle_weights[k];
    }

    // Normalising particle weights
    std::vector <double> normalised_particle_weights = Particle_Weight_Normalisation(particle_weights);

    // Create weighted distribution for resampling, resampling and then generating the required outputs
    std::vector <int> rows = weighted_sampling_with_replacement(N, normalised_particle_weights, 1, N);

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


//' @export
// [[Rcpp::export]]
Rcpp::List full_output_particle_filter(int N, std::vector <double> rainfall, std::vector <int> obsData,
                                       Rcpp::NumericVector fitted_parameters, Rcpp::NumericVector static_parameters,
                                       Rcpp::String mortality_density_function, Rcpp::String rainfall_relationship,
                                       Rcpp::String rainfall_effect, Rcpp::String decline_type,
                                       double sampling_point, Rcpp::StringVector offset_month_vector, Rcpp::StringVector sampling_month_vector, Rcpp::String likelihood_choice) {

  // Initialising Number of Datapoints and the Initial States Vector
  int number_of_datapoints = obsData.size();
  std::vector <int> initial_states = test_initial_state_sample(fitted_parameters, static_parameters, mortality_density_function);

  // Filling the particles vector with the initial states
  Rcpp::NumericMatrix particles(N, 4);
  for (int x = 0; x < N; x++) {
    particles(x, 0) = initial_states[0];
    particles(x, 1) = initial_states[1];
    particles(x, 2) = initial_states[2];
    particles(x, 3) = initial_states[3];
  }

  // Calculating the Size/Value of the Offset Term
  int offset = 0;
  for(int o = 0; o <  offset_month_vector.size(); o++) {
    if (sampling_month_vector[o] == "January" | sampling_month_vector[o] == "March" | sampling_month_vector[o] == "May" | sampling_month_vector[o] == "July" | sampling_month_vector[o] == "August" | sampling_month_vector[o] == "October" | sampling_month_vector[o] == "December") {
      offset = offset + (31 / static_parameters[0]);
    }

    else if(sampling_month_vector[o] == "April" | sampling_month_vector[o] == "June" | sampling_month_vector[o] == "September" | sampling_month_vector[o] == "November") {
      offset = offset + (30 / static_parameters[0]);
    }

    else if(sampling_month_vector[o] == "February") {
      offset = offset + (28 / static_parameters[0]);
    }

  }
  Rcpp::Rcout << "The offset is " << offset << std::endl;

  // Specifying the TIMEPOINTS to run the model between
  //    - RETURNS A VECTOR LENGTH #TIMEPOINTS + 1, DON'T FORGET THE +1 - E.G. 20 DATAPOINTS = 21 TIMEPOINTS TO RUN MODEL BETWEEN
  Rcpp::NumericVector timepoints(number_of_datapoints + 1);
  for(int p = 0; p < number_of_datapoints + 1; p++) {

    if (p == 0) {
      timepoints[p] = 0;
    }

    else if (p == 1) {

      if (sampling_month_vector[p - 1] == "January" | sampling_month_vector[p - 1] == "March" | sampling_month_vector[p - 1] == "May" | sampling_month_vector[p - 1] == "July" | sampling_month_vector[p - 1] == "August" | sampling_month_vector[p - 1] == "October" | sampling_month_vector[p - 1] == "December") {
        timepoints[p] =  round(sampling_point * (31 / static_parameters[0]));
      }

      else if(sampling_month_vector[p - 1] == "April" | sampling_month_vector[p - 1] == "June" | sampling_month_vector[p - 1] == "September" | sampling_month_vector[p - 1] == "November") {
        timepoints[p] =  round(sampling_point * (30 / static_parameters[0]));
      }

      else if(sampling_month_vector[p - 1] == "February") {
        timepoints[p] =  round(sampling_point * (28 / static_parameters[0]));
      }

    }

    else if (p > 1) {

      if (sampling_month_vector[p - 1] == "January" | sampling_month_vector[p - 1] == "March" | sampling_month_vector[p - 1] == "May" | sampling_month_vector[p - 1] == "July" | sampling_month_vector[p - 1] == "August" | sampling_month_vector[p - 1] == "October" | sampling_month_vector[p - 1] == "December") {
        timepoints[p] = timepoints[p - 1] + 31 / static_parameters[0];
      }

      else if(sampling_month_vector[p - 1] == "April" | sampling_month_vector[p - 1] == "June" | sampling_month_vector[p - 1] == "September" | sampling_month_vector[p - 1] == "November") {
        timepoints[p] = timepoints[p - 1] + 30 / static_parameters[0];
      }

      else if(sampling_month_vector[p - 1] == "February") {
        timepoints[p] = timepoints[p - 1] + 28 / static_parameters[0];
      }
    }
  }
  int total_length_in_steps = timepoints[number_of_datapoints]; // total number of days the model is run for

  // Pre-allocating storage for the MODEL OUTPUTS
  Rcpp::NumericMatrix E_after(N, total_length_in_steps + 100); // plus 100 to deal with random edge cases
  Rcpp::NumericMatrix L_after(N, total_length_in_steps + 100); // becomes NAs in model outputs
  Rcpp::NumericMatrix P_after(N, total_length_in_steps + 100);
  Rcpp::NumericMatrix M_after(N, total_length_in_steps + 100);
  Rcpp::NumericMatrix loglikelihood(number_of_datapoints, N);
  Rcpp::NumericMatrix M_values(number_of_datapoints, N);
  Rcpp::NumericMatrix particle_weights_proper(number_of_datapoints, N);

  double final_log_likelihood = 0;

  // Running the PARTICLE FILTER
  for (int i = 0; i < number_of_datapoints; i++) {  // need to check that it's -2 and not -1. Think I'm right

    // Pre-allocating storage for TEMPORARY OUTPUTS
    std::vector <double> particle_weights(N);
    Rcpp::NumericMatrix final_output_particles_before_weighting(N, 4);
    Rcpp::NumericMatrix temp_E_Output(N, timepoints[i + 1] - timepoints[i]);
    Rcpp::NumericMatrix temp_L_Output(N, timepoints[i + 1] - timepoints[i]);
    Rcpp::NumericMatrix temp_P_Output(N, timepoints[i + 1] - timepoints[i]);
    Rcpp::NumericMatrix temp_M_Output(N, timepoints[i + 1] - timepoints[i]);

    // Iterating over the NUMBER OF PARTICLES we're using
    for (int j = 0; j < N; j++) {

      // Initialising E0, L0, P0 and M0
      fitted_parameters[24] = particles(j, 0); // E  // need to make sure this is doing what I want it to!
      fitted_parameters[25] = particles(j, 1); // L
      fitted_parameters[26] = particles(j, 2); // P
      fitted_parameters[27] = particles(j, 3); // M

      // RUNS THE MODEL for a SINGLE PARTICLE between the TWO TIMEPOINTS
      Rcpp::List particle_output = general_mosquito_population_model(timepoints[i], timepoints[i + 1], fitted_parameters, static_parameters, rainfall, mortality_density_function, rainfall_relationship, rainfall_effect, decline_type);

      // Extracts the FULL OUTPUT for each compartment for that particle
      Rcpp::NumericVector E_output = particle_output["E_Output"];
      Rcpp::NumericVector L_output = particle_output["L_Output"];
      Rcpp::NumericVector P_output = particle_output["P_Output"];
      Rcpp::NumericVector M_output = particle_output["M_Output"];

      // Adds the FULL OUTPUT to the TEMPORARY OUTPUT matrices
      for (int l = 0; l < (timepoints[i + 1] - timepoints[i]); l++) {
        temp_E_Output(j, l) = E_output[l]; // need to know howw to add NumericVectors/std::vectors to NumericMatrix's
        temp_L_Output(j, l) = L_output[l];
        temp_P_Output(j, l) = P_output[l];
        temp_M_Output(j, l) = M_output[l];
      }

      // Extracts the FINAL OUTPUT VALUE for each compartment for that particle
      final_output_particles_before_weighting(j, 0) = E_output[(timepoints[i + 1] - timepoints[i]) - 1];
      final_output_particles_before_weighting(j, 1) = L_output[(timepoints[i + 1] - timepoints[i]) - 1];
      final_output_particles_before_weighting(j, 2) = P_output[(timepoints[i + 1] - timepoints[i]) - 1];
      final_output_particles_before_weighting(j, 3) = M_output[(timepoints[i + 1] - timepoints[i]) - 1];

      // Calculate the PARTICLE WEIGHTS // MAKE SURE TO ENSURE POP_FRAC AND OVERDISP ARE AT THE CORRECT POSITIONS IN THE VECTOR TO BE INCLUDED HERE
      double current_M = final_output_particles_before_weighting(j, 3);
      M_values(i , j) = current_M;

      if (likelihood_choice == "poisson") {
        particle_weights[j] = Poisson(obsData[i], current_M, fitted_parameters[10]);
      }
      else if (likelihood_choice == "negative binomial") {
        particle_weights[j] = Negative_Binomial(obsData[i], current_M, fitted_parameters[9], fitted_parameters[10]); // double check I've got the overdispersion and population fraction the right way around!!
      }

      particle_weights_proper(i, j) = particle_weights[j];
    }

    // Adding the unnormalised particle weights to the loglikelihood vector
    for (int k = 0; k < N; k++) {
      loglikelihood(i, k) = particle_weights[k];
    }

    // Normalising particle weights
    std::vector <double> normalised_particle_weights = Particle_Weight_Normalisation(particle_weights);

    // Create weighted distribution for resampling, resampling and then generating the required outputs
    std::vector <int> rows = weighted_sampling_with_replacement(N, normalised_particle_weights, 1, N);

    for (int d = 0; d < N; d++) {
      particles(d, 0) = final_output_particles_before_weighting(rows[d], 0);
      particles(d, 1) = final_output_particles_before_weighting(rows[d], 1);
      particles(d, 2) = final_output_particles_before_weighting(rows[d], 2);
      particles(d, 3) = final_output_particles_before_weighting(rows[d], 3);

      int z = 0;

      for (int n = timepoints[i]; n < timepoints[i + 1]; n++) {
        E_after(d, n) = temp_E_Output(rows[d], z);
        L_after(d, n) = temp_L_Output(rows[d], z);
        P_after(d, n) = temp_P_Output(rows[d], z);
        M_after(d, n) = temp_M_Output(rows[d], z);
        z = z + 1;
      }
    }
  }

  for (int l = 0; l < number_of_datapoints; l++) {
    final_log_likelihood = final_log_likelihood + Rcpp::mean(loglikelihood(l, Rcpp::_));
  }

  return Rcpp::List::create(Rcpp::Named("E_Output") = E_after,
                            Rcpp::Named("L_Output") = L_after,
                            Rcpp::Named("P_Output") = P_after,
                            Rcpp::Named("M_Output") = M_after,
                            Rcpp::Named("Timepoints") = timepoints,
                            Rcpp::Named("Loglikelihood") = loglikelihood,
                            Rcpp::Named("M_Values") = M_values,
                            Rcpp::Named("particle weights") = particle_weights_proper,
                            Rcpp::Named("initial states") = initial_states,
                            Rcpp::Named("final log lik") = final_log_likelihood);
}


// k Observed data point
// n Simulated data point
// r overdispersion parameter
// p fraction of population

//' @export
// [[Rcpp::export]]
double Negative_Binomial(double k, double n, double r, double p) {
  double m = p * n;
  double res = (lgamma(r + k) - (lgamma(k + 1) + lgamma(r))) + k * (log(m) - (log(r + m))) + r * (log(r) - log(r + m));
  return(res);
}

//' @export
// [[Rcpp::export]]
double Poisson(double observed_data, double model_output, double population_fraction) {

  // bool log = true;
  // double res = Rcpp::dpois(observed_data, model_output * population_fraction, log);
  double res = R::dpois(observed_data, model_output * population_fraction, 1);
  return(res);

}

//' @export
// [[Rcpp::export]]
std::vector <double> Particle_Weight_Normalisation(std::vector <double> particle_weights) {

  // Declaring a vector for the processed particle weights
  std::vector <double> particle_weights_normal_scale(particle_weights.size(), 0);
  std::vector <double> particle_weights_normal_scale_normalised(particle_weights.size(), 0);

  // Checking the passed values, and adapting them if they're out of relevant range
  for (int i = 0; i < particle_weights.size(); i++) {
    double weight = particle_weights[i];
    if (weight < -16) {
      weight = -16;
    }
    particle_weights_normal_scale[i] = exp(weight);
  }

  // Accumulate takes the initial point, the end point and the initial value you want to sum to, from and together.
  double summed_particle_weights = accumulate(particle_weights_normal_scale.begin(), particle_weights_normal_scale.end(), 0.0);

  // Loop to return a normalised set of particle weights, not on the log scale
  for (int i = 0; i < particle_weights_normal_scale.size(); i++) {
    particle_weights_normal_scale_normalised[i] = particle_weights_normal_scale[i] / summed_particle_weights;
  }

  return(particle_weights_normal_scale_normalised);
}


int rbinomial1(int trials, double p) {
  return(R::rbinom(trials,p));
}

//' @export
// [[Rcpp::export]]
std::vector<int> weighted_sampling_with_replacement(int n, std::vector<double> prob, double p_tot, int K)
{
  int k;
  int draw;
  std::vector<int> output;

  for (k = 0; k < K - 1; k++) {

    draw = rbinomial1(n, (prob[k] / p_tot));
    n -= draw;
    while (draw > 0) {
      output.emplace_back(k);
      draw--;
    }

    if (n <= 0) {
      return(output);
    }
    else {
      p_tot -= prob[k];
    }
  }
  while (n > 0) {
    output.emplace_back(K - 1);
    n--;
  }
  return(output);
}


