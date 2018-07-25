#include "MCMC_Functions.hpp"

//' @export
// [[Rcpp::export]]
Rcpp::List runMCMC_seq_props(int start_sd_adaptation, // Time to start SD adaptation - relates to SD adaptation
                             int end_sd_adaptation,
                             std::vector <double> acceptance_ratio, // Desired acceptance ratio - relates to SD adaptation
                             std::vector<double> sd_proposals, // SD for proposal function
                             int number_of_iterations, // Number of MCMC iterations
                             Rcpp::NumericVector model_parameters, // Initial parameter values for values that are going to be fitted
                             Rcpp::NumericVector static_parameters, // Values for static parameter values
                             int N, // Number of particles - relates to Particle Filter/Model
                             std::vector <double> rainfall, // Rainfall recorded - relates to Particle Filter/Model
                             std::vector <int> obsData, // Observed data - relates to Particle Filter/Model
                             int number_of_datapoints, // Number of datapoints -  relates to Particle Filter/Model
                             int data_timeframe, // Time between datapoints - relates to Particle Filter/Model
                             Rcpp::String density_function, // Which density function regulating mortality - relates to Particle Filter/Model
                             Rcpp::String prior_choice,
                             Rcpp::LogicalVector fitted_yn) {

  // Storage for the MCMC chains
  Rcpp::NumericMatrix MCMC_chain_output(number_of_iterations + 1, Rcpp::sum(fitted_yn));

  double accepted_variable;

  // Storage for various things required to check the stuff
  Rcpp::NumericMatrix accepted_checker(number_of_iterations + 1, Rcpp::sum(fitted_yn));
  Rcpp::NumericMatrix acceptance_ratio_tracker(number_of_iterations + 1, Rcpp::sum(fitted_yn));
  Rcpp::NumericMatrix proposal_sd_tracker(number_of_iterations + 1, Rcpp::sum(fitted_yn));
  Rcpp::NumericVector current_acceptance_ratio(Rcpp::sum(fitted_yn));

  Rcpp::StringVector parameter_names = fitted_yn.names();
  Rcpp::StringVector naming(Rcpp::sum(fitted_yn));
  Rcpp::NumericVector indexer(Rcpp::sum(fitted_yn));

  // Assigning initial values to the first row of the MCMC output
  int p = 0;
  for (int d = 0; d < fitted_yn.size(); d++) {

    if (fitted_yn[d] == TRUE) {
      Rcpp::String to_be_added = parameter_names[d];
      naming[p] = to_be_added;
      indexer[p] = d;
      MCMC_chain_output(0, p) = model_parameters[to_be_added];
      p = p + 1;

    }
  }
  colnames(MCMC_chain_output) = naming;

  // Calculating the posterior likelihood for the initial parameter values
  // posterior usually takes a new value for a single parameter, replaces current value of that parameter in
  // current parameters, then calculates the posterior. By using parameter_values[0] as my "new" parameter,
  // allows me to calculate the posterior for my initial parameter values.
  double current_posterior_likelihood = posterior_seq_proposals(N, rainfall, obsData, number_of_datapoints,
                                                                data_timeframe, model_parameters, static_parameters, density_function,
                                                                model_parameters[0], 0, prior_choice, fitted_yn); // whole bunch of inputs here, see if I can streamline this

  // Initialising variables required for the standard deviation tuning
  std::vector <double> acceptances(fitted_yn.size(), 0); // double check this generates a vector with the right dimensions and contents
  std::vector <double> rejections(fitted_yn.size(), 0); // double check this generates a vector with the right dimensions and contents

  // Creating a vector to store and sequentially update model parameters as the MCMC iterates
  Rcpp::NumericVector model_parameters_for_MCMC = Rcpp::NumericVector::create(Rcpp::Named("dE") = model_parameters[0], Rcpp::Named("dL") = model_parameters[1],
                                                                              Rcpp::Named("dP") = model_parameters[2], Rcpp::Named("muE0") = model_parameters[3],
                                                                              Rcpp::Named("muL0") = model_parameters[4], Rcpp::Named("muP") = model_parameters[5],
                                                                              Rcpp::Named("muM") = model_parameters[6], Rcpp::Named("lambda") = model_parameters[7],
                                                                              Rcpp::Named("tau") = model_parameters[8], Rcpp::Named("beta") = model_parameters[9],
                                                                              Rcpp::Named("overdisp") = model_parameters[10], Rcpp::Named("pop_frac") = model_parameters[11],
                                                                              Rcpp::Named("E") = 0.0, Rcpp::Named("L") = 0.0,
                                                                              Rcpp::Named("P") = 0.0, Rcpp::Named("M") = 0.0);

  // Iterating over each parameter for a specified number of iterations
  for (int i = 0; i < number_of_iterations; i++) {

    for (int j = 0; j < Rcpp::sum(fitted_yn); j++) { // NOTE: MEANS THAT THE ORDER MATTERS CURRENTLY RE FLEXIBILITY IN TERMS OF WHICH PARAMETERS ARE FITTED

      // Creating a vector of fitted parameter values to pass into the posterior function
      // Changes the model parameters in the vector to their current values in the chain, whilst leaving static model parameters unchanged
      // Note again that I can turn off parameters, but only in sequence right to left (can't turn off specific ones in middle of vector currently)
      for (int k = 0; k < indexer.size(); k++) {
        int index_for_adding = indexer[k];
        model_parameters_for_MCMC[index_for_adding] = MCMC_chain_output(i, k);
        }

      int index = indexer[j];
      model_parameters_for_MCMC.names() = Rcpp::StringVector::create("dE", "dL", "dP", "muE0", "muL0", "muP", "muM", "lambda", "tau", "beta", "overdisp", "pop_frac", "E", "L", "P", "M");

      // Generating the proposed parameter value
      double current_parameter_value = MCMC_chain_output(i, j);
      double proposed_parameter_value = seq_proposal_function(sd_proposals[index], current_parameter_value);

      // Assessing the likelihood of the proposed parameter value compared to the current parameter value
      double proposed_posterior_likelihood = posterior_seq_proposals(N, rainfall, obsData, number_of_datapoints,
                                                                     data_timeframe, model_parameters_for_MCMC, static_parameters, density_function,
                                                                     proposed_parameter_value, index, prior_choice, fitted_yn);

      double likelihood_ratio = exp(proposed_posterior_likelihood - current_posterior_likelihood);

      // Deciding whether to accept or reject, updating chain and current_posterior_likelihood as appropriate
      if(R::runif(0, 1) < likelihood_ratio) { // double check this is correct
        MCMC_chain_output(i + 1, j) = proposed_parameter_value;
        current_posterior_likelihood = proposed_posterior_likelihood;
        acceptances[j] = acceptances[j] + 1;
        accepted_checker(i, j) = 2;
        accepted_variable = 1;
      }
      else {
        MCMC_chain_output(i + 1, j) = current_parameter_value;
        current_posterior_likelihood = current_posterior_likelihood; // remove this as it's superfluous- just to tell myself what's going on
        rejections[j] = rejections[j] + 1;
        accepted_checker(i, j) = 1;
        accepted_variable = 0;
      }

      // Tracking the acceptance ratio
      current_acceptance_ratio[j] = (acceptances[j] / (i+1));
      acceptance_ratio_tracker(i, j) = current_acceptance_ratio[j];

      // Adapting the standard deviation of parameter proposals to ensure a decent acceptance ratio
      if (i >= start_sd_adaptation & i < end_sd_adaptation) {

        sd_proposals[index] = seq_proposal_SD_adapter(accepted_variable, i, start_sd_adaptation, sd_proposals[index]);
        proposal_sd_tracker(i, j) = sd_proposals[index];

      }

    }

  }

  return(Rcpp::List::create(Rcpp::Named("MCMC_Output") = MCMC_chain_output,
                            Rcpp::Named("Proposal_Sd_final") = sd_proposals,
                            Rcpp::Named("Acceptance_Ratio_Final") = current_acceptance_ratio,
                            Rcpp::Named("Acceptance_Checker") = accepted_checker,
                            Rcpp::Named("Acceptances") = acceptances,
                            Rcpp::Named("Rejections") = rejections,
                            Rcpp::Named("Acceptance Ratio Tracker") = acceptance_ratio_tracker,
                            Rcpp::Named("Proposal SD Tracker") = proposal_sd_tracker,
                            Rcpp::Named("naming") = naming,
                            Rcpp::Named("indexer") = indexer));
}

