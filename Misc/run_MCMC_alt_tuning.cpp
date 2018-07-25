#include "MCMC_Functions.hpp"
#include "Proposal_Functions.hpp"

//' @export
// [[Rcpp::export]]
Rcpp::List runMCMC_alt_tune(int start_sd_adaptation, // Time to start SD adaptation - relates to SD adaptation
                            int end_sd_adaptation,
                            std::vector <double> max_sd, // Maximum SD allowed for parameter - relates to SD adaptation
                            std::vector <double> acceptance_ratio, // Desired acceptance ratio - relates to SD adaptation
                            std::vector<double> sd_proposals, // SD for proposal function
                            int number_of_iterations, // Number of MCMC iterations
                            Rcpp::StringVector parameters_to_be_fitted, // String of parameters (their names) to be fitted
                            Rcpp::NumericVector fitted_parameters, // Initial parameter values for values that are going to be fitted
                            Rcpp::NumericVector static_parameters, // Values for static parameter values
                            int N, // Number of particles - relates to Particle Filter/Model
                            std::vector <double> rainfall, // Rainfall recorded - relates to Particle Filter/Model
                            std::vector <int> obsData, // Observed data - relates to Particle Filter/Model
                            int number_of_datapoints, // Number of datapoints -  relates to Particle Filter/Model
                            int data_timeframe, // Time between datapoints - relates to Particle Filter/Model
                            Rcpp::String density_function,
                            Rcpp::String prior_choice,
                            Rcpp::LogicalVector fitted_yn) { // Which density function regulating mortality - relates to Particle Filter/Model

  // Storage for the MCMC chains
  Rcpp::NumericMatrix MCMC_chain_output(number_of_iterations + 1, parameters_to_be_fitted.size());
  colnames(MCMC_chain_output) = parameters_to_be_fitted;

  // Storage for various things required to check the stuff
  Rcpp::NumericMatrix accepted_checker(number_of_iterations + 1, parameters_to_be_fitted.size());
  Rcpp::NumericMatrix acceptance_ratio_tracker(number_of_iterations + 1, parameters_to_be_fitted.size());
  Rcpp::NumericMatrix proposal_sd_tracker(number_of_iterations + 1, parameters_to_be_fitted.size());
  Rcpp::NumericVector current_acceptance_ratio(parameters_to_be_fitted.size());

  // Assigning initial values to the first row of the MCMC output
  for (int d = 0; d < parameters_to_be_fitted.size(); d++) {
    Rcpp::String name_parameter_to_be_added = parameters_to_be_fitted[d];
    MCMC_chain_output(0, d) = fitted_parameters[name_parameter_to_be_added];
  }

  // Calculating the posterior likelihood for the initial parameter values
  // posterior usually takes a new value for a single parameter, replaces current value of that parameter in
  // current parameters, then calculates the posterior. By using parameter_values[0] as my "new" parameter,
  // allows me to calculate the posterior for my initial parameter values.
  double current_posterior_likelihood = posterior(N, rainfall, obsData, number_of_datapoints,
                                                  data_timeframe, fitted_parameters, static_parameters, density_function,
                                                  fitted_parameters[0], 0, prior_choice, fitted_yn); // whole bunch of inputs here, see if I can streamline this

  // Initialising variables required for the standard deviation tunning
  std::vector <double> acceptances(parameters_to_be_fitted.size(), 0); // double check this generates a vector with the right dimensions and contents
  std::vector <double> rejections(parameters_to_be_fitted.size(), 0); // double check this generates a vector with the right dimensions and contents

  // Iterating over each parameter for a specified number of iterations
  for (int i = 0; i < number_of_iterations; i++) {

    for (int j = 0; j < parameters_to_be_fitted.size(); j++) {

      // Creating a vector of fitted parameter values to pass into the posterior function
      Rcpp::NumericVector parameters_for_fitting(parameters_to_be_fitted.size() + 4);
      for (int k = 0; k < parameters_for_fitting.size(); k++) {
        if (k < parameters_to_be_fitted.size()) {
          parameters_for_fitting[k] = MCMC_chain_output(i, k);
        }
        else if (k == parameters_to_be_fitted.size()) {
          parameters_for_fitting[k] = 0; //E
        }
        else if (k == parameters_to_be_fitted.size() + 1) {
          parameters_for_fitting[k] = 0; //L
        }
        else if (k == parameters_to_be_fitted.size() + 2) {
          parameters_for_fitting[k] = 0; //P
        }
        else if (k == parameters_to_be_fitted.size() + 3) {
          parameters_for_fitting[k] = 0; //M
        }
      }

      // Generating the proposed parameter value
      double current_parameter_value = MCMC_chain_output(i, j);
      double proposed_parameter_value = proposal_function(sd_proposals[j], current_parameter_value);

      // Assessing the likelihood of the proposed parameter value compared to the current parameter value
      double proposed_posterior_likelihood = posterior(N, rainfall, obsData, number_of_datapoints,
                                                       data_timeframe, parameters_for_fitting, static_parameters, density_function,
                                                       proposed_parameter_value, 0, prior_choice, fitted_yn);

      double likelihood_ratio = exp(proposed_posterior_likelihood - current_posterior_likelihood);

      // Deciding whether to accept or reject, updating chain and current_posterior_likelihood as appropriate
      if(R::runif(0, 1) < likelihood_ratio) { // double check this is correct
        MCMC_chain_output(i + 1, j) = proposed_parameter_value;
        current_posterior_likelihood = proposed_posterior_likelihood;
        acceptances[j] = acceptances[j] + 1;
        accepted_checker(i, j) = 2;
      }
      else {
        MCMC_chain_output(i + 1, j) = current_parameter_value;
        current_posterior_likelihood = current_posterior_likelihood; // remove this as it's superfluous- just to tell myself what's going on
        rejections[j] = rejections[j] + 1;
        accepted_checker(i, j) = 1;
      }

      // Tracking the acceptance ratio
      current_acceptance_ratio[j] = (acceptances[j] / (i+1));
      acceptance_ratio_tracker(i, j) = current_acceptance_ratio[j];

      // Adapting the standard deviation of parameter proposals to ensure a decent acceptance ratio
      if (i >= start_sd_adaptation & i < end_sd_adaptation) {

        sd_proposals[j] = proposal_SD_adapter_mark_two(sd_proposals[j], current_acceptance_ratio[j], max_sd[j], 0.98, i, start_sd_adaptation);
        proposal_sd_tracker(i, j) = sd_proposals[j];

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
                            Rcpp::Named("Proposal SD Tracker") = proposal_sd_tracker));
}
