// [[Rcpp::depends(RcppArmadillo)]]
#include "MCMC_Functions.hpp"
#include "MVN_Sampler.hpp"

//' @export
// [[Rcpp::export]]
Rcpp::List testerr(int start_sd_adaptation, // Time to start covariance matrix adaptation
                               int end_sd_adaptation, // Time to stop covariance matrix adaptation
                               int number_of_iterations, // Number of MCMC iterations
                               std::vector <double> initial_sds, // Initial SDs to fill the covariance matrix diag with
                               Rcpp::NumericVector model_parameters, // Variable model parameter values - can be fitted
                               Rcpp::NumericVector static_parameters, // Static model parameter values - not fitted
                               int N, // Number of particles
                               std::vector <double> rainfall, // Rainfall recorded
                               std::vector <int> obsData, // Observed data
                               int number_of_datapoints, // Number of datapoints
                               int data_timeframe, // Time between datapoints
                               Rcpp::String density_function, // Density function regulating mortality
                               Rcpp::String prior_choice, // Prior choice
                               Rcpp::LogicalVector fitted_yn) {  // Which parameters to be fitted

  // Storage for the MCMC chains, tracking acceptances/rejections etc
  Rcpp::NumericMatrix MCMC_chain_output(number_of_iterations + 1, Rcpp::sum(fitted_yn));
  double acceptances = 0;
  double rejections = 0;
  double current_acceptance_ratio = 0;
  Rcpp::NumericVector acceptance_ratio_tracker(number_of_iterations + 1);
  Rcpp::NumericVector accepted_checker(number_of_iterations + 1);

  // Storage for things that generate flexibility in choosing which parameters to fit
  Rcpp::StringVector parameter_names = fitted_yn.names();
  Rcpp::StringVector naming(Rcpp::sum(fitted_yn));
  Rcpp::NumericVector indexer(Rcpp::sum(fitted_yn));

  // Assigning initial values to the first row of the MCMC output and naming that output's columns
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

  // Defining variables involved in the Adaptive MCMC part CHECK WHETHER SIZE NEEDS TO BE INITIALISED
  // double accepted_variable;
  arma::mat current_covariance_matrix(Rcpp::sum(fitted_yn), Rcpp::sum(fitted_yn));
  for (int x = 0; x < Rcpp::sum(fitted_yn); x++) {
    int index = indexer[x];
    for (int y = 0; y < Rcpp::sum(fitted_yn); y++) {
      if (y == x) {
        current_covariance_matrix(x, y) = initial_sds[index];
      }
      else {
        current_covariance_matrix(x, y) = 0;
      }
    }
  }
  arma::mat initial_covariance_matrix_checker = current_covariance_matrix;

  //NOT SURE THIS IS WHAT I WANT TO BE ASSIGNING AS THE INITIAL VALUES FOR MU TBH.
  // this needs to be a row vector at any rate
  arma::mat current_mu(1, Rcpp::sum(fitted_yn));
  for (int x = 0; x < Rcpp::sum(fitted_yn); x++) {
    int specific_index = indexer[x];
    current_mu(0, x) = model_parameters[specific_index];
  }
  arma::mat initial_mu_checker = current_mu;
  double current_scaling_factor = 1;


  // Calculating the posterior likelihood for the initial parameter values
  double current_posterior_likelihood = posterior_joint_proposals(N, rainfall, obsData, number_of_datapoints,
                                                                  data_timeframe, model_parameters, static_parameters, density_function,
                                                                  prior_choice, fitted_yn); // whole bunch of inputs here, see if I can streamline this


  // Creating a vector to store and sequentially update model parameters as the MCMC iterates
  Rcpp::NumericVector model_parameters_for_MCMC = Rcpp::NumericVector::create(Rcpp::Named("dE") = 10, Rcpp::Named("dL") = model_parameters[1],
                                                                              Rcpp::Named("dP") = model_parameters[2], Rcpp::Named("muE0") = model_parameters[3],
                                                                              Rcpp::Named("muL0") = model_parameters[4], Rcpp::Named("muP") = model_parameters[5],
                                                                              Rcpp::Named("muM") = model_parameters[6], Rcpp::Named("lambda") = model_parameters[7],
                                                                              Rcpp::Named("tau") = model_parameters[8], Rcpp::Named("beta") = model_parameters[9],
                                                                              Rcpp::Named("overdisp") = model_parameters[10], Rcpp::Named("pop_frac") = model_parameters[11],
                                                                              Rcpp::Named("E") = 0.0, Rcpp::Named("L") = 0.0,
                                                                              Rcpp::Named("P") = 0.0, Rcpp::Named("M") = 0.0);

  arma::mat current_parameter_values(1, Rcpp::sum(fitted_yn));
  for (int u = 0; u < Rcpp::sum(fitted_yn); u++) {
    current_parameter_values(0, u) = MCMC_chain_output(0, u);
  }
  arma::mat proposed_parameter_values = mvrnormArma(current_parameter_values, current_covariance_matrix);

  // Filling input parameter vector with the output from the proposal function. Fixed parameters remain unchanged.
  for (int k = 0; k < indexer.size(); k++) {
    int index_for_adding = indexer[k];
    model_parameters_for_MCMC[index_for_adding] = proposed_parameter_values[k];
  }

  model_parameters_for_MCMC.names() = Rcpp::StringVector::create("dE", "dL", "dP", "muE0", "muL0", "muP", "muM", "lambda", "tau", "beta", "overdisp", "pop_frac", "E", "L", "P", "M");


  // Assessing the likelihood of the proposed parameter value compared to the current parameter value
  double proposed_posterior_likelihood = posterior_joint_proposals(N, rainfall, obsData, number_of_datapoints,
                                                                   data_timeframe, model_parameters_for_MCMC, static_parameters,
                                                                   density_function, prior_choice, fitted_yn);
  double likelihood_ratio = exp(proposed_posterior_likelihood - current_posterior_likelihood);

  // Deciding whether to accept or reject, updating chain and current_posterior_likelihood as appropriate
  if(R::runif(0, 1) < likelihood_ratio) {

    // Taking arma mat output proposed parameter values and sequentially adding the results to next row of
    // the MCMC chain (if these parameter values were accepted)
    for (int u = 0; u < proposed_parameter_values.size(); u++) {
      MCMC_chain_output(1, u) = proposed_parameter_values[u];
    }
    current_posterior_likelihood = proposed_posterior_likelihood;
    acceptances = acceptances + 1;
    accepted_checker[0] = 2;
  }
  else {

    // Taking arma mat output current parameter values and sequentially adding these same parameter values
    // to next row of the MCMC chain (as the proposed parameters were rejected)
    for (int u = 0; u < current_parameter_values.size(); u++) {
      MCMC_chain_output(1, u) = current_parameter_values[u];
    }
    current_posterior_likelihood = current_posterior_likelihood; // remove this as it's superfluous- just to tell myself what's going on
    rejections = rejections + 1;
    accepted_checker[0] = 1;
  }

  // Tracking the acceptance ratio
  current_acceptance_ratio = (acceptances / (1));
  acceptance_ratio_tracker[0] = current_acceptance_ratio;

  // Adapting the standard deviation of parameter proposals to ensure a decent acceptance ratio

    // Considered using current_parameter_values here but that isn't updated (if the proposed values are accepted)
    // until the start of the next iteration. The MCMC Output matrix IS however and so I use that to fill a temporary
    // vector containing what either are currently (if proposed were rejected) the parameter values OR what are
    // going to be assigned to current_parameter_values in the beginning of the next iteration (if proposed were accepted).
    arma::mat latest_parameter_values_but_col = MCMC_chain_output(0, Rcpp::_);
    arma::mat latest_parameter_values = latest_parameter_values_but_col.t();


    Rcpp::List adapter_output = joint_proposal_SD_adapter(1, 0, 1, current_scaling_factor, current_mu, latest_parameter_values, // technically parameter values for t + 1 as the acceptance/rejection step precedes calling this function
                                                          current_covariance_matrix);

    // // Currently having to take the output and assign it to Rcpp types before adding those elements to arma::mat/vec types.
    // // Ask Rich if there's a nicer way of doing this.
    // Rcpp::NumericVector new_mu = adapter_output["New_Mu"];
    // Rcpp::NumericMatrix new_covariance_matrix = adapter_output["New_Covariance_Matrix"];
    //
    // for (int u = 0; u < current_parameter_values.size(); u++) {
    //   current_mu[u] = new_mu[u];
    // }
    // for (int u = 0; u < current_parameter_values.size(); u++) {
    //   for (int v = 0; v < current_parameter_values.size(); v++) {
    //     current_covariance_matrix(u, v) = new_covariance_matrix(u, v);
    //   }
    // }
    // current_scaling_factor = adapter_output["New_Scaling_Factor"];


  return(Rcpp::List::create(Rcpp::Named("MCMC_Output") = MCMC_chain_output,
                            Rcpp::Named("Acceptance_Ratio_Final") = current_mu,
                            Rcpp::Named("CovMat") = initial_covariance_matrix_checker,
                            Rcpp::Named("ModParams")= model_parameters_for_MCMC,
                            Rcpp::Named("CurrParams") = current_parameter_values,
                            Rcpp::Named("Proposed") = latest_parameter_values));
}
