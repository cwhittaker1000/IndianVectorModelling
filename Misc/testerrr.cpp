#include <Rcpp.h>
using namespace Rcpp;

//' @export
// [[Rcpp::export]]
Rcpp::List test_poo_runMCMC(int start_sd_adaptation, // Time to start SD adaptation - relates to SD adaptation
                       int end_sd_adaptation,
                       std::vector <double> max_sd, // Maximum SD allowed for parameter - relates to SD adaptation
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

  return(Rcpp::List::create(Rcpp::Named("MCMC_Output") = MCMC_chain_output,
                            Rcpp::Named("Proposal_Sd_final") = sd_proposals,
                            Rcpp::Named("naming") = naming,
                            Rcpp::Named("indexer") = indexer));
}
