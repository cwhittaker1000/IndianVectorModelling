#ifndef MCMC_FUNCTIONS_H
#define MCMC_FUNCTIONS_H

// Things that need to be included in the file mosquito_population_model.cpp
#include <RcppArmadillo.h>
#include <vector>

// Functions Contained Within File mosquito_population_model.cpp
double prior_seq_proposals(Rcpp::NumericVector parameter_values,  Rcpp::String prior_choice, Rcpp::LogicalVector fitted_yn);
double likelihood_function(int N, std::vector <double> rainfall, std::vector <int> obsData, int number_of_datapoints, int data_timeframe, Rcpp::NumericVector fitted_parameters, Rcpp::NumericVector static_parameters, Rcpp::String density_function);
double posterior_seq_proposals(int N, std::vector <double> rainfall, std::vector <int> obsData, int number_of_datapoints,
                               int data_timeframe, Rcpp::NumericVector fitted_parameters, Rcpp::NumericVector static_parameters,
                               Rcpp::String density_function, double parameter_value, int parameter_index,  Rcpp::String prior_choice, Rcpp::LogicalVector fitted_yn);
double seq_proposal_SD_adapter(double accepted_variable, double current_iteration, double iteration_cooling_began,
                               double current_sd);
double seq_proposal_function(double sd, double current_parameter_value);
Rcpp::List joint_proposal_SD_adapter(double accepted_variable, double current_iteration, double iteration_cooling_began,
                                     double current_scaling_factor, arma::mat mu_previous, arma::mat current_parameter_values, // technically parameter values for t + 1 as the acceptance/rejection step precedes calling this function
                                     arma::mat current_covariance_matrix);
double posterior_joint_proposals(int N, std::vector <double> rainfall, std::vector <int> obsData, int number_of_datapoints,
                                 int data_timeframe, Rcpp::NumericVector fitted_parameters, Rcpp::NumericVector static_parameters,
                                 Rcpp::String density_function, Rcpp::String prior_choice, Rcpp::LogicalVector fitted_yn);


#endif
