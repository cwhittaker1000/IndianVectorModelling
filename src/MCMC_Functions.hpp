#ifndef MCMC_FUNCTIONS_H
#define MCMC_FUNCTIONS_H

// Things that need to be included in the file mosquito_population_model.cpp
#include <RcppArmadillo.h>
#include <vector>

// Functions Contained Within File mosquito_population_model.cpp
double prior(Rcpp::NumericVector parameter_values, Rcpp::LogicalVector fitted_yn, Rcpp::String likelihood_choice);

double likelihood(int N, std::vector <int> obsData,
                  Rcpp::NumericVector fitted_parameters, Rcpp::NumericVector static_parameters,
                  std::vector <double> rainfall,
                  Rcpp::String mortality_density_function, Rcpp::String rainfall_relationship, Rcpp::String rainfall_effect, Rcpp::String decline_type,
                  double sampling_point, Rcpp::StringVector offset_month_vector, Rcpp::StringVector sampling_month_vector,
                  Rcpp::String likelihood_choice, Rcpp::String calc_inside_mosquito_model);

double posterior(int N, std::vector <int> obsData,
                 Rcpp::NumericVector fitted_parameters, Rcpp::NumericVector static_parameters,
                 std::vector <double> rainfall,
                 Rcpp::LogicalVector fitted_yn,
                 Rcpp::String mortality_density_function, Rcpp::String rainfall_relationship, Rcpp::String rainfall_effect, Rcpp::String decline_type,
                 double sampling_point, Rcpp::StringVector offset_month_vector, Rcpp::StringVector sampling_month_vector,
                 Rcpp::String likelihood_choice, Rcpp::String calc_inside_mosquito_model);

Rcpp::List proposal_SD_adapter(double accepted_variable, double current_iteration, double iteration_cooling_began,
                               double current_scaling_factor, arma::mat mu_previous, arma::mat current_parameter_values, // technically parameter values for t + 1 as the acceptance/rejection step precedes calling this function
                              arma::mat current_covariance_matrix);

arma::mat mvrnormArma(arma::mat mu, arma::mat sigma);

#endif
