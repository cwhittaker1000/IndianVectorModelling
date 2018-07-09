#ifndef MCMC_FUNCTIONS_H
#define MCMC_FUNCTIONS_H

// Things that need to be included in the file mosquito_population_model.cpp
#include <Rcpp.h>
#include <vector>

// Functions Contained Within File mosquito_population_model.cpp
double prior(Rcpp::NumericVector parameter_values,  Rcpp::String prior_choice);
double likelihood_function(int N, std::vector <double> rainfall, std::vector <int> obsData, int number_of_datapoints, int data_timeframe, Rcpp::NumericVector fitted_parameters, Rcpp::NumericVector static_parameters, Rcpp::String density_function);
double posterior(int N, std::vector <double> rainfall, std::vector <int> obsData, int number_of_datapoints,
                 int data_timeframe, Rcpp::NumericVector fitted_parameters, Rcpp::NumericVector static_parameters,
                 Rcpp::String density_function, double parameter_value, Rcpp::String parameter_index,  Rcpp::String prior_choice);



#endif
