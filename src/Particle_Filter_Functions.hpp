#ifndef PARTICLE_FILTER_FUNCTIONS_H
#define PARTICLE_FILTER_FUNCTIONS_H

// Things that need to be included in the file mosquito_population_model.cpp
#include <RcppArmadillo.h>
#include <vector>
#include <string>
#include <iostream>

// Functions Contained Within File mosquito_population_model.cpp
double min_output_particle_filter(int N, std::vector <double> rainfall, std::vector <int> obsData, int number_of_datapoints, int data_timeframe, Rcpp::NumericVector fitted_parameters, Rcpp::NumericVector static_parameters, Rcpp::String density_function);
Rcpp::List full_output_particle_filter(int N, std::vector <double> rainfall, std::vector <int> obsData, int number_of_datapoints, int data_timeframe, Rcpp::NumericVector fitted_parameters, Rcpp::NumericVector static_parameters, Rcpp::String density_function);
double Negative_Binomial(double k, double n, double r, double p);
std::vector <double> Particle_Weight_Normalisation(std::vector <double> particle_weights);
std::vector <int> initial_state_sample(Rcpp::NumericVector fitted_parameters, Rcpp::NumericVector static_parameters, double initial_K);
std::vector<int> weighted_sampling_with_replacement(int n, std::vector<double> prob, double p_tot, int K);


#endif
