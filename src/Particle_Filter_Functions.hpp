#ifndef PARTICLE_FILTER_FUNCTIONS_H
#define PARTICLE_FILTER_FUNCTIONS_H

// Things that need to be included in the file mosquito_population_model.cpp
#include <RcppArmadillo.h>
#include <vector>
#include <string>
#include <iostream>

// Functions Contained Within File mosquito_population_model.cpp
double min_output_particle_filter(int N, std::vector <double> rainfall, std::vector <int> obsData,
                                  Rcpp::NumericVector fitted_parameters, Rcpp::NumericVector static_parameters,
                                  Rcpp::String mortality_density_function, Rcpp::String rainfall_relationship,
                                  Rcpp::String rainfall_effect, Rcpp::String decline_type,
                                  double sampling_point, Rcpp::StringVector offset_month_vector, Rcpp::StringVector sampling_month_vector, Rcpp::String likelihood_choice);

Rcpp::List full_output_particle_filter(int N, std::vector <double> rainfall, std::vector <int> obsData,
                                       Rcpp::NumericVector fitted_parameters, Rcpp::NumericVector static_parameters,
                                       Rcpp::String mortality_density_function, Rcpp::String rainfall_relationship,
                                       Rcpp::String rainfall_effect, Rcpp::String decline_type,
                                       double sampling_point, Rcpp::StringVector offset_month_vector, Rcpp::StringVector sampling_month_vector, Rcpp::String likelihood_choice);

double Negative_Binomial(double k, double n, double r, double p);

double Poisson(double observed_data, double model_output, double population_fraction);

std::vector <double> Particle_Weight_Normalisation(std::vector <double> particle_weights);

std::vector <int> initial_state_sample(Rcpp::NumericVector fitted_parameters, Rcpp::NumericVector static_parameters, double initial_K);

std::vector<int> weighted_sampling_with_replacement(int n, std::vector<double> prob, double p_tot, int K);


#endif
