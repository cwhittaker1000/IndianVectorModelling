#ifndef MIN_OUTPUT_PARTICLE_FILTER_H
#define MIN_OUTPUT_PARTICLE_FILTER_H

// Things that need to be included in the file mosquito_population_model.cpp
#include <Rcpp.h>
#include <vector>


// Functions Contained Within File mosquito_population_model.cpp
double min_output_particle_filter(int N, std::vector <double> rainfall, std::vector <int> obsData, int number_of_datapoints, int data_timeframe, Rcpp::NumericVector fitted_parameters, Rcpp::NumericVector static_parameters, Rcpp::String density_function);


#endif
