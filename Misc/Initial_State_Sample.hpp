#ifndef INITIAL_STATE_SAMPLE_H
#define INITIAL_STATE_SAMPLE_H

// Things that need to be included in the file initial_state_sample.cpp
#include <Rcpp.h>

// Functions Contained Within File initial_state_sample.cpp
std::vector <int> initial_state_sample(Rcpp::NumericVector fitted_parameters, Rcpp::NumericVector static_parameters, double initial_K);

#endif
