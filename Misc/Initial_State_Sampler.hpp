#ifndef INITIAL_STATE_SAMPLER_H
#define INITIAL_STATE_SAMPLER_H


#include "Mosquito_Model.hpp"
//#include "Particle_Filter_Functions.hpp"

// Functions Contained Within File mosquito_population_model.cpp
std::vector <int> initial_state_sample(Rcpp::NumericVector fitted_parameters,
                                       Rcpp::NumericVector static_parameters,
                                       Rcpp::String mortality_density_function);

#endif
