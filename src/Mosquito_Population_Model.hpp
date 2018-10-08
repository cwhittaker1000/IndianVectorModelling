#ifndef MOSQUITO_POPULATION_MODEL_H
#define MOSQUITO_POPULATION_MODEL_H

// Things that need to be included in the file mosquito_population_model.cpp
#include <RcppArmadillo.h>
#include <vector>
#include <map>
#include <iostream>
#include <random>
#include <string>
#include <fstream>
#include <numeric>
#include <boost/regex.hpp>
#include <boost/math/distributions/uniform.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/random/binomial_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/math/distributions/normal.hpp>
#include <stdio.h>
#include <stdlib.h>

// Functions Contained Within File mosquito_population_model.cpp
Rcpp::List mosquito_population_model(int start_time, int end, Rcpp::NumericVector fitted_parameters,
                                     Rcpp::NumericVector static_parameters, std::vector<double> rainfall,
                                     Rcpp::String density_function, Rcpp::String rainfall_relationship, Rcpp::String rainfall_effect);


#endif
