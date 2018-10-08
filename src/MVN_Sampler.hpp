#ifndef MVN_SAMPLER_H
#define MVN_SAMPLER_H

// Things that need to be included in the file mosquito_population_model.cpp
#include <RcppArmadillo.h>
#include <vector>
// [[Rcpp::depends(RcppArmadillo)]]


// Functions Contained Within File mosquito_population_model.cpp
arma::mat mvrnormArma(arma::mat mu, arma::mat sigma);


#endif
