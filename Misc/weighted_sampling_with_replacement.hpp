#ifndef WEIGHTED_SAMPLING_WITH_REPLACEMENT_H
#define WEIGHTED_SAMPLING_WITH_REPLACEMENT_H

// Things that need to be included in the file weight_sampling_with_replacement.cpp
#include <Rcpp.h>
#include <vector>

// Functions Contained Within File weight_sampling_with_replacement.cpp
std::vector<int> weighted_sampling(int n, std::vector<double> prob, double p_tot, int K);

#endif
