#ifndef PROPOSAL_FUNCTIONS_H
#define PROPOSAL_FUNCTIONS_H

// Things that need to be included in the file mosquito_population_model.cpp
#include <Rcpp.h>
#include <vector>


// Functions Contained Within File mosquito_population_model.cpp
double proposal_function(double sd, double current_parameter_value);
double proposal_SD_adapter(double current_sd, double acceptance_ratio, double current_acceptance_ratio, double max_sd);
double proposal_SD_adapter_mark_two(double current_sd, double current_acceptance_ratio, double max_sd, double cooling_factor,
                                    double current_iteration, double iteration_cooling_began);



#endif
