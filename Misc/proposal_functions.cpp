#include "Proposal_Functions.hpp"

//' @export
// [[Rcpp::export]]
double proposal_function(double sd, double current_parameter_value) {
  double ran = R::rnorm(0.0, sd);
  double proposed_value = ran + current_parameter_value;
  if (proposed_value <= 0) {
    proposed_value = 0;
  }
  return proposed_value;
}

//' @export
// [[Rcpp::export]]
double proposal_SD_adapter(double current_sd, double acceptance_ratio, double current_acceptance_ratio, double max_sd){


  if (current_acceptance_ratio == 1) {
    current_acceptance_ratio = 0.99;
  }
  if (current_acceptance_ratio == 0) {
    current_acceptance_ratio = 0.01;
  }

  double res = (current_sd * R::qnorm5(acceptance_ratio/2, 0, 1, TRUE, FALSE)) / R::qnorm5(current_acceptance_ratio/2, 0, 1, TRUE, FALSE);

  if (res > max_sd) {
    res = max_sd;
  }

  if (res <= 0) {
    res = 0.001;
  }

  return res;
}


//' @export
// [[Rcpp::export]]
double proposal_SD_adapter_mark_two(double current_sd, double current_acceptance_ratio, double max_sd, double cooling_factor,
                                    double current_iteration, double iteration_cooling_began) {

  double iterations_since_cooling_began = current_iteration - iteration_cooling_began;

  if (current_acceptance_ratio == 1) {
    current_acceptance_ratio = 0.99;
  }
  if (current_acceptance_ratio == 0) {
    current_acceptance_ratio = 0.01;
  }

  double res = current_sd * exp(pow(cooling_factor, iterations_since_cooling_began) * (current_acceptance_ratio - 0.234));

  if (res > max_sd) {
    res = max_sd;
  }

  if (res <= 0) {
    res = 0.001;
  }

  return res;
}


//' @export
// [[Rcpp::export]]
double proposal_SD_adapter_mark_three(double accepted_variable, double current_iteration, double iteration_cooling_began, double current_sd) {

  double iterations_since_cooling_began = current_iteration - iteration_cooling_began;

  double scaling_factor = pow((iterations_since_cooling_began + 1), -0.6);

  double log_updated_sd = log(current_sd) + scaling_factor * (accepted_variable - 0.25);

  double updated_sd = exp(log_updated_sd);

  return(updated_sd);

}
