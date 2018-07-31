#include "Particle_Filter_Functions.hpp"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' @export
// [[Rcpp::export]]
double prior_seq_proposals(Rcpp::NumericVector parameter_values, Rcpp::String prior_choice, Rcpp::LogicalVector fitted_yn) {

  double res = 0;

  if (prior_choice == "uninformative") {

    if (fitted_yn["dE"]) {
      res = res + R::dunif(parameter_values["dE"], 0.0001, 1, TRUE); //dE
    }
    if (fitted_yn["dL"]) {
      res = res + R::dunif(parameter_values["dE"], 0.0001, 1, TRUE); //dE
    }
    if (fitted_yn["dP"]) {
      res = res + R::dunif(parameter_values["dP"], 0.0001, 3, TRUE); //dP
    }
    if (fitted_yn["muE0"]) {
      res = res + R::dunif(parameter_values["muE0"], 0.0001, 1, TRUE); //muE0
    }
    if (fitted_yn["muL0"]) {
      res = res + R::dunif(parameter_values["muL0"], 0.0001, 1, TRUE); //muL0
    }
    if (fitted_yn["muP"]) {
      res = res + R::dunif(parameter_values["muP"], 0.0001, 1, TRUE); //muP
    }
    if (fitted_yn["muM"]) {
      res = res + R::dunif(parameter_values["muM"], 0.0001, 1, TRUE); //muM
    }
    if (fitted_yn["lambda"]) {
      res = res + R::dunif(parameter_values["lambda"], 1, 10, TRUE); //lambda
    }
    if (fitted_yn["tau"]) {
      res = res + R::dunif(parameter_values["tau"], 1, 30, TRUE); //tau DOESN'T THIS HAVE TO BE DISCRETE REALLY?
    }
    if (fitted_yn["beta"]) {
      res = res + R::dunif(parameter_values["beta"], 1, 20, TRUE); //beta
    }
    if (fitted_yn["overdisp"]) {
      res = res + R::dunif(parameter_values["overdisp"], 0.001, 1, TRUE); //overdispersion
    }
    if (fitted_yn["pop_frac"]) {
      res = res + R::dunif(parameter_values["pop_frac"], 0.001, 1, TRUE); //pop_frac
    }
  }

  if (prior_choice == "informative") {

    if (fitted_yn["dE"]) {
      res = res + R::dnorm4(parameter_values["dE"], 0.150602, 0.04, TRUE); //dE
    }
    if (fitted_yn["dL"]) {
      res = res + R::dnorm4(parameter_values["dL"], 0.268812, 0.06, TRUE); //dL
    }
    if (fitted_yn["dP"]) {
      res = res + R::dnorm4(parameter_values["dP"], 1, 0.1, TRUE); //dP
    }
    if (fitted_yn["muE0"]) {
      res = res + R::dnorm4(parameter_values["muE0"], 0.035, 0.0056, TRUE); //muE0
    }
    if (fitted_yn["muL0"]) {
      res = res + R::dnorm4(parameter_values["muL0"], 0.035, 0.0056, TRUE); //muL0
    }
    if (fitted_yn["muP"]) {
      res = res + R::dnorm4(parameter_values["muP"], 0.25, 0.05, TRUE); //muP
    }
    if (fitted_yn["muM"]) {
      res = res + R::dnorm4(parameter_values["muM"], 0.091, 0.005, TRUE); //muM
    }
    if (fitted_yn["lambda"]) {
      res = res + R::dnorm4(parameter_values["lambda"], 13.06, 4, TRUE); //lambda
    }
    if (fitted_yn["tau"]) {
      res = res + R::dnorm4(parameter_values["tau"], 7, 3, TRUE); //tau DOESN'T THIS HAVE TO BE DISCRETE REALLY?
    }
    if (fitted_yn["beta"]) {
      res = res + R::dnorm4(parameter_values["beta"], 15, 7, TRUE); //beta
    }
    if (fitted_yn["overdisp"]) {
      res = res + R::dunif(parameter_values["overdisp"], 0.001, 10, TRUE); //overdispersion DECIDE WHETHER TO HAVE THIS NORMAL DIST
    }
    if (fitted_yn["pop_frac"]) {
      res = res + R::dunif(parameter_values["pop_frac"], 0.001, 1, TRUE); //pop_frac DECIDE WHETHER TO HAVE THIS NORMAL DIST
    }
  }

  return(res);
}

//' @export
// [[Rcpp::export]]
double likelihood_function(int N, std::vector <double> rainfall, std::vector <int> obsData, int number_of_datapoints, int data_timeframe, Rcpp::NumericVector fitted_parameters, Rcpp::NumericVector static_parameters, Rcpp::String density_function) {

  double loglikelihood = min_output_particle_filter(N, rainfall, obsData, number_of_datapoints, data_timeframe, fitted_parameters, static_parameters, density_function);
  return(loglikelihood);

}

//' @export
// [[Rcpp::export]]
double posterior_seq_proposals(int N, std::vector <double> rainfall, std::vector <int> obsData, int number_of_datapoints,
                               int data_timeframe, Rcpp::NumericVector fitted_parameters, Rcpp::NumericVector static_parameters,
                               Rcpp::String density_function, double parameter_value, int parameter_index, Rcpp::String prior_choice, Rcpp::LogicalVector fitted_yn) {

  double posterior_likelihood;
  fitted_parameters[parameter_index] = parameter_value;

  if (prior_seq_proposals(fitted_parameters, prior_choice, fitted_yn) < -5000) {
    posterior_likelihood = -10000;
    return(posterior_likelihood);
  }

  else {
    posterior_likelihood = likelihood_function(N, rainfall, obsData, number_of_datapoints, data_timeframe, fitted_parameters, static_parameters, density_function)
    + prior_seq_proposals(fitted_parameters, prior_choice, fitted_yn);
    return(posterior_likelihood);
  }

}


//' @export
// [[Rcpp::export]]
double seq_proposal_function(double sd, double current_parameter_value) {
  double ran = R::rnorm(0.0, sd);
  double proposed_value = ran + current_parameter_value;
  if (proposed_value <= 0) {
    proposed_value = 0;
  }
  return proposed_value;
}


//' @export
// [[Rcpp::export]]
double seq_proposal_SD_adapter(double accepted_variable, double current_iteration, double iteration_cooling_began, double current_sd) {

  double iterations_since_cooling_began = current_iteration - iteration_cooling_began;
  double scaling_factor = pow((iterations_since_cooling_began + 1), -0.6);
  double log_updated_sd = log(current_sd) + scaling_factor * (accepted_variable - 0.25);
  double updated_sd = exp(log_updated_sd);

  return(updated_sd);

}


// BOTH MU_PREVIOUS AND CURRENT_PARAMETERS NEED TO BE ROW VECTORS

//' @export
// [[Rcpp::export]]
Rcpp::List joint_proposal_SD_adapter(double accepted_variable, double current_iteration, double iteration_cooling_began,
                                 double current_scaling_factor, arma::mat mu_previous, arma::mat current_parameter_values, // technically parameter values for t + 1 as the acceptance/rejection step precedes calling this function
                                 arma::mat current_covariance_matrix) {

  double iterations_since_cooling_began = current_iteration - iteration_cooling_began;
  double cooldown = pow((iterations_since_cooling_began + 1), -0.6);

  arma::mat new_correlation_matrix = ((1 - cooldown) * current_covariance_matrix) + (cooldown * (((current_parameter_values - mu_previous).t()) * (current_parameter_values - mu_previous)));
  arma::mat new_mu = ((1 - cooldown) * mu_previous) + (cooldown * current_parameter_values);
  double log_new_scaling_factor = log(current_scaling_factor) + cooldown * (accepted_variable - 0.25);
  double new_scaling_factor = exp(log_new_scaling_factor);
  arma::mat new_covariance_matrix = new_scaling_factor * new_correlation_matrix;

  return(Rcpp::List::create(Rcpp::Named("New_Covariance_Matrix") = new_covariance_matrix,
                            Rcpp::Named("New_Mu") = new_mu,
                            Rcpp::Named("New_Scaling_Factor") = new_scaling_factor,
                            Rcpp::Named("Cooldown") = cooldown));

}



//' @export
// [[Rcpp::export]]
double posterior_joint_proposals(int N, std::vector <double> rainfall, std::vector <int> obsData, int number_of_datapoints,
                                 int data_timeframe, Rcpp::NumericVector fitted_parameters, Rcpp::NumericVector static_parameters,
                                 Rcpp::String density_function, Rcpp::String prior_choice, Rcpp::LogicalVector fitted_yn) {

  double posterior_likelihood;

  if (prior_seq_proposals(fitted_parameters, prior_choice, fitted_yn) < -5000) {
    posterior_likelihood = -10000;
    return(posterior_likelihood);
  }

  else {
    posterior_likelihood = likelihood_function(N, rainfall, obsData, number_of_datapoints, data_timeframe, fitted_parameters, static_parameters, density_function)
    + prior_seq_proposals(fitted_parameters, prior_choice, fitted_yn);
    return(posterior_likelihood);
  }

}


