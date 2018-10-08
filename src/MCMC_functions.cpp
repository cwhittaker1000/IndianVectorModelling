#include "Particle_Filter_Functions.hpp"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' @export
// [[Rcpp::export]]
double prior(Rcpp::NumericVector parameter_values, Rcpp::String prior_choice, Rcpp::LogicalVector fitted_yn, Rcpp::String likelihood_choice) {

  double res = 0;
  bool overdisp = true;

  if (likelihood_choice == "poisson") {
    overdisp = false;
  }

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
    if (overdisp == true) {
      res = res + R::dunif(parameter_values["overdisp"], 0.001, 1, TRUE); //overdispersion
    }
    if (fitted_yn["pop_frac"]) {
      res = res + R::dunif(parameter_values["pop_frac"], 0.001, 1, TRUE); //pop_frac
    }
    if (fitted_yn["scaling_factor"]) {
      res = res + R::dunif(parameter_values["scaling_factor"], 0.001, 20, TRUE); //pop_frac BETA DISTRIBUTION USED BASED ON MRR DATA
    }
    if (fitted_yn["z"]) {
      res = res + R::dunif(parameter_values["z"], 0.001, 1000, TRUE); //pop_frac BETA DISTRIBUTION USED BASED ON MRR DATA
    }
    if (fitted_yn["K_static"]) {
      res = res + R::dunif(parameter_values["K_static"], 0, 200, TRUE); //z Uniform distribution used
    }
    if (fitted_yn["K_max"]) {
      res = res + R::dunif(parameter_values["K_max"], 0, 500, TRUE); //z Uniform distribution used
    }
    if (fitted_yn["hill_1"]) {
      res = res + R::dunif(parameter_values["hill_1"], 0.0001, 200, TRUE); //z Uniform distribution used
    }
    if (fitted_yn["hill_2"]) {
      res = res + R::dunif(parameter_values["hill_2"], -200, 200, TRUE); //z Uniform distribution used
    }
  }

  if (prior_choice == "informative") {

    if (fitted_yn["dE"]) {
      res = res + R::dnorm4(parameter_values["dE"], 6.64, 0.933, TRUE); //dE - Michael's Posterior
    }
    if (fitted_yn["dL"]) {
      res = res + R::dnorm4(parameter_values["dL"], 3.72, 0.862, TRUE); //dL - Michael's Posterior
    }
    if (fitted_yn["dP"]) {
      res = res + R::dnorm4(parameter_values["dP"], 0.64, 0.291, TRUE); //dP - Michael's Posterior
    }
    if (fitted_yn["muE0"]) {
      res = res + R::dnorm4(parameter_values["muE0"], 0.034, 0.0056, TRUE); //muE0 - Michael's Posterior
    }
    if (fitted_yn["muL0"]) {
      res = res + R::dnorm4(parameter_values["muL0"], 0.035, 0.0056, TRUE); //muL0 - Michael's Posterior
    }
    if (fitted_yn["muP"]) {
      res = res + R::dnorm4(parameter_values["muP"], 0.25, 0.03571, TRUE); //muP - Michael's Posterior
    }
    if (fitted_yn["muM"]) {
      res = res + R::dnorm4(parameter_values["muM"], 0.096, 0.00459, TRUE); //muM - Michael's Posterior
    }
    if (fitted_yn["lambda"]) {
      res = res + R::dnorm4(parameter_values["lambda"], 13.25, 5, TRUE); //lambda - NOT Michael's Posterior
    }
    if (fitted_yn["tau"]) {
      res = res + R::dnorm4(parameter_values["tau"], 4, 0.7653, TRUE); //tau DOESN'T THIS HAVE TO BE DISCRETE REALLY? YES. CHANGE THIS. ALSO VALUES ARE FOR Michael's Posterior
    }
    if (fitted_yn["beta"]) {
      res = res + R::dnorm4(parameter_values["beta"], 21.19, 4.9082, TRUE); //beta - Michael's Posterior
    }
    if (overdisp == true) {
      res = res + R::dnorm4(parameter_values["overdisp"], 1, 0.5, TRUE); //overdispersion
    }
    if (fitted_yn["pop_frac"]) {
      res = res + R::dbeta(parameter_values["pop_frac"], 0.2750417, 5.8075736, TRUE); //pop_frac BETA DISTRIBUTION USED BASED ON MRR DATA
    }
    if (fitted_yn["scaling_factor"]) {
      res = res + R::dunif(parameter_values["scaling_factor"], 1, 500, TRUE); //pop_frac BETA DISTRIBUTION USED BASED ON MRR DATA
    }
    if (fitted_yn["z"]) {
      res = res + R::dunif(parameter_values["z"], 1, 1000, TRUE); //pop_frac BETA DISTRIBUTION USED BASED ON MRR DATA
    }
    if (fitted_yn["K_static"]) {
      res = res + R::dunif(parameter_values["K_static"], 0, 200, TRUE); //z Uniform distribution used
    }
    if (fitted_yn["K_max"]) {
      res = res + R::dunif(parameter_values["K_max"], 0, 500, TRUE); //z Uniform distribution used
    }
    if (fitted_yn["hill_1"]) {
      res = res + R::dunif(parameter_values["hill_1"], 0.0001, 200, TRUE); //z Uniform distribution used
    }
    if (fitted_yn["hill_2"]) {
      res = res + R::dunif(parameter_values["hill_2"], -200, 200, TRUE); //z Uniform distribution used
    }
  }

  return(res);
}

//' @export
// [[Rcpp::export]]
double likelihood_function(int N, std::vector <double> rainfall, std::vector <int> obsData,
                           Rcpp::NumericVector fitted_parameters, Rcpp::NumericVector static_parameters, Rcpp::String density_function,
                           double sampling_point, Rcpp::StringVector month_vector, Rcpp::String rainfall_relationship,
                           Rcpp::String rainfall_effect, Rcpp::String likelihood_choice) {

  double loglikelihood = min_output_particle_filter(N, rainfall, obsData, fitted_parameters,
                                                    static_parameters, density_function, sampling_point,
                                                    month_vector, rainfall_relationship, rainfall_effect, likelihood_choice);
  return(loglikelihood);

}

//' @export
// [[Rcpp::export]]
double posterior_joint_proposals(int N, std::vector <double> rainfall, std::vector <int> obsData,
                                 Rcpp::NumericVector fitted_parameters, Rcpp::NumericVector static_parameters,
                                 Rcpp::String density_function, Rcpp::String prior_choice, Rcpp::LogicalVector fitted_yn, double sampling_point,
                                 Rcpp::StringVector month_vector, Rcpp::String rainfall_relationship, Rcpp::String rainfall_effect,
                                 Rcpp::String likelihood_choice) {

  double posterior_likelihood;

  if (prior(fitted_parameters, prior_choice, fitted_yn, likelihood_choice) < -5000) {
    posterior_likelihood = -10000;
    return(posterior_likelihood);
  }

  else {
    posterior_likelihood = likelihood_function(N, rainfall, obsData, fitted_parameters, static_parameters,
                                               density_function, sampling_point, month_vector, rainfall_relationship,
                                               rainfall_effect, likelihood_choice)
    + prior(fitted_parameters, prior_choice, fitted_yn, likelihood_choice);
    return(posterior_likelihood);
  }

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






