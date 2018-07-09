#include <Rcpp.h>
#include "Min_Output_Particle_Filter.hpp"

// [[Rcpp::export]]
double prior(Rcpp::NumericVector parameter_values, Rcpp::String prior_choice) {

  double res = 0;

  if (prior_choice == "uninformative") {

    res = res + R::dunif(parameter_values["dE"], 0.0001, 1, TRUE); //dE
    res = res + R::dunif(parameter_values["dL"], 0.0001, 1, TRUE); //dL
    res = res + R::dunif(parameter_values["dP"], 0.0001, 3, TRUE); //dP
    res = res + R::dunif(parameter_values["muE0"], 0.0001, 1, TRUE); //muE0
    res = res + R::dunif(parameter_values["muL0"], 0.0001, 1, TRUE); //muL0
    res = res + R::dunif(parameter_values["muP"], 0.0001, 1, TRUE); //muP
    res = res + R::dunif(parameter_values["muM"], 0.0001, 1, TRUE); //muM
    res = res + R::dunif(parameter_values["lambda"], 1, 10, TRUE); //lambda
    res = res + R::dunif(parameter_values["tau"], 1, 30, TRUE); //tau DOESN'T THIS HAVE TO BE DISCRETE REALLY?
    res = res + R::dunif(parameter_values["beta"], 1, 20, TRUE); //beta
    res = res + R::dunif(parameter_values["overdisp"], 1, 20, TRUE); //overdispersion
    res = res + R::dunif(parameter_values["pop_frac"], 0.001, 1, TRUE); //pop_frac

  }

  if (prior_choice == "informative") {

    res = res + R::dnorm4(parameter_values["dE"], 0.150602, 0.04, TRUE); //dE
    res = res + R::dnorm4(parameter_values["dL"], 0.268812, 0.06, TRUE); //dL
    res = res + R::dnorm4(parameter_values["dP"], 1, 0.1, TRUE); //dP
    res = res + R::dnorm4(parameter_values["muE0"], 0.035, 0.0056, TRUE); //muE0
    res = res + R::dnorm4(parameter_values["muL0"], 0.035, 0.0056, TRUE); //muL0
    res = res + R::dnorm4(parameter_values["muP"], 0.25, 0.05, TRUE); //muP
    res = res + R::dnorm4(parameter_values["muM"], 0.091, 0.005, TRUE); //muM
    res = res + R::dnorm4(parameter_values["lambda"], 13.06, 4, TRUE); //lambda
    res = res + R::dnorm4(parameter_values["tau"], 7, 3, TRUE); //tau DOESN'T THIS HAVE TO BE DISCRETE REALLY?
    res = res + R::dnorm4(parameter_values["beta"], 15, 7, TRUE); //beta
    res = res + R::dunif(parameter_values["overdisp"], 0.001, 10, TRUE); //overdispersion DECIDE WHETHER TO HAVE THIS NORMAL DIST
    res = res + R::dunif(parameter_values["pop_frac"], 0.001, 1, TRUE); //pop_frac DECIDE WHETHER TO HAVE THIS NORMAL DIST

  }

  return(res);
}

// [[Rcpp::export]]
double likelihood_function(int N, std::vector <double> rainfall, std::vector <int> obsData, int number_of_datapoints, int data_timeframe, Rcpp::NumericVector fitted_parameters, Rcpp::NumericVector static_parameters, Rcpp::String density_function) {

  double loglikelihood = min_output_particle_filter(N, rainfall, obsData, number_of_datapoints, data_timeframe, fitted_parameters, static_parameters, density_function);
  return(loglikelihood);

}

// [[Rcpp::export]]
double posterior(int N, std::vector <double> rainfall, std::vector <int> obsData, int number_of_datapoints,
                 int data_timeframe, Rcpp::NumericVector fitted_parameters, Rcpp::NumericVector static_parameters,
                 Rcpp::String density_function, double parameter_value, Rcpp::String parameter_index, Rcpp::String prior_choice) {

  double posterior_likelihood;
  fitted_parameters[parameter_index] = parameter_value;

  if (prior(fitted_parameters, prior_choice) < -5000) {
    posterior_likelihood = -10000;
    return(posterior_likelihood);
  }

  else {
    posterior_likelihood = likelihood_function(N, rainfall, obsData, number_of_datapoints, data_timeframe, fitted_parameters, static_parameters, density_function)
    + prior(fitted_parameters, prior_choice);
    return(posterior_likelihood);
  }

}
