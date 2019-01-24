#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

//' @export
// [[Rcpp::export]]
Rcpp::List kloop(double offset, double total_length_in_steps,
                 double temp_tau_static_prior, double temp_tau_rainfall_prior) {

  std::vector<double> Offset_Period_Exponential_Static;
  std::vector<double> Offset_Period_Exponential_Rainfall;
  std::vector<double> Tau_Static_Exponential_Normalisation_Factor;
  std::vector<double> Tau_Rainfall_Exponential_Normalisation_Factor;

  for (int s = 0; s <= (offset + total_length_in_steps); s++) {

    double calc_placeholder_static = -((offset + total_length_in_steps) - s)/temp_tau_static_prior;
    double calc_placeholder_rainfall = -((offset + total_length_in_steps) - s)/temp_tau_rainfall_prior;

    Offset_Period_Exponential_Static.insert(Offset_Period_Exponential_Static.begin(), exp(calc_placeholder_static));
    Offset_Period_Exponential_Rainfall.insert(Offset_Period_Exponential_Rainfall.begin(), exp(calc_placeholder_rainfall));

  }

  for (int x = 0; x <= (offset + total_length_in_steps); x++) {

    double tau_static_temp_exp_norm_factor = (1.0 / (temp_tau_static_prior * (1 - exp(- (x + 1) / temp_tau_static_prior))));
    double tau_rainfall_temp_exp_norm_factor = (1.0 / (temp_tau_rainfall_prior * (1 - exp(-(x + 1) / temp_tau_rainfall_prior))));

    Tau_Static_Exponential_Normalisation_Factor.emplace_back(tau_static_temp_exp_norm_factor);
    Tau_Rainfall_Exponential_Normalisation_Factor.emplace_back(tau_rainfall_temp_exp_norm_factor);

  }

  return Rcpp::List::create(Rcpp::Named("Offset_Period_Exponential_Static") = Offset_Period_Exponential_Static,
                            Rcpp::Named("Offset_Period_Exponential_Rainfall") = Offset_Period_Exponential_Rainfall,
                            Rcpp::Named("Tau_Static_Exponential_Normalisation_Factor") = Tau_Static_Exponential_Normalisation_Factor,
                            Rcpp::Named("Tau_Rainfall_Exponential_Normalisation_Factor") = Tau_Rainfall_Exponential_Normalisation_Factor);
}
