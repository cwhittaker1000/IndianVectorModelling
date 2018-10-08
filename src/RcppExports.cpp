// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// mosquito_population_model_fluv
Rcpp::List mosquito_population_model_fluv(int start_time, int end, Rcpp::NumericVector fitted_parameters, Rcpp::NumericVector static_parameters, std::vector<double> rainfall, Rcpp::String mortality_density_function, Rcpp::String rainfall_relationship, Rcpp::String rainfall_effect);
RcppExport SEXP _IndianVectorModelling_mosquito_population_model_fluv(SEXP start_timeSEXP, SEXP endSEXP, SEXP fitted_parametersSEXP, SEXP static_parametersSEXP, SEXP rainfallSEXP, SEXP mortality_density_functionSEXP, SEXP rainfall_relationshipSEXP, SEXP rainfall_effectSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type start_time(start_timeSEXP);
    Rcpp::traits::input_parameter< int >::type end(endSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type fitted_parameters(fitted_parametersSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type static_parameters(static_parametersSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type rainfall(rainfallSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type mortality_density_function(mortality_density_functionSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type rainfall_relationship(rainfall_relationshipSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type rainfall_effect(rainfall_effectSEXP);
    rcpp_result_gen = Rcpp::wrap(mosquito_population_model_fluv(start_time, end, fitted_parameters, static_parameters, rainfall, mortality_density_function, rainfall_relationship, rainfall_effect));
    return rcpp_result_gen;
END_RCPP
}
// Hill_Function
double Hill_Function(double rainfall, double K_max, double a, double b);
RcppExport SEXP _IndianVectorModelling_Hill_Function(SEXP rainfallSEXP, SEXP K_maxSEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type rainfall(rainfallSEXP);
    Rcpp::traits::input_parameter< double >::type K_max(K_maxSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(Hill_Function(rainfall, K_max, a, b));
    return rcpp_result_gen;
END_RCPP
}
// test_initial_state_sample
std::vector <int> test_initial_state_sample(Rcpp::NumericVector fitted_parameters, Rcpp::NumericVector static_parameters, Rcpp::String mortality_density_function, double initial_K);
RcppExport SEXP _IndianVectorModelling_test_initial_state_sample(SEXP fitted_parametersSEXP, SEXP static_parametersSEXP, SEXP mortality_density_functionSEXP, SEXP initial_KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type fitted_parameters(fitted_parametersSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type static_parameters(static_parametersSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type mortality_density_function(mortality_density_functionSEXP);
    Rcpp::traits::input_parameter< double >::type initial_K(initial_KSEXP);
    rcpp_result_gen = Rcpp::wrap(test_initial_state_sample(fitted_parameters, static_parameters, mortality_density_function, initial_K));
    return rcpp_result_gen;
END_RCPP
}
// prior
double prior(Rcpp::NumericVector parameter_values, Rcpp::String prior_choice, Rcpp::LogicalVector fitted_yn, Rcpp::String likelihood_choice);
RcppExport SEXP _IndianVectorModelling_prior(SEXP parameter_valuesSEXP, SEXP prior_choiceSEXP, SEXP fitted_ynSEXP, SEXP likelihood_choiceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type parameter_values(parameter_valuesSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type prior_choice(prior_choiceSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type fitted_yn(fitted_ynSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type likelihood_choice(likelihood_choiceSEXP);
    rcpp_result_gen = Rcpp::wrap(prior(parameter_values, prior_choice, fitted_yn, likelihood_choice));
    return rcpp_result_gen;
END_RCPP
}
// likelihood_function
double likelihood_function(int N, std::vector <double> rainfall, std::vector <int> obsData, Rcpp::NumericVector fitted_parameters, Rcpp::NumericVector static_parameters, Rcpp::String density_function, double sampling_point, Rcpp::StringVector month_vector, Rcpp::String rainfall_relationship, Rcpp::String rainfall_effect, Rcpp::String likelihood_choice);
RcppExport SEXP _IndianVectorModelling_likelihood_function(SEXP NSEXP, SEXP rainfallSEXP, SEXP obsDataSEXP, SEXP fitted_parametersSEXP, SEXP static_parametersSEXP, SEXP density_functionSEXP, SEXP sampling_pointSEXP, SEXP month_vectorSEXP, SEXP rainfall_relationshipSEXP, SEXP rainfall_effectSEXP, SEXP likelihood_choiceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< std::vector <double> >::type rainfall(rainfallSEXP);
    Rcpp::traits::input_parameter< std::vector <int> >::type obsData(obsDataSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type fitted_parameters(fitted_parametersSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type static_parameters(static_parametersSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type density_function(density_functionSEXP);
    Rcpp::traits::input_parameter< double >::type sampling_point(sampling_pointSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type month_vector(month_vectorSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type rainfall_relationship(rainfall_relationshipSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type rainfall_effect(rainfall_effectSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type likelihood_choice(likelihood_choiceSEXP);
    rcpp_result_gen = Rcpp::wrap(likelihood_function(N, rainfall, obsData, fitted_parameters, static_parameters, density_function, sampling_point, month_vector, rainfall_relationship, rainfall_effect, likelihood_choice));
    return rcpp_result_gen;
END_RCPP
}
// posterior_joint_proposals
double posterior_joint_proposals(int N, std::vector <double> rainfall, std::vector <int> obsData, Rcpp::NumericVector fitted_parameters, Rcpp::NumericVector static_parameters, Rcpp::String density_function, Rcpp::String prior_choice, Rcpp::LogicalVector fitted_yn, double sampling_point, Rcpp::StringVector month_vector, Rcpp::String rainfall_relationship, Rcpp::String rainfall_effect, Rcpp::String likelihood_choice);
RcppExport SEXP _IndianVectorModelling_posterior_joint_proposals(SEXP NSEXP, SEXP rainfallSEXP, SEXP obsDataSEXP, SEXP fitted_parametersSEXP, SEXP static_parametersSEXP, SEXP density_functionSEXP, SEXP prior_choiceSEXP, SEXP fitted_ynSEXP, SEXP sampling_pointSEXP, SEXP month_vectorSEXP, SEXP rainfall_relationshipSEXP, SEXP rainfall_effectSEXP, SEXP likelihood_choiceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< std::vector <double> >::type rainfall(rainfallSEXP);
    Rcpp::traits::input_parameter< std::vector <int> >::type obsData(obsDataSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type fitted_parameters(fitted_parametersSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type static_parameters(static_parametersSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type density_function(density_functionSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type prior_choice(prior_choiceSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type fitted_yn(fitted_ynSEXP);
    Rcpp::traits::input_parameter< double >::type sampling_point(sampling_pointSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type month_vector(month_vectorSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type rainfall_relationship(rainfall_relationshipSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type rainfall_effect(rainfall_effectSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type likelihood_choice(likelihood_choiceSEXP);
    rcpp_result_gen = Rcpp::wrap(posterior_joint_proposals(N, rainfall, obsData, fitted_parameters, static_parameters, density_function, prior_choice, fitted_yn, sampling_point, month_vector, rainfall_relationship, rainfall_effect, likelihood_choice));
    return rcpp_result_gen;
END_RCPP
}
// joint_proposal_SD_adapter
Rcpp::List joint_proposal_SD_adapter(double accepted_variable, double current_iteration, double iteration_cooling_began, double current_scaling_factor, arma::mat mu_previous, arma::mat current_parameter_values, arma::mat current_covariance_matrix);
RcppExport SEXP _IndianVectorModelling_joint_proposal_SD_adapter(SEXP accepted_variableSEXP, SEXP current_iterationSEXP, SEXP iteration_cooling_beganSEXP, SEXP current_scaling_factorSEXP, SEXP mu_previousSEXP, SEXP current_parameter_valuesSEXP, SEXP current_covariance_matrixSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type accepted_variable(accepted_variableSEXP);
    Rcpp::traits::input_parameter< double >::type current_iteration(current_iterationSEXP);
    Rcpp::traits::input_parameter< double >::type iteration_cooling_began(iteration_cooling_beganSEXP);
    Rcpp::traits::input_parameter< double >::type current_scaling_factor(current_scaling_factorSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mu_previous(mu_previousSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type current_parameter_values(current_parameter_valuesSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type current_covariance_matrix(current_covariance_matrixSEXP);
    rcpp_result_gen = Rcpp::wrap(joint_proposal_SD_adapter(accepted_variable, current_iteration, iteration_cooling_began, current_scaling_factor, mu_previous, current_parameter_values, current_covariance_matrix));
    return rcpp_result_gen;
END_RCPP
}
// mosquito_population_model
Rcpp::List mosquito_population_model(int start_time, int end, Rcpp::NumericVector fitted_parameters, Rcpp::NumericVector static_parameters, std::vector<double> rainfall, Rcpp::String mortality_density_function, Rcpp::String rainfall_relationship, Rcpp::String rainfall_effect);
RcppExport SEXP _IndianVectorModelling_mosquito_population_model(SEXP start_timeSEXP, SEXP endSEXP, SEXP fitted_parametersSEXP, SEXP static_parametersSEXP, SEXP rainfallSEXP, SEXP mortality_density_functionSEXP, SEXP rainfall_relationshipSEXP, SEXP rainfall_effectSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type start_time(start_timeSEXP);
    Rcpp::traits::input_parameter< int >::type end(endSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type fitted_parameters(fitted_parametersSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type static_parameters(static_parametersSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type rainfall(rainfallSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type mortality_density_function(mortality_density_functionSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type rainfall_relationship(rainfall_relationshipSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type rainfall_effect(rainfall_effectSEXP);
    rcpp_result_gen = Rcpp::wrap(mosquito_population_model(start_time, end, fitted_parameters, static_parameters, rainfall, mortality_density_function, rainfall_relationship, rainfall_effect));
    return rcpp_result_gen;
END_RCPP
}
// mvrnormArma
arma::mat mvrnormArma(arma::mat mu, arma::mat sigma);
RcppExport SEXP _IndianVectorModelling_mvrnormArma(SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(mvrnormArma(mu, sigma));
    return rcpp_result_gen;
END_RCPP
}
// min_output_particle_filter
double min_output_particle_filter(int N, std::vector <double> rainfall, std::vector <int> obsData, Rcpp::NumericVector fitted_parameters, Rcpp::NumericVector static_parameters, Rcpp::String density_function, double sampling_point, Rcpp::StringVector month_vector, Rcpp::String rainfall_relationship, Rcpp::String rainfall_effect, Rcpp::String likelihood_choice);
RcppExport SEXP _IndianVectorModelling_min_output_particle_filter(SEXP NSEXP, SEXP rainfallSEXP, SEXP obsDataSEXP, SEXP fitted_parametersSEXP, SEXP static_parametersSEXP, SEXP density_functionSEXP, SEXP sampling_pointSEXP, SEXP month_vectorSEXP, SEXP rainfall_relationshipSEXP, SEXP rainfall_effectSEXP, SEXP likelihood_choiceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< std::vector <double> >::type rainfall(rainfallSEXP);
    Rcpp::traits::input_parameter< std::vector <int> >::type obsData(obsDataSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type fitted_parameters(fitted_parametersSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type static_parameters(static_parametersSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type density_function(density_functionSEXP);
    Rcpp::traits::input_parameter< double >::type sampling_point(sampling_pointSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type month_vector(month_vectorSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type rainfall_relationship(rainfall_relationshipSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type rainfall_effect(rainfall_effectSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type likelihood_choice(likelihood_choiceSEXP);
    rcpp_result_gen = Rcpp::wrap(min_output_particle_filter(N, rainfall, obsData, fitted_parameters, static_parameters, density_function, sampling_point, month_vector, rainfall_relationship, rainfall_effect, likelihood_choice));
    return rcpp_result_gen;
END_RCPP
}
// full_output_particle_filter
Rcpp::List full_output_particle_filter(int N, std::vector <double> rainfall, std::vector <int> obsData, Rcpp::NumericVector fitted_parameters, Rcpp::NumericVector static_parameters, Rcpp::String density_function, double sampling_point, Rcpp::StringVector month_vector, Rcpp::String rainfall_relationship, Rcpp::String rainfall_effect, Rcpp::String likelihood_choice);
RcppExport SEXP _IndianVectorModelling_full_output_particle_filter(SEXP NSEXP, SEXP rainfallSEXP, SEXP obsDataSEXP, SEXP fitted_parametersSEXP, SEXP static_parametersSEXP, SEXP density_functionSEXP, SEXP sampling_pointSEXP, SEXP month_vectorSEXP, SEXP rainfall_relationshipSEXP, SEXP rainfall_effectSEXP, SEXP likelihood_choiceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< std::vector <double> >::type rainfall(rainfallSEXP);
    Rcpp::traits::input_parameter< std::vector <int> >::type obsData(obsDataSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type fitted_parameters(fitted_parametersSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type static_parameters(static_parametersSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type density_function(density_functionSEXP);
    Rcpp::traits::input_parameter< double >::type sampling_point(sampling_pointSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type month_vector(month_vectorSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type rainfall_relationship(rainfall_relationshipSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type rainfall_effect(rainfall_effectSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type likelihood_choice(likelihood_choiceSEXP);
    rcpp_result_gen = Rcpp::wrap(full_output_particle_filter(N, rainfall, obsData, fitted_parameters, static_parameters, density_function, sampling_point, month_vector, rainfall_relationship, rainfall_effect, likelihood_choice));
    return rcpp_result_gen;
END_RCPP
}
// Negative_Binomial
double Negative_Binomial(double k, double n, double r, double p);
RcppExport SEXP _IndianVectorModelling_Negative_Binomial(SEXP kSEXP, SEXP nSEXP, SEXP rSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type k(kSEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(Negative_Binomial(k, n, r, p));
    return rcpp_result_gen;
END_RCPP
}
// Poisson
double Poisson(double observed_data, double model_output, double population_fraction);
RcppExport SEXP _IndianVectorModelling_Poisson(SEXP observed_dataSEXP, SEXP model_outputSEXP, SEXP population_fractionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type observed_data(observed_dataSEXP);
    Rcpp::traits::input_parameter< double >::type model_output(model_outputSEXP);
    Rcpp::traits::input_parameter< double >::type population_fraction(population_fractionSEXP);
    rcpp_result_gen = Rcpp::wrap(Poisson(observed_data, model_output, population_fraction));
    return rcpp_result_gen;
END_RCPP
}
// Particle_Weight_Normalisation
std::vector <double> Particle_Weight_Normalisation(std::vector <double> particle_weights);
RcppExport SEXP _IndianVectorModelling_Particle_Weight_Normalisation(SEXP particle_weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector <double> >::type particle_weights(particle_weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(Particle_Weight_Normalisation(particle_weights));
    return rcpp_result_gen;
END_RCPP
}
// initial_state_sample
std::vector <int> initial_state_sample(Rcpp::NumericVector fitted_parameters, Rcpp::NumericVector static_parameters, double initial_K);
RcppExport SEXP _IndianVectorModelling_initial_state_sample(SEXP fitted_parametersSEXP, SEXP static_parametersSEXP, SEXP initial_KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type fitted_parameters(fitted_parametersSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type static_parameters(static_parametersSEXP);
    Rcpp::traits::input_parameter< double >::type initial_K(initial_KSEXP);
    rcpp_result_gen = Rcpp::wrap(initial_state_sample(fitted_parameters, static_parameters, initial_K));
    return rcpp_result_gen;
END_RCPP
}
// weighted_sampling_with_replacement
std::vector<int> weighted_sampling_with_replacement(int n, std::vector<double> prob, double p_tot, int K);
RcppExport SEXP _IndianVectorModelling_weighted_sampling_with_replacement(SEXP nSEXP, SEXP probSEXP, SEXP p_totSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type prob(probSEXP);
    Rcpp::traits::input_parameter< double >::type p_tot(p_totSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(weighted_sampling_with_replacement(n, prob, p_tot, K));
    return rcpp_result_gen;
END_RCPP
}
// runMCMC_joint_props
Rcpp::List runMCMC_joint_props(int start_sd_adaptation, int end_sd_adaptation, int number_of_iterations, std::vector <double> initial_sds, Rcpp::NumericVector model_parameters, Rcpp::NumericVector static_parameters, int N, std::vector <double> rainfall, std::vector <int> obsData, Rcpp::String density_function, Rcpp::String prior_choice, Rcpp::LogicalVector fitted_yn, double sampling_point, Rcpp::StringVector month_vector, Rcpp::String rainfall_relationship, Rcpp::String rainfall_effect, Rcpp::String likelihood_choice, int print_output);
RcppExport SEXP _IndianVectorModelling_runMCMC_joint_props(SEXP start_sd_adaptationSEXP, SEXP end_sd_adaptationSEXP, SEXP number_of_iterationsSEXP, SEXP initial_sdsSEXP, SEXP model_parametersSEXP, SEXP static_parametersSEXP, SEXP NSEXP, SEXP rainfallSEXP, SEXP obsDataSEXP, SEXP density_functionSEXP, SEXP prior_choiceSEXP, SEXP fitted_ynSEXP, SEXP sampling_pointSEXP, SEXP month_vectorSEXP, SEXP rainfall_relationshipSEXP, SEXP rainfall_effectSEXP, SEXP likelihood_choiceSEXP, SEXP print_outputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type start_sd_adaptation(start_sd_adaptationSEXP);
    Rcpp::traits::input_parameter< int >::type end_sd_adaptation(end_sd_adaptationSEXP);
    Rcpp::traits::input_parameter< int >::type number_of_iterations(number_of_iterationsSEXP);
    Rcpp::traits::input_parameter< std::vector <double> >::type initial_sds(initial_sdsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type model_parameters(model_parametersSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type static_parameters(static_parametersSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< std::vector <double> >::type rainfall(rainfallSEXP);
    Rcpp::traits::input_parameter< std::vector <int> >::type obsData(obsDataSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type density_function(density_functionSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type prior_choice(prior_choiceSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type fitted_yn(fitted_ynSEXP);
    Rcpp::traits::input_parameter< double >::type sampling_point(sampling_pointSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type month_vector(month_vectorSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type rainfall_relationship(rainfall_relationshipSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type rainfall_effect(rainfall_effectSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type likelihood_choice(likelihood_choiceSEXP);
    Rcpp::traits::input_parameter< int >::type print_output(print_outputSEXP);
    rcpp_result_gen = Rcpp::wrap(runMCMC_joint_props(start_sd_adaptation, end_sd_adaptation, number_of_iterations, initial_sds, model_parameters, static_parameters, N, rainfall, obsData, density_function, prior_choice, fitted_yn, sampling_point, month_vector, rainfall_relationship, rainfall_effect, likelihood_choice, print_output));
    return rcpp_result_gen;
END_RCPP
}
// testmvrnormArma
arma::mat testmvrnormArma(arma::mat mu, arma::mat sigma);
RcppExport SEXP _IndianVectorModelling_testmvrnormArma(SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(testmvrnormArma(mu, sigma));
    return rcpp_result_gen;
END_RCPP
}
// two_species_mosquito_population_model_with_comp_k_part
Rcpp::List two_species_mosquito_population_model_with_comp_k_part(int start_time, int end, Rcpp::NumericVector fitted_parameters_species_one, Rcpp::NumericVector fitted_parameters_species_two, Rcpp::NumericVector static_parameters, std::vector<double> rainfall, Rcpp::String mortality_density_function_species_one, Rcpp::String mortality_density_function_species_two, Rcpp::String rainfall_relationship_species_one, Rcpp::String rainfall_relationship_species_two, double competition_parameter_two_on_one, double competition_parameter_one_on_two, double species_one_static_K, double species_two_static_K, double threshold_dd_1, double threshold_dd_2, double max_K_1, double max_K_2);
RcppExport SEXP _IndianVectorModelling_two_species_mosquito_population_model_with_comp_k_part(SEXP start_timeSEXP, SEXP endSEXP, SEXP fitted_parameters_species_oneSEXP, SEXP fitted_parameters_species_twoSEXP, SEXP static_parametersSEXP, SEXP rainfallSEXP, SEXP mortality_density_function_species_oneSEXP, SEXP mortality_density_function_species_twoSEXP, SEXP rainfall_relationship_species_oneSEXP, SEXP rainfall_relationship_species_twoSEXP, SEXP competition_parameter_two_on_oneSEXP, SEXP competition_parameter_one_on_twoSEXP, SEXP species_one_static_KSEXP, SEXP species_two_static_KSEXP, SEXP threshold_dd_1SEXP, SEXP threshold_dd_2SEXP, SEXP max_K_1SEXP, SEXP max_K_2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type start_time(start_timeSEXP);
    Rcpp::traits::input_parameter< int >::type end(endSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type fitted_parameters_species_one(fitted_parameters_species_oneSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type fitted_parameters_species_two(fitted_parameters_species_twoSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type static_parameters(static_parametersSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type rainfall(rainfallSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type mortality_density_function_species_one(mortality_density_function_species_oneSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type mortality_density_function_species_two(mortality_density_function_species_twoSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type rainfall_relationship_species_one(rainfall_relationship_species_oneSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type rainfall_relationship_species_two(rainfall_relationship_species_twoSEXP);
    Rcpp::traits::input_parameter< double >::type competition_parameter_two_on_one(competition_parameter_two_on_oneSEXP);
    Rcpp::traits::input_parameter< double >::type competition_parameter_one_on_two(competition_parameter_one_on_twoSEXP);
    Rcpp::traits::input_parameter< double >::type species_one_static_K(species_one_static_KSEXP);
    Rcpp::traits::input_parameter< double >::type species_two_static_K(species_two_static_KSEXP);
    Rcpp::traits::input_parameter< double >::type threshold_dd_1(threshold_dd_1SEXP);
    Rcpp::traits::input_parameter< double >::type threshold_dd_2(threshold_dd_2SEXP);
    Rcpp::traits::input_parameter< double >::type max_K_1(max_K_1SEXP);
    Rcpp::traits::input_parameter< double >::type max_K_2(max_K_2SEXP);
    rcpp_result_gen = Rcpp::wrap(two_species_mosquito_population_model_with_comp_k_part(start_time, end, fitted_parameters_species_one, fitted_parameters_species_two, static_parameters, rainfall, mortality_density_function_species_one, mortality_density_function_species_two, rainfall_relationship_species_one, rainfall_relationship_species_two, competition_parameter_two_on_one, competition_parameter_one_on_two, species_one_static_K, species_two_static_K, threshold_dd_1, threshold_dd_2, max_K_1, max_K_2));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_IndianVectorModelling_mosquito_population_model_fluv", (DL_FUNC) &_IndianVectorModelling_mosquito_population_model_fluv, 8},
    {"_IndianVectorModelling_Hill_Function", (DL_FUNC) &_IndianVectorModelling_Hill_Function, 4},
    {"_IndianVectorModelling_test_initial_state_sample", (DL_FUNC) &_IndianVectorModelling_test_initial_state_sample, 4},
    {"_IndianVectorModelling_prior", (DL_FUNC) &_IndianVectorModelling_prior, 4},
    {"_IndianVectorModelling_likelihood_function", (DL_FUNC) &_IndianVectorModelling_likelihood_function, 11},
    {"_IndianVectorModelling_posterior_joint_proposals", (DL_FUNC) &_IndianVectorModelling_posterior_joint_proposals, 13},
    {"_IndianVectorModelling_joint_proposal_SD_adapter", (DL_FUNC) &_IndianVectorModelling_joint_proposal_SD_adapter, 7},
    {"_IndianVectorModelling_mosquito_population_model", (DL_FUNC) &_IndianVectorModelling_mosquito_population_model, 8},
    {"_IndianVectorModelling_mvrnormArma", (DL_FUNC) &_IndianVectorModelling_mvrnormArma, 2},
    {"_IndianVectorModelling_min_output_particle_filter", (DL_FUNC) &_IndianVectorModelling_min_output_particle_filter, 11},
    {"_IndianVectorModelling_full_output_particle_filter", (DL_FUNC) &_IndianVectorModelling_full_output_particle_filter, 11},
    {"_IndianVectorModelling_Negative_Binomial", (DL_FUNC) &_IndianVectorModelling_Negative_Binomial, 4},
    {"_IndianVectorModelling_Poisson", (DL_FUNC) &_IndianVectorModelling_Poisson, 3},
    {"_IndianVectorModelling_Particle_Weight_Normalisation", (DL_FUNC) &_IndianVectorModelling_Particle_Weight_Normalisation, 1},
    {"_IndianVectorModelling_initial_state_sample", (DL_FUNC) &_IndianVectorModelling_initial_state_sample, 3},
    {"_IndianVectorModelling_weighted_sampling_with_replacement", (DL_FUNC) &_IndianVectorModelling_weighted_sampling_with_replacement, 4},
    {"_IndianVectorModelling_runMCMC_joint_props", (DL_FUNC) &_IndianVectorModelling_runMCMC_joint_props, 18},
    {"_IndianVectorModelling_testmvrnormArma", (DL_FUNC) &_IndianVectorModelling_testmvrnormArma, 2},
    {"_IndianVectorModelling_two_species_mosquito_population_model_with_comp_k_part", (DL_FUNC) &_IndianVectorModelling_two_species_mosquito_population_model_with_comp_k_part, 18},
    {NULL, NULL, 0}
};

RcppExport void R_init_IndianVectorModelling(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
