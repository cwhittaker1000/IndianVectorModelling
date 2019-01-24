// [[Rcpp::depends(RcppArmadillo)]]
#include "MCMC_Functions.hpp"
#include "MVN_Sampler.hpp"

//' @export
// [[Rcpp::export]]
Rcpp::List runMCMC_joint_props(int N, // Number of particles,
                               int start_sd_adaptation, // Time to start covariance matrix adaptation
                               int end_sd_adaptation, // Time to stop covariance matrix adaptation
                               int number_of_iterations, // Number of MCMC iterations
                               std::vector <double> initial_sds, // Initial SDs to fill the covariance matrix diag with
                               Rcpp::NumericVector model_parameters, // Variable model parameter values - can be fitted
                               Rcpp::NumericVector static_parameters, // Static model parameter values - not fitted
                               Rcpp::LogicalVector fitted_yn, // Which parameters to be fitted
                               std::vector <double> rainfall, // Rainfall recorded
                               std::vector <int> obsData, // Observed data
                               Rcpp::String mortality_density_function, // Density function regulating mortality
                               Rcpp::String rainfall_relationship,
                               Rcpp::String rainfall_effect,
                               Rcpp::String decline_type,
                               double sampling_point,
                               Rcpp::StringVector offset_month_vector,
                               Rcpp::StringVector sampling_month_vector,
                               Rcpp::String likelihood_choice,
                               int print_output,
                               Rcpp::String calc_inside_mosquito_model) {

  // Removing overdispersion parameter if choice of likelihood is Poisson
  if (likelihood_choice == "poisson") {
    fitted_yn["overdisp"] = false;
  }

  // Storage for the MCMC chains, tracking acceptances/rejections etc
  Rcpp::NumericMatrix MCMC_chain_output(number_of_iterations + 1, Rcpp::sum(fitted_yn));
  Rcpp::NumericVector cholesky_decomp_failures(number_of_iterations + 1);
  Rcpp::List covariance_matrix_failures(number_of_iterations + 1);
  double acceptances = 0;
  double rejections = 0;
  Rcpp::NumericVector mu_tracker(Rcpp::sum(fitted_yn));

  double current_acceptance_ratio = 0;
  Rcpp::NumericVector total_acceptance_ratio_tracker(number_of_iterations + 1);
  Rcpp::NumericVector accepted_checker(number_of_iterations + 1);
  Rcpp::NumericVector accepted_status_checker(number_of_iterations + 1);

  double acceptances_adapt_only = 0;
  double cooldown;
  Rcpp::NumericVector post_adapt_acceptance_ratio_tracker(number_of_iterations + 1);

  // Storage for things that generate flexibility in choosing which parameters to fit
  Rcpp::StringVector parameter_names = fitted_yn.names();
  Rcpp::StringVector naming(Rcpp::sum(fitted_yn));
  Rcpp::NumericVector indexer(Rcpp::sum(fitted_yn));
  Rcpp::NumericMatrix latest_parameter_values_output(number_of_iterations + 1, Rcpp::sum(fitted_yn));

  // Assigning initial values to the first row of the MCMC output and naming that output's columns
  int p = 0;
  for (int d = 0; d < fitted_yn.size(); d++) {
    if (fitted_yn[d] == TRUE) {
      Rcpp::String to_be_added = parameter_names[d];
      naming[p] = to_be_added;
      indexer[p] = d;
      MCMC_chain_output(0, p) = model_parameters[to_be_added];
      p = p + 1;
    }
  }
  colnames(MCMC_chain_output) = naming;

  // Defining variables involved in the Adaptive MCMC part
  double accepted_variable;
  arma::mat current_covariance_matrix(Rcpp::sum(fitted_yn), Rcpp::sum(fitted_yn));
  for (int x = 0; x < Rcpp::sum(fitted_yn); x++) {

    int index = indexer[x];

    for (int y = 0; y < Rcpp::sum(fitted_yn); y++) {
      if (y == x) {
        current_covariance_matrix(x, y) = initial_sds[index];
      }
      else {
        current_covariance_matrix(x, y) = 0;
      }
    }
  }
  arma::mat initial_covariance_matrix_checker = current_covariance_matrix;

  //NOT SURE THIS IS WHAT I WANT TO BE ASSIGNING AS THE INITIAL VALUES FOR MU TBH.
  // this needs to be a ROW VECTOR at any rate
  // Update: Have moved the current_mu assignment to when adaptation starts i.e. will use chain
  // values at that point
  arma::mat current_mu(1, Rcpp::sum(fitted_yn));
  // for (int x = 0; x < Rcpp::sum(fitted_yn); x++) {
  //   int specific_index = indexer[x];
  //   current_mu(0, x) = model_parameters[specific_index];
  // }
  // arma::mat initial_mu_checker = current_mu;
  // Rcpp::Rcout << "The current mu is" << current_mu << std::endl;
  double current_scaling_factor = 1;

  // for (int x = 0; x < Rcpp::sum(fitted_yn); x++) {
  //   int specific_index = indexer[x];
  //   current_mu(0, x) = model_parameters[specific_index];
  // }
  // arma::mat initial_mu_checker = current_mu;

  // Calculating the posterior likelihood for the initial parameter values
  double current_posterior_likelihood = posterior(N, obsData,
                                                  model_parameters, static_parameters,
                                                  rainfall,
                                                  fitted_yn,
                                                  mortality_density_function, rainfall_relationship,
                                                  rainfall_effect, decline_type,
                                                  sampling_point, offset_month_vector,
                                                  sampling_month_vector, likelihood_choice,
                                                  calc_inside_mosquito_model); // whole bunch of inputs here, see if I can streamline this

  // Creating a vector to store and sequentially update model parameters as the MCMC iterates
  Rcpp::NumericVector model_parameters_for_MCMC(28);

  // Actual Model Parameters - See General_Mosquito_Model for Description of Which Parameter is at Which Position
  model_parameters_for_MCMC[0] = model_parameters[0]; model_parameters_for_MCMC[1] = model_parameters[1]; model_parameters_for_MCMC[2] = model_parameters[2];
  model_parameters_for_MCMC[3] = model_parameters[3]; model_parameters_for_MCMC[4] = model_parameters[4]; model_parameters_for_MCMC[5] = model_parameters[5];
  model_parameters_for_MCMC[6] = model_parameters[6]; model_parameters_for_MCMC[7] = model_parameters[7]; model_parameters_for_MCMC[8] = model_parameters[8];
  model_parameters_for_MCMC[9] = model_parameters[9]; model_parameters_for_MCMC[10] = model_parameters[10]; model_parameters_for_MCMC[11] = model_parameters[11];
  model_parameters_for_MCMC[12] = model_parameters[12]; model_parameters_for_MCMC[13] = model_parameters[13]; model_parameters_for_MCMC[14] = model_parameters[14];
  model_parameters_for_MCMC[15] = model_parameters[15]; model_parameters_for_MCMC[16] = model_parameters[16]; model_parameters_for_MCMC[17] = model_parameters[17];
  model_parameters_for_MCMC[18] = model_parameters[18]; model_parameters_for_MCMC[19] = model_parameters[19]; model_parameters_for_MCMC[20] = model_parameters[20];
  model_parameters_for_MCMC[21] = model_parameters[21]; model_parameters_for_MCMC[22] = model_parameters[22]; model_parameters_for_MCMC[23] = model_parameters[23];

  // E, L, P, M and Offset - these values are calculated within the Particle Filter function and the Values Added In There. 0 is just a placeholder!
  model_parameters_for_MCMC[24] = 0; model_parameters_for_MCMC[25] = 0; model_parameters_for_MCMC[26] = 0; model_parameters_for_MCMC[27] = 0;

  // Rcpp::NumericVector model_parameters_for_MCMC = Rcpp::NumericVector::create(Rcpp::Named("dE") = model_parameters[0], Rcpp::Named("dL") = model_parameters[1],
  //                                                                             Rcpp::Named("dP") = model_parameters[2], Rcpp::Named("muE0") = model_parameters[3],
  //                                                                             Rcpp::Named("muL0") = model_parameters[4], Rcpp::Named("muP") = model_parameters[5],
  //                                                                             Rcpp::Named("muM") = model_parameters[6], Rcpp::Named("lambda") = model_parameters[7],
  //                                                                             Rcpp::Named("tau") = model_parameters[8], Rcpp::Named("beta") = model_parameters[9],
  //                                                                             Rcpp::Named("overdisp") = model_parameters[10], Rcpp::Named("pop_frac") = model_parameters[11],
  //                                                                             Rcpp::Named("scaling_factor") = model_parameters[12], Rcpp::Named("z") = model_parameters[13],
  //                                                                             Rcpp::Named("K_static") = model_parameters[14], Rcpp::Named("K_max") = model_parameters[15],
  //                                                                             Rcpp::Named("hill_1") = model_parameters[16], Rcpp::Named("hill_2") = model_parameters[17],
  //                                                                             Rcpp::Named("E") = 0.0, Rcpp::Named("L") = 0.0,
  //                                                                             Rcpp::Named("P") = 0.0, Rcpp::Named("M") = 0.0);
  int cov_fail_counter = 0;
  // Iterating over each parameter for a specified number of iterations
  for (int i = 0; i < number_of_iterations; i++) {

    arma::vec eigenvalues = arma::eig_sym(current_covariance_matrix);
    if (eigenvalues.min() <= 0.0) {
      cholesky_decomp_failures[i] = 1;
      covariance_matrix_failures[cov_fail_counter] = current_covariance_matrix;
      cov_fail_counter = cov_fail_counter + 1;
    }

    // Generating the proposed parameter values
    // CURRENT_PARAMETER_VALUES needs to be a ROW VECTOR!!!!
    arma::mat current_parameter_values(1, Rcpp::sum(fitted_yn));
    for (int u = 0; u < Rcpp::sum(fitted_yn); u++) {
      current_parameter_values(0, u) = MCMC_chain_output(i, u);
    }

    //Rcpp::Rcout << "The iteration number is " << i << std::endl;
    //Rcpp::Rcout << "The current parameter values are" << current_parameter_values << std::endl;
    //Rcpp::Rcout << "The current covariance matrix is " << current_covariance_matrix << std::endl;
    //Rcpp::Rcout << "The current lowest eigenvalue is " << eigenvalues.min() << std::endl;
    //Rcpp::Rcout << "The number of failures so far is" << cov_fail_counter << std::endl;
    arma::mat proposed_parameter_values = mvrnormArma(current_parameter_values, current_covariance_matrix);
    //Rcpp::Rcout << "I've actually made it past the MVN call- weird!" << std::endl;

    // Filling input parameter vector with the output from the proposal function. Fixed parameters remain unchanged.
    for (int k = 0; k < indexer.size(); k++) {
      int index_for_adding = indexer[k];
      model_parameters_for_MCMC[index_for_adding] = proposed_parameter_values[k];
    }

    // This might be superfluous
    Rcpp::StringVector names(28);
    names[0] = "dE"; names[1] = "dL"; names[2] = "dP"; names[3] = "muE0"; names[4] = "muL0"; names[5] = "muP"; names[6] = "muM";
    names[7] = "lambda"; names[8] = "beta"; names[9] = "overdisp"; names[10] = "pop_frac"; names[11] = "z";
    names[12] = "tau_rainfall"; names[13] = "scaling_factor_rainfall"; names[14] = "K_Max_Rain"; names[15] = "Hill_Rainfall_1"; names[16] = "Hill_Rainfall_2";
    names[17] = "tau_static"; names[18] = "scaling_factor_static"; names[19] = "K_Max_Static"; names[20] = "Washout_Threshold"; names[21] = "washout_exp_scaling_factor";
    names[22] = "washout_hill_one"; names[23] = "washout_hill_two";
    names[24] = "E"; names[25] = "L"; names[26] = "P"; names[27] = "M";
    model_parameters_for_MCMC.names() = names;

    // Manually handling instances when proposed overdispersion value is less than or equal to 0
    // DOUBLE CHECKED AND FAIRLY CERTAIN THE BELOW ISSUE IS SOLVED :)
    // CURRENTLY THIS CHANGES THE MODEL INPUT BUT DOESN'T CHANGE WHAT GETS ADDED TO THE MCMC OUTPUT
    // IF THE PROPOSED VALUE IS ACCEPTED- THAT IS STILL WHAT WAS INITIALLY PROPOSED AND CONTAINED
    // WITHIN "PROPOSED_PARAMETERS". THINK THE SOLUTION IS TO USE INDEXER AND ADD TO MCMC_OUTPUT
    // FROM MODEL_PARAMETERS_FOR_MCMC DIRECTLY.
    // think I've solved this by changing line 206 such that the parameter value is added
    // from the model_parameters_for_MCMC using indexer rather than directly from the proposed_parameters
    // vector which actually isn't ever altered from the set of parameters proposed.
    // Consider adding a for loop in here that basically loops over the entire vector, and if it sees any instances
    // of 0, changes them to a very very small number.
    if (model_parameters_for_MCMC["dP"] <= 0) {
      model_parameters_for_MCMC["dP"] = 0.01;
    }
    if (model_parameters_for_MCMC["muE0"] <= 0) {
      model_parameters_for_MCMC["muE0"] = 0.001;
    }
    if (model_parameters_for_MCMC["muL0"] <= 0) {
      model_parameters_for_MCMC["muL0"] = 0.001;
    }
    if (model_parameters_for_MCMC["muP"] <= 0) {
      model_parameters_for_MCMC["muP"] = 0.001;
    }
    if (model_parameters_for_MCMC["muM"] <= 0) {
      model_parameters_for_MCMC["muM"] = 0.001;
    }
    if (model_parameters_for_MCMC["overdisp"] <= 0) {
      model_parameters_for_MCMC["overdisp"] = 0.001;
    }
    if (model_parameters_for_MCMC["pop_frac"] <= 0) {
      model_parameters_for_MCMC["pop_frac"] = 0.001;
    }
    if (model_parameters_for_MCMC["lambda"] <= 0) {
      model_parameters_for_MCMC["lambda"] = 0.001;
    }
    if (model_parameters_for_MCMC["z"] < 1) {
      model_parameters_for_MCMC["z"] = 1;
    }
    if (model_parameters_for_MCMC["scaling_factor_rainfall"] <= 0) {
      model_parameters_for_MCMC["scaling_factor_rainfall"] = 0.1;
    }
    if (model_parameters_for_MCMC["scaling_factor_static"] <= 0) {
      model_parameters_for_MCMC["scaling_factor_static"] = 0.1;
    }
    if (model_parameters_for_MCMC["tau_rainfall"] < 1) {
      model_parameters_for_MCMC["tau_rainfall"] = 1;
    }
    if (model_parameters_for_MCMC["tau_static"] < 1) {
      model_parameters_for_MCMC["tau_static"] = 1;
    }
    if (model_parameters_for_MCMC["K_Max_Static"] < 0) {
      model_parameters_for_MCMC["K_Max_Static"] = 0;
    }

    // Assessing the likelihood of the proposed parameter value compared to the current parameter value
    double proposed_posterior_likelihood = posterior(N, obsData,
                                                     model_parameters_for_MCMC, static_parameters,
                                                     rainfall, fitted_yn,
                                                     mortality_density_function, rainfall_relationship, rainfall_effect, decline_type,
                                                     sampling_point, offset_month_vector, sampling_month_vector,
                                                     likelihood_choice, calc_inside_mosquito_model);

    double likelihood_ratio = exp(proposed_posterior_likelihood - current_posterior_likelihood);

    // Deciding whether to accept or reject, updating chain and current_posterior_likelihood as appropriate
    if(R::runif(0, 1) < likelihood_ratio) {

      // Taking arma mat output proposed parameter values and sequentially adding the results to next row of
      // the MCMC chain (if these parameter values were accepted)
      for (int u = 0; u < proposed_parameter_values.size(); u++) {
        int index_for_adding = indexer[u];
        MCMC_chain_output(i + 1, u) = model_parameters_for_MCMC[index_for_adding];
      }
      current_posterior_likelihood = proposed_posterior_likelihood;
      acceptances = acceptances + 1;
      accepted_checker[i] = 2;
      accepted_variable = 1;
    }
    else {

      // Taking arma mat output current parameter values and sequentially adding these same parameter values
      // to next row of the MCMC chain (as the proposed parameters were rejected)
        for (int u = 0; u < current_parameter_values.size(); u++) {
          MCMC_chain_output(i + 1, u) = current_parameter_values[u];
        }
        current_posterior_likelihood = current_posterior_likelihood; // remove this as it's superfluous- just to tell myself what's going on
        rejections = rejections + 1;
        accepted_checker[i] = 1;
        accepted_variable = 0;
      }

      // Tracking the acceptance ratio
      current_acceptance_ratio = (acceptances / (i + 1));
      total_acceptance_ratio_tracker[i] = current_acceptance_ratio;

      accepted_status_checker[i] = accepted_variable;

      // Adapting the standard deviation of parameter proposals to ensure a decent acceptance ratio
      if (i > start_sd_adaptation & i < end_sd_adaptation) { // breaks if I change this to i >= start_sd_adaptation and have no idea why!

        if (i == start_sd_adaptation + 1) {
          for (int x = 0; x < Rcpp::sum(fitted_yn); x++) {
            //int specific_index = indexer[x];
            current_mu(0, x) = MCMC_chain_output(i + 1, x);
          }
          arma::mat initial_mu_checker = current_mu;
          if (print_output == 1) {
            Rcpp::Rcout << "Success- current mu at time of adaptation is " << current_mu << std::endl;
          }
        }
        // Considered using current_parameter_values here but that isn't updated (if the proposed values are accepted)
        // until the start of the next iteration. The MCMC Output matrix IS however and so I use that to fill a temporary
        // vector containing what either are currently (if proposed were rejected) the parameter values OR what are
        // going to be assigned to current_parameter_values in the beginning of the next iteration (if proposed were accepted).
        arma::mat latest_parameter_values_but_col = MCMC_chain_output(i + 1, Rcpp::_); //double check this should be i + 1 and not i
        arma::mat latest_parameter_values = latest_parameter_values_but_col.t();

        // Tracking acceptance ratio only when adaptation starts
        if (accepted_variable == 1) {
          acceptances_adapt_only = acceptances_adapt_only + 1;
        }

        post_adapt_acceptance_ratio_tracker[i] = acceptances_adapt_only/(i - start_sd_adaptation);

        for (int h = 0; h < Rcpp::sum(fitted_yn); h++) {
          latest_parameter_values_output(i, h) = latest_parameter_values[h];
        }

        Rcpp::List adapter_output = proposal_SD_adapter(accepted_variable, i, start_sd_adaptation,
                                                        current_scaling_factor, current_mu, latest_parameter_values, // technically parameter values for t + 1 as the acceptance/rejection step precedes calling this function
                                                        current_covariance_matrix);

        // Currently having to take the output and assign it to Rcpp types before adding those elements to arma::mat/vec types.
        // Ask Rich if there's a nicer way of doing this.
        Rcpp::NumericVector new_mu = adapter_output["New_Mu"];
        Rcpp::NumericMatrix new_covariance_matrix = adapter_output["New_Covariance_Matrix"];
        cooldown = adapter_output["Cooldown"];

        for (int u = 0; u < Rcpp::sum(fitted_yn); u++) {
          current_mu[u] = new_mu[u];
        }
        for (int u = 0; u < current_parameter_values.size(); u++) {
          for (int v = 0; v < current_parameter_values.size(); v++) {
            current_covariance_matrix(u, v) = new_covariance_matrix(u, v);
          }
        }
        current_scaling_factor = adapter_output["New_Scaling_Factor"];

      }

      if (print_output == 1) {
        if (i == 1 | i == 2 | i == 3 | i == 4 | i == 5 | i == 6 | i == 7 | i == 8 | i == 9 | i == 10 | i % 100 == 0) {
          Rcpp::Rcout << current_covariance_matrix << "This is the covariance matrix" << std::endl;
          Rcpp::Rcout << "The overall acceptance ratio is " << current_acceptance_ratio << std::endl;
          Rcpp::Rcout << "The iteration number is " << i << std::endl;
          if (i > start_sd_adaptation & i < end_sd_adaptation) {
            Rcpp::Rcout << "The current scaling factor is " << current_scaling_factor << std::endl;
            Rcpp::Rcout << "The current mu is " << current_mu << std::endl;
            Rcpp::Rcout << "The cooldown factor is " << cooldown << std::endl;
            Rcpp::Rcout << "The acceptance ratio since starting to adapt is " << post_adapt_acceptance_ratio_tracker[i] << "\n" << std::endl;
            Rcpp::Rcout << "Number of cholesky failuress is " << cov_fail_counter << std::endl;
            Rcpp::Rcout << "The iteration number is " << i << std::endl;
          }
          else {
            Rcpp::Rcout << "\n " << std::endl;
          }
        }
      }
  }

  return(Rcpp::List::create(Rcpp::Named("MCMC_Output") = MCMC_chain_output,
                            Rcpp::Named("Acceptance_Ratio_Total") = current_acceptance_ratio,
                            Rcpp::Named("Acceptance_Ratio_Post_Adapt") = post_adapt_acceptance_ratio_tracker,
                            Rcpp::Named("Acceptance_Checker") = accepted_checker,
                            Rcpp::Named("Acceptances_Post_Adapt") = acceptances_adapt_only,
                            Rcpp::Named("Acceptances_Total")= acceptances,
                            Rcpp::Named("Accepted_Variable_Checker") = accepted_status_checker,
                            Rcpp::Named("initial_covariance_matrix") = initial_covariance_matrix_checker,
                            Rcpp::Named("final_covariance_matrix") = current_covariance_matrix,
                            Rcpp::Named("Mu_Output") = latest_parameter_values_output,
                            Rcpp::Named("Iteration of Cholesky Failures") = cholesky_decomp_failures,
                            Rcpp::Named("Covariance Failures") = covariance_matrix_failures,
                            Rcpp::Named("Number of Chol Failures") = cov_fail_counter));
}
