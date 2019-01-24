/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
///                                                                                                   ///
///  Running the Adaptive Metropolis Hastings MCMC With a Particle Filter                             ///
///                                                                                                   ///
///  Charlie Whittaker                                                                                ///
///  Imperial College London                                                                          ///
///  charles.whittaker16@imperial.ac.uk                                                               ///
///                                                                                                   ///
///  The code below specifies the function run_particle_MCMC. This function utilises functions        ///
///  contained in the files:                                                                          ///
///                                                                                                   ///
///       1) Mosquito_Model.cpp - which encodes a stochastic model of mosquito population dynamics.   ///
///       2) Particle_Filter_Functions.cpp - which encodes a Particle Filter that produces an         ///
///          estimate of the mosquito model likelihood given some data and a parameter set.           ///
///       3) MCMC_functions.cpp - which encode the various functions required to run a                ///
///          Metrpolois-Hastings based MCMC, as well as a function allowing sampled from a            ///
///          multivariate normal distribution (required for parameter proposals). It also encodes     ///
///          an algorithm that allows sequential and continuous updating of the covariance matrix     ///
///          that parameterises the multivariate normal distribution.                                 ///
///                                                                                                   ///                                      ///
///                                                                                                   ///
///  The code that comprises run_particle_MCMC is structured as follows:                              ///
///                                                                                                   ///
///       1) PARAMETER SETTING & INITIALISATION                                                       ///
///       2) IMPLEMENTING THE ACTUAL METROPOLIS-HASTINGS ALGORITHM- PROPOSING NEW PARAMETERS,         ///
///          EVALUATING THE POSTERIOR DENSITY OF THESE PARAMETERS, AND THEN DECIDING TO ACCEPT OR     ///
///          REJECT THESE PARAMETERS ACCORDING TO THE LIKELIHOOD RATIO I.E. IN A WAY THAT IS          ///
///          PROPORTIONAL TO THEIR POSTERIOR DENSITY COMPARED TO THE POSTERIOR DENSITY OF THE         ///
///          PREVIOUS PARAMETER SET.                                                                  ///
///       3) ADAPTING THE COVARIANCE MATRIX CURRENTLY BEING USED TO PROPOSE PARAMETERS BASED          ///
///          DEPENDING ON WHETHER OR NOT THE PARAMETERS PROPOSED PREVIOUSLY WERE ACCEPTED OR NOT.     ///
///       4) (optional) PRINTING OUT VARIOUS USEFUL MCMC DIAGNOSTICS REALTIME WHILST THE MODEL IS     ///
///          RUNNING                                                                                  ///
///                                                                                                   ///
///  Below the code for the Mosquito Population Model, there is code for two other functions,         ///
///  including a Hill Function (returning some inputs put through a Hill Function) and also the       ///
///  function used to produce Initial Values (given a set of parameter values) to input into the      ///
///  model.                                                                                           ///
///                                                                                                   ///
///                                                                                                   ///
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

// Specifying the Includes and Depends Required
#include "MCMC_Functions.hpp"
// [[Rcpp::depends(RcppArmadillo)]]

//' @export
// [[Rcpp::export]]
Rcpp::List run_particle_MCMC(int N, // Number of particles,
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
                             Rcpp::String rainfall_relationship, // How present and past rainfall is weighted (mean, linear, exponential)
                             Rcpp::String rainfall_effect, // How weighted rainfall is incorporated- raw or Hill
                             Rcpp::String decline_type,  // How dryout occurs- exponential or Hill based decline
                             double sampling_point,  // When in the month mosquito sampling happens
                             Rcpp::StringVector offset_month_vector, // Months preceding when we have data (offset)
                             Rcpp::StringVector sampling_month_vector, // Months in which sampling took place
                             Rcpp::String likelihood_choice, // Negative binomial or Poisson likelihood
                             int print_output,  // WHehter or not to print out key diagnostics whilst running
                             Rcpp::String calc_inside_mosquito_model) {  // For exponential weighting relationship, where to calculate the exponential weighting factors

  /////////////////////////////////////////////////////////////////////////////////
  //                                                                             //
  //  1.  Initialising Variables and Preparing Various Aspects to Run the MCMC   //
  //                                                                             //
  /////////////////////////////////////////////////////////////////////////////////

  // This section of code initialises all the variables, storage vectors, storage matrices required to run the MCMC, and
  // additionally does any general housekeeping to get everything in order before the MCMC is ran. This includes adapting
  // the fitted_yn vector to set the overdispersion parameter to be false if the Poisson likelihood is chosen, creates a
  // logical vector specifying which model parameters are to be fitted,

  // Removing Overdispersion Parameter if Choice of Likelihood is Poisson
  if (likelihood_choice == "poisson") {
    fitted_yn["overdisp"] = false;
  }

  // Initialising the Objects, Vectors and Matrices Required to Store Outputs of Interest
  Rcpp::NumericMatrix MCMC_chain_output(number_of_iterations + 1, Rcpp::sum(fitted_yn));  // stores MCMC output
  Rcpp::NumericVector cholesky_decomp_failures(number_of_iterations + 1); // tracks Cholesky Decomposition failures (matrices that can't be inverted)
  Rcpp::List covariance_matrix_failures(number_of_iterations + 1); // stores covariances matrices for which the Cholesky Decomposition fails
  Rcpp::NumericVector post_adapt_acceptance_ratio_tracker(number_of_iterations + 1); // Stores acceptance ratio at each point following starting adaptation
  Rcpp::NumericVector total_acceptance_ratio_tracker(number_of_iterations + 1); // Tracks the overall acceptance ratio at each point
  Rcpp::NumericVector accepted_status_checker(number_of_iterations + 1); // Tracks whether each proposal was rejected or not
  double current_acceptance_ratio = 0; // Tracks current acceptance ratio
  double acceptances = 0; // Tracks number of acceptances
  double rejections = 0; // Tracks number of rejections
  int cov_fail_counter = 0; // Counter to track number of times covariance matrix cannot be Cholesky decomposed
  double acceptances_adapt_only = 0; // Tracks number of acceptances since beginning adaptation
  Rcpp::NumericMatrix latest_parameter_values_output(number_of_iterations + 1, Rcpp::sum(fitted_yn)); // compare to MCMC chans to make sure right stuff is being put in

  // Generates an indexer from which to access fitted parameters from full parameters vector, fills the names vector with
  // parameters that are to be fitted, and fills the first row of the MCMC output
  int p = 0;
  Rcpp::StringVector parameter_names = fitted_yn.names(); // vector of the names of all parameters in the logical vector fitted_yn
  Rcpp::StringVector naming(Rcpp::sum(fitted_yn)); // vector to contain only the names of those parameters that are being fitted
  Rcpp::NumericVector indexer(Rcpp::sum(fitted_yn)); // inf
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
  double current_scaling_factor = 1;
  arma::mat current_mu(1, Rcpp::sum(fitted_yn)); // initialised later on with the MCMC chain values at the iteration number when adaptation starts
  arma::mat current_covariance_matrix(Rcpp::sum(fitted_yn), Rcpp::sum(fitted_yn)); // fills the covariance matrix diagonalswith the sds vector
  for (int x = 0; x < Rcpp::sum(fitted_yn); x++) {                                 // for those parameters that are being fitted
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
  double accepted_variable; // tracks whether most recent parameters were accepted or rejected
  double cooldown; // The cooldown factor in the covariance matrix adaptation algorithm
  arma::mat initial_covariance_matrix_checker = current_covariance_matrix; // checks the covariance matrix is being initialised properly

  // Calculating the Posterior Density for the Initial Parameter Values
  double current_posterior_likelihood = posterior(N, obsData,
                                                  model_parameters, static_parameters,
                                                  rainfall,
                                                  fitted_yn,
                                                  mortality_density_function, rainfall_relationship,
                                                  rainfall_effect, decline_type,
                                                  sampling_point, offset_month_vector,
                                                  sampling_month_vector, likelihood_choice,
                                                  calc_inside_mosquito_model);

  // Creating a vector to store and sequentially update model parameters as the MCMC iterates
  Rcpp::NumericVector model_parameters_for_MCMC(model_parameters.size()); // NOTE: this needs to be the length of the model_parameters vector. Change 29 to model_parameters.size()

  // Actual Model Parameters - See General_Mosquito_Model for Description of Which Parameter is at Which Position
  // THIS MIGHT BE SUPERFLUOUS
  model_parameters_for_MCMC[0] = model_parameters[0]; model_parameters_for_MCMC[1] = model_parameters[1]; model_parameters_for_MCMC[2] = model_parameters[2];
  model_parameters_for_MCMC[3] = model_parameters[3]; model_parameters_for_MCMC[4] = model_parameters[4]; model_parameters_for_MCMC[5] = model_parameters[5];
  model_parameters_for_MCMC[6] = model_parameters[6]; model_parameters_for_MCMC[7] = model_parameters[7]; model_parameters_for_MCMC[8] = model_parameters[8];
  model_parameters_for_MCMC[9] = model_parameters[9]; model_parameters_for_MCMC[10] = model_parameters[10]; model_parameters_for_MCMC[11] = model_parameters[11];
  model_parameters_for_MCMC[12] = model_parameters[12]; model_parameters_for_MCMC[13] = model_parameters[13]; model_parameters_for_MCMC[14] = model_parameters[14];
  model_parameters_for_MCMC[15] = model_parameters[15]; model_parameters_for_MCMC[16] = model_parameters[16]; model_parameters_for_MCMC[17] = model_parameters[17];
  model_parameters_for_MCMC[18] = model_parameters[18]; model_parameters_for_MCMC[19] = model_parameters[19]; model_parameters_for_MCMC[20] = model_parameters[20];
  model_parameters_for_MCMC[21] = model_parameters[21]; model_parameters_for_MCMC[22] = model_parameters[22]; model_parameters_for_MCMC[23] = model_parameters[23];
  model_parameters_for_MCMC[24] = model_parameters[24];

  // E, L, P, M and Offset - these values are calculated within the Particle Filter function and the Values Added In There. 0 is just a placeholder!
  model_parameters_for_MCMC[25] = 0; model_parameters_for_MCMC[26] = 0; model_parameters_for_MCMC[27] = 0; model_parameters_for_MCMC[28] = 0;

  ///////////////////////////////////////////////////////////////////////
  //                                                                   //
  //  2.  Implementing Metropolis-Hastings Particle Filter Algorithm   //
  //                                                                   //
  ///////////////////////////////////////////////////////////////////////

  // This section of code implements the Particle MCMC algorithm. Specifically, it proposes a new set of parameters based on the
  // values currently in the MCMC chain. It then evaluates the Posterior Density of these new values and then either accepts or
  // rejects these values based upon the ratio of the Posterior Density for the current and proposed parameter values. The
  // next row of the MCMC chain is then filled with either the proposed values, or the current values, depending on whether the
  // proposed values were accepted or rejected.

  // Running the MCMC number_of_iterations times
  for (int i = 0; i < number_of_iterations; i++) {

    // Checking whether any of the covariance matrix eigenvalues are < 0 and hence the Cholesky Decomposition will fail
    arma::vec eigenvalues = arma::eig_sym(current_covariance_matrix);
    if (eigenvalues.min() <= 0.0) {
      cholesky_decomp_failures[i] = 1;
      covariance_matrix_failures[cov_fail_counter] = current_covariance_matrix;
      cov_fail_counter = cov_fail_counter + 1;
    }

    // Generating the Proposed Parameter Values
    arma::mat current_parameter_values(1, Rcpp::sum(fitted_yn));
    for (int u = 0; u < Rcpp::sum(fitted_yn); u++) {
      current_parameter_values(0, u) = MCMC_chain_output(i, u); // Note: current_parameter_values needs to be a row vector
    }
    //Rcpp::Rcout << "The iteration number is " << i << std::endl;
    //Rcpp::Rcout << "The current parameter values are" << current_parameter_values << std::endl;
    //Rcpp::Rcout << "The current covariance matrix is " << current_covariance_matrix << std::endl;
    //Rcpp::Rcout << "The current lowest eigenvalue is " << eigenvalues.min() << std::endl;
    //Rcpp::Rcout << "The number of failures so far is" << cov_fail_counter << std::endl;
    arma::mat proposed_parameter_values = mvrnormArma(current_parameter_values, current_covariance_matrix);
    //Rcpp::Rcout << "I've actually made it past the MVN call- weird!" << std::endl;

    // Filling the Parameter Vector with the Output from the Proposal Function. Fixed Parameters Remain Unchanged.
    for (int k = 0; k < indexer.size(); k++) {
      int index_for_adding = indexer[k];
      model_parameters_for_MCMC[index_for_adding] = proposed_parameter_values[k];
    }

    // Naming the Parameter Vector contents so that it plays nicely with all the other model/particle filter functions (Might be superfluous)
    Rcpp::StringVector names(29); // has to be number of parameters - i.e. whatever index M is + 1 (because C++ indexing starts from 0)
    names[0] = "dE"; names[1] = "dL"; names[2] = "dP"; names[3] = "muE0"; names[4] = "muL0"; names[5] = "muP"; names[6] = "muM";
    names[7] = "lambda"; names[8] = "beta"; names[9] = "overdisp"; names[10] = "pop_frac"; names[11] = "z";
    names[12] = "tau_rainfall"; names[13] = "scaling_factor_rainfall"; names[14] = "K_Max_Rain"; names[15] = "Hill_Rainfall_1"; names[16] = "Hill_Rainfall_2";
    names[17] = "tau_static"; names[18] = "scaling_factor_static"; names[19] = "K_Max_Static"; names[20] = "Washout_Occurrence", names[21] = "Washout_Threshold";
    names[22] = "washout_exp_scaling_factor"; names[23] = "washout_hill_one"; names[24] = "washout_hill_two";
    names[25] = "E"; names[26] = "L"; names[27] = "P"; names[28] = "M";
    model_parameters_for_MCMC.names() = names;

    // Manually handling instances when proposed values are outside of range that can be handled by the model
        // Consider adding some sort of loop in here that basically loops over the vector, and if it sees any instances
        // of 0 for these parameters, changes them to a very very small number.
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

    // Evaluating the Posterior Density of the Proposed Parameters, the Likelihood Ratio and then Whether to Accept the Proposed Values
    double proposed_posterior_likelihood = posterior(N, obsData,
                                                     model_parameters_for_MCMC, static_parameters,
                                                     rainfall, fitted_yn,
                                                     mortality_density_function, rainfall_relationship, rainfall_effect, decline_type,
                                                     sampling_point, offset_month_vector, sampling_month_vector,
                                                     likelihood_choice, calc_inside_mosquito_model);
    double likelihood_ratio = exp(proposed_posterior_likelihood - current_posterior_likelihood);

    // If accepted, add the proposed parameters to the next row of the MCMC chain
    if(R::runif(0, 1) < likelihood_ratio) {
      for (int u = 0; u < proposed_parameter_values.size(); u++) {
        int index_for_adding = indexer[u];
        MCMC_chain_output(i + 1, u) = model_parameters_for_MCMC[index_for_adding];
      }
      current_posterior_likelihood = proposed_posterior_likelihood;
      acceptances = acceptances + 1;
      accepted_variable = 1;
    }
    // If rejected, add the current parameters to the next row of the MCMC chain
    else {
      for (int u = 0; u < current_parameter_values.size(); u++) {
        MCMC_chain_output(i + 1, u) = current_parameter_values[u];
      }
      current_posterior_likelihood = current_posterior_likelihood; // remove this as it's superfluous- just to tell myself what's going on
      rejections = rejections + 1;
      accepted_variable = 0;
    }

    // Tracking the Acceptance Ratio and Other Related Diagnostics
    current_acceptance_ratio = (acceptances / (i + 1));
    total_acceptance_ratio_tracker[i] = current_acceptance_ratio;
    accepted_status_checker[i] = accepted_variable;


    //////////////////////////////////////////////////////////////////////
    //                                                                  //
    //  3.  Adapting the Covariance Matrix Used to Propose Parameters   //
    //                                                                  //
    //////////////////////////////////////////////////////////////////////

    // This section of code uses the function proposal_SD_adapter to adapt the shape and size of the covariance matrix that
    // parameterises the multivariate normal distribution from which new parameter values (or more accurately, vectors) are
    // proposed. This adaptation is carried out in order to improve the efficiency of the MCMC sampling, and to achieve an
    // acceptance ratio near to 0.25 (which previous work has shown to be optimal for efficient parameter space exploration).
    //The adaptation algorithm comes from the Supplementary Information of Johnstone et al., 2016 JMCC
    // (further details available here: https://www.sciencedirect.com/science/article/pii/S0022282815301231). In brief, the
    // covariance matrix is adapted depending on whether the most recently proposed set of parameters were accepted (in which
    // case accepted_variable = 1) or rejected (in which case accepted_variable = 0). The extent of the changes that get made
    // to the covariance matrix over time decrease- that is, the covariance matrix is changed a lot at the beginning of the
    // MCMC, and progressively less so as iterations increase. This is to satisfy a number of Markov-Chain properties and ensure
    // the chain converges to the actual Posterior Distribution.

    // Adapting the Covariance Matrix of Parameter Proposals
    if (i > start_sd_adaptation & i < end_sd_adaptation) { // Note: The MCMC breaks if I change this to i >= start_sd_adaptation and have no idea why!

      // Note: initially had some issues with the algorithm. In the paper, before running the MCMC, they find a MLE (or similar),
      // and then initiate their MCMC from there, using the MLE parameter values as their current_mu. I'm not doing this, and an
      // issue I initially had if starting adaptation some 1000s of iterations in, but using a current_mu specified by my initial
      // parameter values. This issue typically manifested in the form of the Cholesky Decompisition failing. Not convinced that
      // was the sole reason for the failure, but to get round this, I've encoded the below such that when adaptation first begins
      // (i.e. i == start_adaptation + 1) current_mu is specified THEN using the current parameter values in the MCMC chain.
      if (i == start_sd_adaptation + 1) {
        for (int x = 0; x < Rcpp::sum(fitted_yn); x++) {
          current_mu(0, x) = MCMC_chain_output(i + 1, x);
        }
        // arma::mat initial_mu_checker = current_mu; Was used in checking the adaptation algorithm
        if (print_output == 1) {
          Rcpp::Rcout << "Success- current mu at time of adaptation is " << current_mu << std::endl; // prints relevant output if required
        }
      }
      // Considered using current_parameter_values here but that isn't updated (if the proposed values are accepted)
      // until the start of the next iteration. The MCMC Output matrix IS however and so I use that to fill a temporary
      // vector containing what either are currently (if proposed were rejected) the parameter values OR what are
      // going to be assigned to current_parameter_values in the beginning of the next iteration (if proposed were accepted).
      arma::mat latest_parameter_values_but_col = MCMC_chain_output(i + 1, Rcpp::_); //double check this should be i + 1 and not i
      arma::mat latest_parameter_values = latest_parameter_values_but_col.t();

      // Tracking the Acceptance Ratio and Related Diagnostics
      if (accepted_variable == 1) {
        acceptances_adapt_only = acceptances_adapt_only + 1; // total acceptances since adaptation began
      }
      post_adapt_acceptance_ratio_tracker[i] = acceptances_adapt_only/(i - start_sd_adaptation); // acceptance ratio since adaptation began
      for (int h = 0; h < Rcpp::sum(fitted_yn); h++) {
        latest_parameter_values_output(i, h) = latest_parameter_values[h]; // tracking latest parameter values, ensure they match what they should be
      }

      // Adaptation and Alteration of the Covariance Matrix and Related Variables
      Rcpp::List adapter_output = proposal_SD_adapter(accepted_variable, i, start_sd_adaptation,
                                                      current_scaling_factor, current_mu, latest_parameter_values, // technically parameter values for t + 1 as the acceptance/rejection step precedes calling this function
                                                      current_covariance_matrix);

      Rcpp::NumericVector new_mu = adapter_output["New_Mu"]; // adapter_outputs's new_mu, which is then assigned to current_mu
      for (int u = 0; u < Rcpp::sum(fitted_yn); u++) {
        current_mu[u] = new_mu[u];
      }
      Rcpp::NumericMatrix new_covariance_matrix = adapter_output["New_Covariance_Matrix"]; // adapter_output's new covariance matrix, which is assigned to be the current covariance matrix
      for (int u = 0; u < current_parameter_values.size(); u++) {
        for (int v = 0; v < current_parameter_values.size(); v++) {
          current_covariance_matrix(u, v) = new_covariance_matrix(u, v);
        }
      }
      cooldown = adapter_output["Cooldown"]; // reqired for the next iteration when proposal_SD_adapter is called again
      current_scaling_factor = adapter_output["New_Scaling_Factor"]; // required for the next iteraton when proposal_SD_adapter is called again

    }

    //////////////////////////////////////////////////////////////////////
    //                                                                  //
    //  4.  (optional) Print Out Various Useful MCMC Diagnostics        //
    //                                                                  //
    //////////////////////////////////////////////////////////////////////

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

