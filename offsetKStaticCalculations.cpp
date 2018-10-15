if (rainfall_relationship == "mean") {

  // Calculating Raw Rainfall Sum
  std::vector<double>::const_iterator first = rF.begin() + t + offset - tau_with_dt_static;
  std::vector<double>::const_iterator last = rF.begin() + t + offset;
  std::vector<double> rFx_static(first, last);
  rFsum = std::accumulate(rFx_static.begin(), rFx_static.end(), 0.0);
  // Calculating Whether Washout Occurs
  double rainfall_average_static_calc = rFsum / tau_with_dt_static;
  if (rainfall_average_static_calc > Washout_Threshold) {
    K_Static = 0;
    marker = 0;
  }
  else {
    if (marker == 0) { // First Time Washout Stops, Carrying Capacity Is At Max
      K_Static = K_Max_Static;
      marker_stop = t;
      marker = 1;
    }
    else { // Then Carrying Capacity Begins to Decline
      if (decline_type == "exponential") {
        double calculation = -washout_exp_scaling_factor * (t - marker_stop);
        K_Static = K_Max_Static * exp(calculation);
      }
      else {
        K_Static = Hill_Function((t - marker_stop), K_Static, washout_hill_one, washout_hill_two);
      }
    }
  }
  k_static_output[t] = K_Static;
}

else if (rainfall_relationship == "linear") {
  // Calculating Linearly Weighted Rainfall Sum
  double rFsum = 0.0;
  int counter = 1;
  for (int r = (t + offset- tau_with_dt_static + 1); r <= t + offset; r++) {
    rFsum = rFsum + (((t - tau_with_dt_static + counter) - t + tau_with_dt_static) * rF[r]);
    counter = counter + 1;
  }
  // Calculate Whether Washout Occurs
  double rainfall_average = scaling_factor * ((2.0 / pow(tau_with_dt_static, 2)) * rFsum); // unsure if this is what I want, double check
  if (rainfall_average > Washout_Threshold) {
    K_Static = 0;
    marker = 0;
  }
  else {
    if (marker == 0) { // First Time Washout Stops, Carrying Capacity Is At Max
      K_Static = K_Max_Static;
      marker_stop = t;
      marker = 1;
    }
    else { // Then Carrying Capacity Begins to Decline
      if (decline_type == "exponential") {
        double calculation = -washout_exp_scaling_factor * (t - marker_stop);
        K_Static = K_Max_Static * exp(calculation);
      }
      else {
        K_Static = Hill_Function((t - marker_stop), K_Static, washout_hill_one, washout_hill_two);
      }
    }
  }
  k_static_output[t] = K_Static;
}

else if (rainfall_relationship == "exponential") {
  // Calculating Exponentially Weighted Rainfall Sum
  double temp_tau = tau_with_dt_static;
  double rFsum = 0.0;
  double calc;
  for (int r = 0; r <= t + offset; r++) {
    calc = (-(t + offset - r))/temp_tau;
    rFsum = rFsum + exp(calc) * rF[r];
  }
  // Calculating Whether Washout Occurs
  double rainfall_average = scaling_factor * (1.0 / (temp_tau * (1 - exp(-(t + offset + 1) / temp_tau)))) * rFsum;
  if (rainfall_average > Washout_Threshold) {
    K_Static = 0;
    marker = 0;
  }
  else {
    if (marker == 0) { // First Time Washout Stops, Carrying Capacity Is At Max
      K_Static = K_Max_Static;
      marker_stop = t;
      marker = 1;
    }
    else { // Then Carrying Capacity Begins to Decline
      if (decline_type == "exponential") {
        double calculation = -washout_exp_scaling_factor * (t - marker_stop);
        K_Static = K_Max_Static * exp(calculation);
      }
      else {
        K_Static = Hill_Function((t - marker_stop), K_Static, washout_hill_one, washout_hill_two);
      }
    }
  }
  k_static_output[t] = K_Static;
}
