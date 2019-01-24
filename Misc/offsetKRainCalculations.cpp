if (rainfall_relationship == "mean") {

  if (rainfall_effect == "raw") {
    // Calculating Raw Rainfall Average
    std::vector<double>::const_iterator first = rF.begin() + t + offset- tau_with_dt_rain;
    std::vector<double>::const_iterator last = rF.begin() + t + offset;
    std::vector<double> rFx(first, last);
    rFsum = std::accumulate(rFx.begin(), rFx.end(), 0.0); // Don't forget to do 0.0!!! Otherwise accumulator will produce an int!!
    K_Rain = scaling_factor * (1.0 + ((1.0 / tau_with_dt_rain) * rFsum));
    k_rain_output[t] = K_Rain;
  }
  else if (rainfall_effect == "hill") {
    // Calculating Hill Function'd Raw Average Rainfall
    std::vector<double>::const_iterator first = rF.begin() + t + offset- tau_with_dt_rain;
    std::vector<double>::const_iterator last = rF.begin() + t + offset;
    std::vector<double> rFx(first, last);
    rFsum = std::accumulate(rFx.begin(), rFx.end(), 0.0); // Don't forget to do 0.0!!! Otherwise accumulator will produce an int!!
    double mean_rainfall = 1.0 + ((1.0 / tau_with_dt_rain) * rFsum);
    double hill_output = Hill_Function(mean_rainfall, K_Max_Hill_Rainfall, Hill_Rainfall_1, Hill_Rainfall_2);
    K_Rain = scaling_factor * hill_output;
    k_rain_output[t] = K_Rain;
  }
}

else if (rainfall_relationship == "linear") {

  if (rainfall_effect == "raw") {
    // Calculating Linearly Weighted Average Rainfall
    double rFsum = 0.0;
    int counter = 1;
    for (int r = (t + offset - tau_with_dt_rain + 1); r <= t + offset; r++) {
      rFsum = rFsum + (((t - tau_with_dt_rain + counter) - t + tau_with_dt_rain) * rF[r]);
      counter = counter + 1;
    }
    K_Rain = scaling_factor * ((2.0 / pow(tau_with_dt_rain, 2)) * rFsum);
    k_rain_output[t] = K_Rain;
  }

  else if (rainfall_effect == "hill") {
    // Calculating Hill Function'd, Linearly Weighted Average Rainfall
    double rFsum = 0.0;
    int counter = 1;
    for (int r = (t + offset - tau_with_dt_rain + 1); r <= t + offset; r++) {
      rFsum = rFsum + (((t - tau_with_dt_rain + counter) - t + tau_with_dt_rain) * rF[r]);
      counter = counter + 1;
    }
    double mean_rainfall = ((2.0 / pow(tau_with_dt_rain, 2)) * rFsum);
    double hill_output = Hill_Function(mean_rainfall, K_Max_Hill_Rainfall, Hill_Rainfall_1, Hill_Rainfall_2);
    K_Rain = scaling_factor * hill_output;
    k_rain_output[t] = K_Rain;
  }
}

else if (rainfall_relationship == "exponential") {

  if (rainfall_effect == "raw") {
    // Calculating Exponentially Weighted Rainfall Average
    double temp_tau = tau_with_dt_rain;
    double rFsum = 0.0;
    double calc;
    for (int r = 0; r <= t + offset; r++) { // need to decide whether to start at beginning of rain or t = 0 ??
      calc = (-(t + offset - r))/temp_tau;
      rFsum = rFsum + exp(calc) * rF[r];
    }
    K_Rain = scaling_factor * (1.0 / (temp_tau * (1 - exp(-(t + offset + 1) / temp_tau)))) * rFsum;
    k_rain_output[t] = K_Rain;
  }
  else if (rainfall_effect == "hill") {
    // Calculating Hill Function'd Exponentially Weighted Rainfall Average
    double temp_tau = tau_with_dt_rain;
    double rFsum = 0.0;
    double calc;
    for (int r = 0; r <= t + offset; r++) {
      calc = (-(t + offset - r))/temp_tau;
      rFsum = rFsum + exp(calc) * rF[r];
    }
    double mean_rainfall = (1.0 / (temp_tau * (1 - exp(-(t + offset + 1) / temp_tau)))) * rFsum;
    double hill_output = Hill_Function(mean_rainfall, K_Max_Hill_Rainfall, Hill_Rainfall_1, Hill_Rainfall_2);
    K_Rain = scaling_factor * hill_output;
    k_rain_output[t] = K_Rain;
  }
}
