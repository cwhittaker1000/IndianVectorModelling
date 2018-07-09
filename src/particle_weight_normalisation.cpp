// Specifying all the includes and depends required to run the model
#include "Particle_Weight_Normalisation.hpp"
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
std::vector <double> Particle_Weight_Normalisation(std::vector <double> particle_weights) {

  // Declaring a vector for the processed particle weights
  std::vector <double> particle_weights_normal_scale(particle_weights.size(), 0);
  std::vector <double> particle_weights_normal_scale_normalised(particle_weights.size(), 0);

  // Checking the passed values, and adapting them if they're out of relevant range
  for (int i = 0; i < particle_weights.size(); i++) {
    double weight = particle_weights[i];
    if (weight < -16) {
      weight = -16;
    }
    particle_weights_normal_scale[i] = exp(weight);
  }

  // Accumulate takes the initial point, the end point and the initial value you want to sum to, from and together.
  double summed_particle_weights = accumulate(particle_weights_normal_scale.begin(), particle_weights_normal_scale.end(), 0.0);

  // Loop to return a normalised set of particle weights, not on the log scale
  for (int i = 0; i < particle_weights_normal_scale.size(); i++) {
    particle_weights_normal_scale_normalised[i] = particle_weights_normal_scale[i] / summed_particle_weights;
  }

  return(particle_weights_normal_scale_normalised);
}
