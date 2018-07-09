#include "weighted_sampling_with_replacement.hpp"
#include <Rcpp.h>
#include <vector>
#include <string>
#include <iostream>
// [[Rcpp::plugins(cpp11)]]

int rbinomial1(int trials, double p) {
  return(R::rbinom(trials,p));
}

// [[Rcpp::export]]
std::vector<int> weighted_sampling(int n, std::vector<double> prob, double p_tot, int K)
{
  int k;
  int draw;
  std::vector<int> output;

  for (k = 0; k < K - 1; k++) {

    draw = rbinomial1(n, (prob[k] / p_tot));
    n -= draw;
    while (draw > 0) {
      output.emplace_back(k);
      draw--;
    }

    if (n <= 0) {
      return(output);
    }
    else {
      p_tot -= prob[k];
    }
  }
  while (n > 0) {
    output.emplace_back(K - 1);
    n--;
  }
  return(output);
}
