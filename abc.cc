#include "ancestry.hh"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

int main(int argc, char **argv) {
  if (argc != 3) {
    std::cout << "Call " << argv[0] << " <number of steps> <observed nsfs file>"
              << std::endl;
    return 1;
  }
  std::vector<double> obs_sfs;
  std::ifstream file(argv[2]);
  std::string line;
  std::getline(file, line);
  std::stringstream line_stream(line);
  double value;
  while (line_stream >> value) {
    obs_sfs.push_back(value);
  }
  int n = obs_sfs.size() + 1;
  for (int i = 0; i < n - 1; i++) {
    obs_sfs[i] = log(obs_sfs[i]) - log(1 - obs_sfs[i]);
  }
  gsl_rng *gen = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(gen, time(NULL));
  Ancestry anc(n);
  double c = 1;
  double c_mean = 0;
  double sd = 1;
  double c_new;
  std::vector<double> new_sfs(obs_sfs.size(), -1);
  std::vector<double> tmp_sfs(obs_sfs.size(), -1);
  double accept_rate = 0;
  double accept_target = 0.1;
  double d, tol, eta;
  int accepted;
  int reps = 1000;
  for (int i = 0; i < atoi(argv[1]); i++) {
    c_new = c + gsl_ran_gaussian_ziggurat(gen, sd);
    if (c_new < 0) {
      c_new *= -1.0;
    }
    std::fill(new_sfs.begin(), new_sfs.end(), 0);
    for (int j = 0; j < reps; j++) {
      anc.simulate(c_new, gen, tmp_sfs);
      anc.reset(n);
      for (int k = 0; k < n - 1; k++) {
        new_sfs[k] += tmp_sfs[k] / reps;
      }
    }
    for (int j = 0; j < n - 1; j++) {
      new_sfs[j] += log(new_sfs[j]) - log(1 - new_sfs[j]);
    }
    d = 0;
    for (int j = 0; j < n - 1; j++) {
      d += (obs_sfs[j] - new_sfs[j]) * (obs_sfs[j] - new_sfs[j]);
    }
    d = sqrt(d);
    if (i == 0) {
      tol = d;
      c_mean = (c + c_new) / 2.0;
    }
    if (d <= tol) {
      c = c_new;
      accept_rate = (i * accept_rate + 1.0) / (i + 1.0);
      accepted = 1;
    } else {
      accept_rate = i * accept_rate / (i + 1.0);
      accepted = 0;
    }
    if (i > 0) {
      eta = pow(i + 1, -2.0 / 3.0);
      sd = sqrt((1 - eta) * sd * sd + eta * (c - c_mean) * (c - c_mean));
      c_mean = (1 - eta) * c_mean + eta * c;
      if (i < atoi(argv[1]) / 2) {
        tol = exp(log(tol) + eta * (accept_target - accepted));
      }
    }
    std::cout << c << std::endl;
  }
  gsl_rng_free(gen);
  return 1;
}
