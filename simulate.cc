#include "ancestry.hh"
#include <cmath>
#include <cstdlib>
#include <gsl/gsl_rng.h>
#include <iostream>
#include <vector>

int main(int argc, char **argv) {
  if (argc != 4) {
    std::cout << "Call " << argv[0]
              << " <sample size> <c> <number of realisations>" << std::endl;
    return 1;
  }
  int n = atoi(argv[1]);
  std::vector<double> sfs(n - 1, 0);
  std::vector<double> tmp_sfs(n - 1, 0);
  gsl_rng *gen = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(gen, time(NULL));
  Ancestry anc(n);
  double c = atof(argv[2]);
  int reps = atoi(argv[3]);
  for (int j = 0; j < reps; j++) {
    anc.simulate(c, gen, tmp_sfs);
    anc.reset(n);
    for (int k = 0; k < n - 1; k++) {
      sfs[k] += tmp_sfs[k] / reps;
    }
  }
  for (int j = 0; j < n - 1; j++) {
    sfs[j] = log(sfs[j]) - log(1 - sfs[j]);
    std::cout << sfs[j] << " ";
  }
  std::cout << std::endl;
  gsl_rng_free(gen);
  return 1;
}
