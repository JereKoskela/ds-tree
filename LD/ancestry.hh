#ifndef ANC
#define ANC

#include <cmath>
#include <cstdlib>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <iostream>
#include <libconfig.h++>
#include <vector>

struct Ancestry {

  Ancestry() : coal_time(2) {}

  void simulate(const double c, const double gap, const double s,
                const double a, gsl_rng *gen) {
    double rate = 2 + c;
    double sim_time = gsl_ran_exponential(gen, 1 / rate);
    int loc = gsl_rng_uniform_int(gen, 2);
    coal_time[loc] = sim_time;
    if (gsl_rng_uniform(gen) >
        c * (a + (1 - a) * exp(-2 * s * gap)) / (4 * rate)) {
      sim_time += gsl_ran_exponential(gen, 1 / (1 + c / 2));
    }
    coal_time[(loc + 1) % 2] = sim_time;
    return;
  }

  std::vector<double> coal_time;
};

#endif
