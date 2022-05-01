#ifndef ANC
#define ANC

#include <algorithm>
#include <cmath>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_gamma.h>
#include <libconfig.h++>
#include <vector>

struct Ancestry {

  Ancestry(const int n)
      : sim_time(0), parent(), order(), active(), node_time() {
    setup(n);
  }

  void reset(const int n) {
    sim_time = 0;
    parent.clear();
    order.clear();
    node_time.clear();
    active.clear();
    setup(n);
    return;
  }

  void setup(const int n) {
    for (int j = 0; j < n; j++) {
      parent.push_back(-1);
      order.push_back(1);
      node_time.push_back(0);
      active.push_back(j);
    }
    return;
  }

  void simulate_merger(const int family_size, gsl_rng *gen) {
    int m = 0;
    int running_total = 0;
    int mergers = 0;
    int max_mergers = 1;
    int index, event_order;
    for (int i = 0; i < max_mergers; i++) {
      m = gsl_ran_binomial(gen, 1.0 / (max_mergers - i),
                           family_size - running_total);
      if (m > 1) {
        mergers++;
        event_order = 0;
        for (int j = 0; j < m; j++) {
          index = gsl_rng_uniform_int(gen, active.size());
          event_order += order[active[index]];
          parent[active[index]] = parent.size();
          active.erase(active.begin() + index);
        }
        node_time.push_back(sim_time);
        parent.push_back(-1);
        order.push_back(event_order);
      }
      running_total += m;
    }
    for (int i = 0; i < mergers; i++) {
      active.push_back(parent.size() - mergers + i);
    }
    return;
  }

  void fill_relative_branch_lengths(std::vector<double> &sfs) {
    double sfs_sum = 0;
    std::fill(sfs.begin(), sfs.end(), 0);
    for (unsigned int i = 0; i < parent.size() - 1; i++) {
      sfs[order[i] - 1] += node_time[parent[i]] - node_time[i];
      sfs_sum += node_time[parent[i]] - node_time[i];
    }
    for (unsigned int i = 0; i < sfs.size(); i++) {
      sfs[i] /= sfs_sum;
    }
    return;
  }

  void simulate(const double c, gsl_rng *gen, std::vector<double> &sfs) {
    double rate, u, tmp, y;
    int n = (int)active.size();
    int family_size = 0;
    while (n > 1) {
      rate = gsl_sf_choose(n, 2) * (1 + c / 2);
      sim_time += gsl_ran_exponential(gen, 1 / rate);
      if (gsl_rng_uniform(gen) < 1 / (1 + c / 2)) {
        family_size = 2;
      } else {
        y = sqrt(gsl_rng_uniform(gen));
        if (y > 1e-9) {
          u = (n - 1) * log(1 - y) + log(1 + (n - 1) * y);
          u = exp(log(1 - exp(u)) - 2 * log(y) - gsl_sf_lnchoose(n, 2));
        } else {
          u = 0;
          for (int j = 2; j <= n; j += 2) {
            tmp = (j - 1) * exp(gsl_sf_lnchoose(n, j) + (j - 2) * log(y));
            if (tmp / u < 1e-12) {
              break;
            }
            u += tmp;
          }
          for (int j = 3; j <= n; j += 2) {
            tmp = (j - 1) * exp(gsl_sf_lnchoose(n, j) + (j - 2) * log(y));
            if (tmp / u < 1e-12) {
              break;
            }
            u -= tmp;
          }
          u /= gsl_sf_choose(n, 2);
        }
        if (gsl_rng_uniform(gen) < u) {
          do {
            family_size = 2 + gsl_ran_binomial(gen, y, n - 2);
          } while (gsl_rng_uniform(gen) > 1 / gsl_sf_choose(family_size, 2));
        }
      }
      if (family_size > 1) {
        simulate_merger(family_size, gen);
      }
      n = (int)active.size();
      family_size = 0;
    }
    fill_relative_branch_lengths(sfs);
    return;
  }

  double sim_time;
  std::vector<int> parent, order, active;
  std::vector<double> node_time;
};

#endif
