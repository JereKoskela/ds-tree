#include "ancestry.hh"
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <gsl/gsl_rng.h>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

double log_sum_exp(const std::vector<double> &v) {
  double ret = 0;
  double max_elem = *std::max_element(v.begin(), v.end());
  for (unsigned int i = 0; i < v.size(); i++) {
    ret += exp(v[i] - max_elem);
  }
  ret = max_elem + log(ret);
  return ret;
}

double project(const double c, const double lb, const double ub) {
  double ret = c;
  while (ret < lb || ret > ub) {
    if (ret < lb) {
      ret = 2 * lb - ret;
    } else {
      ret = 2 * ub - ret;
    }
  }
  return ret;
}

int main(int argc, char **argv) {
  if (argc != 2) {
    std::cout << "Call " << argv[0] << " <number of steps>" << std::endl;
    return 1;
  }
  int steps = atoi(argv[1]);
  std::vector<double> obs_r_squared, gap;
  std::ifstream file("obs_r_squared.txt");
  std::string line;
  std::getline(file, line);
  std::stringstream line_stream(line);
  double value;
  while (line_stream >> value) {
    obs_r_squared.push_back(value);
  }
  std::ifstream gap_file("gap.txt");
  std::getline(gap_file, line);
  std::stringstream gap_line_stream(line);
  while (gap_line_stream >> value) {
    gap.push_back(log(value));
  }
  unsigned int n = obs_r_squared.size();
  Ancestry anc;
  gsl_rng *gen = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(gen, time(NULL));
  double c = 12;
  double s = 10;
  double m = 1e3;
  double a = 0.1;
  double cMean = 0;
  double sMean = 0;
  double mMean = 0;
  double aMean = 0;
  double sd_c = c / 10;
  double sd_s = s / 10;
  double sd_m = m / 10;
  double sd_a = a / 10;
  double cNew, sNew, mNew, aNew;
  double accept_rate = 0;
  double accept_target = 0.1;
  double d, tol = 0, eta;
  int accepted;
  int reps = 10000;
  std::vector<double> new_r_squared(n);
  std::vector<double> old_r_squared(n);
  std::vector<double> cross_term(n);
  std::vector<double> mean_one(n);
  std::vector<double> mean_two(n);
  std::vector<double> sd_one(n);
  std::vector<double> sd_two(n);
  std::vector<double> log_cross_term(reps);
  std::vector<double> log_mean_one(reps);
  std::vector<double> log_mean_two(reps);
  std::vector<double> log_sd_one(reps);
  std::vector<double> log_sd_two(reps);
  for (unsigned int k = 0; k < n; k++) {
    for (int j = 0; j < reps; j++) {
      anc.simulate(c, gap[k], s, a, gen);
      log_cross_term[j] = -m * (anc.coal_time[0] + anc.coal_time[1]);
      log_mean_one[j] = -m * anc.coal_time[0];
      log_mean_two[j] = -m * anc.coal_time[1];
      log_sd_one[j] =
          -m * anc.coal_time[0] + log(1 - exp(-m * anc.coal_time[0]));
      log_sd_two[j] =
          -m * anc.coal_time[1] + log(1 - exp(-m * anc.coal_time[1]));
    }
    cross_term[k] = log_sum_exp(log_cross_term);
    mean_one[k] = log_sum_exp(log_mean_one);
    mean_two[k] = log_sum_exp(log_mean_two);
    sd_one[k] = log_sum_exp(log_sd_one);
    sd_two[k] = log_sum_exp(log_sd_two);
  }
  for (unsigned int j = 0; j < n; j++) {
    old_r_squared[j] =
        (exp(cross_term[j]) - exp(mean_one[j] + mean_two[j] - log(reps))) /
        exp((sd_one[j] + sd_two[j]) / 2);
  }

  for (int i = 0; i < steps; i++) {
    cNew = c + gsl_ran_gaussian_ziggurat(gen, sd_c);
    cNew = project(cNew, 10, 25);
    sNew = s + gsl_ran_gaussian_ziggurat(gen, sd_s);
    sNew = project(sNew, 0, 1e4);
    mNew = m + gsl_ran_gaussian_ziggurat(gen, sd_m);
    mNew = project(mNew, 0, 1e4);
    aNew = a + gsl_ran_gaussian_ziggurat(gen, sd_a);
    aNew = project(aNew, 0, 1);
    for (unsigned int k = 0; k < n; k++) {
      for (int j = 0; j < reps; j++) {
        anc.simulate(cNew, gap[k], sNew, aNew, gen);
        log_cross_term[j] = -mNew * (anc.coal_time[0] + anc.coal_time[1]);
        log_mean_one[j] = -mNew * anc.coal_time[0];
        log_mean_two[j] = -mNew * anc.coal_time[1];
        log_sd_one[j] =
            -mNew * anc.coal_time[0] + log(1 - exp(-mNew * anc.coal_time[0]));
        log_sd_two[j] =
            -mNew * anc.coal_time[1] + log(1 - exp(-mNew * anc.coal_time[1]));
      }
      cross_term[k] = log_sum_exp(log_cross_term);
      mean_one[k] = log_sum_exp(log_mean_one);
      mean_two[k] = log_sum_exp(log_mean_two);
      sd_one[k] = log_sum_exp(log_sd_one);
      sd_two[k] = log_sum_exp(log_sd_two);
    }
    for (unsigned int j = 0; j < n; j++) {
      new_r_squared[j] = exp(cross_term[j] - (sd_one[j] + sd_two[j]) / 2) -
                         exp(mean_one[j] + mean_two[j] - log(reps) -
                             (sd_one[j] + sd_two[j]) / 2);
    }
    d = 0;
    for (unsigned int j = 0; j < n; j++) {
      d += (obs_r_squared[j] - new_r_squared[j]) *
           (obs_r_squared[j] - new_r_squared[j]) / pow(gap[j], 2);
    }
    d = sqrt(d);
    if (i == 0) {
      tol = d;
      cMean = (c + cNew) / 2;
      sMean = (s + sNew) / 2;
      mMean = (m + mNew) / 2;
      aMean = (a + aNew) / 2;
    }
    if (d <= tol) {
      c = cNew;
      s = sNew;
      m = mNew;
      a = aNew;
      accept_rate = (i * accept_rate + 1) / (i + 1);
      accepted = 1;
      for (unsigned int j = 0; j < n; j++) {
        old_r_squared[j] = new_r_squared[j];
      }
    } else {
      accept_rate = i * accept_rate / (i + 1);
      accepted = 0;
    }
    if (i > 0) {
      eta = pow(i + 1, -2.0 / 3.0);
      sd_c = sqrt((1 - eta) * sd_c * sd_c + eta * (c - cMean) * (c - cMean));
      sd_s = sqrt((1 - eta) * sd_s * sd_s + eta * (s - sMean) * (s - sMean));
      sd_m = sqrt((1 - eta) * sd_m * sd_m + eta * (m - mMean) * (m - mMean));
      sd_a = sqrt((1 - eta) * sd_a * sd_a + eta * (a - aMean) * (a - aMean));
      cMean = (1 - eta) * cMean + eta * c;
      sMean = (1 - eta) * sMean + eta * s;
      mMean = (1 - eta) * mMean + eta * m;
      mMean = (1 - eta) * aMean + eta * a;
      if (i < atoi(argv[1]) / 2) {
        tol = exp(log(tol) + eta * (accept_target - accepted));
      }
    }
    std::cout << c << " " << s << " " << m << " " << a;
    for (unsigned int i = 0; i < n; i++) {
      std::cout << " " << old_r_squared[i];
    }
    std::cout << std::endl;
  }
  gsl_rng_free(gen);
  return 1;
}
