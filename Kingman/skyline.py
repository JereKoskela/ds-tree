import msprime
import tskit
import scipy.special
import numpy
import math
import sys

# recombination probability per site per generation
r = 1e-8

# Single LG simulation
positions = [0, 28e6]
rates = [r]
rate_map = msprime.RateMap(position = positions, rate = rates)

demo = msprime.Demography()
demo.add_population(name = "A", initial_size = 1.8e6)
demo.add_population_parameters_change(time = 5000, initial_size = 0.3e6)
demo.add_population_parameters_change(time = 35000, initial_size = 0.05e6)

# Sample size and population ploidy
# Haploid sample size is n * p
n = 68
p = 2

reps = 1
for i in range(0, reps):
    ts = msprime.sim_ancestry(samples = n,
                              ploidy = p,
                              demography = demo,
                              recombination_rate = rate_map)

    sfs = ts.allele_frequency_spectrum(mode='branch',
                                       span_normalise=False,
                                       polarised=True)
    print(sfs)
