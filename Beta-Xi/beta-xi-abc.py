import msprime
import tskit
import random
import math
import sys
import numpy
import scipy.special

# Desider number of MCMC steps.
steps = 10000

# Replace below with path to observed data file. It reads in the first line.
with open("../South-gl1-nsfs.txt", "r") as f:
    nsfs = numpy.array([float(i) for i in f.readline().strip().split()])
n = len(nsfs) + 1
# Small regulariser to avoid zeros in simulated SFSs
reg = min(nsfs) / 100
logit_sfs = numpy.array([math.log(i + reg) - math.log(1 - i - reg) for i in nsfs])

a_var = 0.1**2
g_var = 0.01**2
a = 1.5
g = 0.05

r = 1e-8
l = 35e6
acc = 0
accept_target = 0.1

a_lb = 1.01
a_ub = 1.99
g_lb = 0

for i in range(1, steps):
    a_new = random.gauss(a, math.sqrt(a_var))
    if a_new < a_lb:
        a_new = 2 * a_lb - a_new
    if a_new > a_ub:
        a_new = 2 * a_ub - a_new
    N = 1e7
    g_new = random.gauss(g, math.sqrt(g_var))
    if g_new < g_lb:
        g_new = 2 * g_lb - g_new
    demography = msprime.Demography.island_model([N], migration_rate = 0, growth_rate = [g_new])
    ts = msprime.sim_ancestry(samples = n / 2,
                              discrete_genome = False,
                              demography = demography,
                              recombination_rate = r,
                              sequence_length = l,
                              model=msprime.BetaCoalescent(alpha = a_new))
    sfs_new = ts.allele_frequency_spectrum(mode='branch',
                                           span_normalise=False,
                                           polarised=True)
    nsfs_new = sfs_new[1:n] / sum(sfs_new)
    logit_sfs_new = numpy.array([math.log(j + reg) - math.log(1 - j - reg) for j in nsfs_new])
    d = math.sqrt(numpy.dot(logit_sfs, logit_sfs_new))
    accepted = 0
    if i == 1:
        tol = d
        a_mean = (a + a_new) / 2
        g_mean = (g + g_new) / 2
    if d <= tol:
        a = a_new
        g = g_new
        acc = ((i - 1) * acc + 1) / i
        accepted = 1
    else:
        acc = (i - 1) * acc / i
    if i > 1:
        eta = 1 / i**(2/3)
        a_var_new = (1 - eta) * a_var + eta * (a - a_mean)**2
        g_var_new = (1 - eta) * g_var + eta * (g - g_mean)**2
        # Regularisation check to ensure positive proposal variances
        if a_var_new > 1e-4 and g_var_new > 1e-6:
            a_mean = (1 - eta) * a_mean + eta * a
            g_mean = (1 - eta) * g_mean + eta * g
            a_var = a_var_new
            g_var = g_var_new
        if i < steps / 2:
            tol = math.exp(math.log(tol) + eta * (accept_target - accepted))
    print(a, g)
