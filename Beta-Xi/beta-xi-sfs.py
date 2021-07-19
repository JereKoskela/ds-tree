import msprime
import tskit

n = 136 # Sample size
r = 1e-8 # Recombination rate
l = 35e6 # Genome length
a = 1.1 # alpha parameter

ts = msprime.sim_ancestry(samples = n / 2,
                          recombination_rate = r,
                          sequence_length = l,
                          model=msprime.BetaCoalescent(alpha = a))
sfs = ts.allele_frequency_spectrum(mode='branch',
                                   span_normalise=False,
                                   polarised=True)
print(sfs / sum(sfs))
