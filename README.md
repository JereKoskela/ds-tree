# ds-tree

This repository contains several tools for inference and model selection under various coalescent models aiming to describe shallow gene genealogies.
The central tools are

1. an implementation of the coalescent of [Durrett and Schweinsberg (2005)](https://www.sciencedirect.com/science/article/pii/S0304414905000608), or more specifically, the infinite chromosome length limit in their Example 2.4, and
2. an adaptive ABC pipeline as described by [Vihola and Franks (2020)](https://academic.oup.com/biomet/article/107/2/381/5721278) for inferring paramenters under the same model.

Further tools include

3. R scripts for post-processing ABC output and illustrating the resulting model fit,
4. a Python script for a corresponding ABC sampler under the diploid Xi-Beta-coalescent with exponential population growth, described in e.g. [Koskela and Wilke Berenguer (2019)](https://www.sciencedirect.com/science/article/pii/S0025556418303523),
5. a Python script for simulating the Kingman coalescent with specified population size histories, and
6. an adaptive ABC pipeline for linkage disequilibrium under a two-locus generalisation of the [Durrett and Schweinsberg (2005)](https://www.sciencedirect.com/science/article/pii/S0304414905000608) coalescent, for sample size 2 only.

All aforementioned Python scripts use [msprime](https://tskit.dev/msprime/docs/stable/intro.html).

# Dependencies

ds-tree has been tested on Kubuntu 18.04 using the g++ 7.5.0 compiler.
To run, it requires a copy of the Gnu Scientific Library (tested on version 2.4).
Parts of the repository (specifically, items 4, and 5 above) also require msprime 1.0 or above (tested on version 1.0), along with its dependencies.
Item 3 above requires R (tested using RStudio 1.4.1106 with R version 4.1.0).

# Compilation

To compile the Durrett-Schweinsberg simulator and its corresponding ABC sampler, call `make` and `make sim` in the project root.
Compilation should take no more than a few seconds.
The LD tool also requires complilation by calling `make` in the LD folder.
The other parts of the repository do not require compilation.

# Tutorial

## Durrett-Schweinsberg coalescent simulator

This simulator returns the logit-transformed mean of a specified number of independent realisations of normalised site frequency spectra under the Durrett-Schweinsberg coalescent with a given sample size and parameter `c` as described in Example 2.4 of [Durrett and Schweinsberg (2005)](https://www.sciencedirect.com/science/article/pii/S0304414905000608).
Run the simulator by compiling, and then calling 
```
./simulate n c r
```
where `n` is the sample size, `c` is the parameter, and `r` is the desired number of replicates.
See included `ds-simulate-spectra.sh` shell script for examples.
The combined serial runtime of `ds-simulate-spectra.sh` on a typical mid-range desktop is around two minutes. 

## Durrett-Schweinsberg ABC

This simulator implements the adaptive ABC MCMC method of [Vihola and Franks (2020)](https://academic.oup.com/biomet/article/107/2/381/5721278) for sampling the posterior distribution of the parameter `c` for a given normalised site frequency spectrum, with an improper uniform prior.
Run the simulator by compiling, and then calling
```
./abc n sfs
```
where `n` is the desired number of MCMC steps and `sfs` is the location of a text file containing the observed normalised spectrum.
The call returns a sampled value of `c` per MCMC step.
Four example spectra (South-gl1-nsfs.txt, South-gl2-nsfs.txt, Thistilfj-gl1-nsfs.txt, and Thistilfj-gl2-nsfs.txt) are provided, along with the `ds-abc.sh` shell script with an example simulation call for each spectrum.
The serial runtime of the script on a mid-range desktop is around 5 hours, with each individual example taking between 60 and 100 minutes. 

**Warning** The numbers of MCMC steps in these examples have been set to 10 000 for the purposes of illustration in a reasonable amount of CPU time.
For throroughly convergent results, we recommend increasing the number of steps to 100 000. 

## Durrett-Schweinsberg postprocessing

The R scripts `ds-abc-graphs.R` and `ds-sfs-graphs.R` are postprocessing scripts for visualising the output from `ds-abc.sh` and `ds-simulate-spectra.sh`.
They assume that the working directory of your R client is set to the location of the output from the `ds-simulate-spectra.sh` and`ds-abc.sh` shell scripts.
These scripts run in no more than a few seconds each.

## Xi-Beta ABC

The `beta-xi-abc.py` Python script implements the ABC MCMC method of [Vihola and Franks (2020)](https://academic.oup.com/biomet/article/107/2/381/5721278) for the Xi-Beta`(2 - a, a)`-coalescent with exponential population growth.
The object of inference is the tuple `(a, g)`, where the former governs the degree of skewness of the family size distribution in the Beta-coalescent, and `g` is the population-rescaled rate of exponential growth.
The (improper) prior distribution is uniform on the 2d positive orthant.
The location of an observed normalised site frequency spectrum needs to be specified on line 13 of the script.
Simulated site frequency spectra are generated using `msprime`.
Using the provided settings, `beta-xi-abc.py` runs in around 4 hours on a midrange desktop.

**Warning** The numbers of MCMC steps in these examples have been set to 10 000 for the purposes of illustration in a reasonable amount of CPU time.
For throroughly convergent results, we recommend increasing the number of steps to 100 000. 

## Beta-Xi simulator

The `beta-xi-sfs.py` is a Python script for simulating a normalised site frequency spectrum under the Beta-Xi-coalescent with exponential population growth using `msprime`.
With the present parameter values, it runs in under a second. 

## Kingman coalescent simulator 

The `skyline.py` and `smcpp.py` Python scripts call `msprime` to simulate site frequency spectra under the Kingman coalescent with respective population size histories specified in the scripts.
Approximate runtimes for these scripts on a mid-range desktop are 7 minutes (`skyline.py`) and 3 minutes (`smcpp.py`).

## Durrett-Schweinsberg LD ABC

The run the ABC pipeline, call `./simulate n` after compilation, where `n` is the desired number of MCMC steps.
The `sim.sh` script contains an example.
Initial values of parameters are hard-coded in `sim.cc`, and the data files `gap.txt` and `obs_r_squared.txt` contain examples of observed data files.
The former specifies desired positions along a genome in physical coordinates (the first position is always at zero), and the latter specifies corresponding observed $r^2$ values for each pairwise comparison between the site at the origin, and the site each given distance away.
The ABC pipeline attempts to match all provided values simultaneously by simulating pairwise realisations of the origin and each provided distance, though the weight given to higher distances decays as an inverse square.

The output of the ABC run is an `n` by `d + 4` array, where `d` is the number of provided distances.
Each row contains fitted values of all `d` $r^2$ values, preempted by four fitted parameter values: `c` as described above, `β / s` in the notation of [Durrett and Schweinsberg (2005)](https://www.sciencedirect.com/science/article/pii/S0304414905000608), a mutation rate `θ`, and a fraction of sweep events which affect the whole chromosome regardless of recombination.

The `R` script `draw_recomb.R` plots histograms, diagnostic trace plots, and fitted $r^2$ profiles from stored ABC output.
