# ds-tree
This repository contains several tools for inference and model selection under various coalescent models aiming to describe shallow gene genealogies.
The central tools are

1. an implementation of the coalescent of [Durrett and Schweinsberg (2005)](https://www.sciencedirect.com/science/article/pii/S0304414905000608), or more specifically, the infinite chromosome length limit in their Example 2.4, and
2. an adaptive ABC pipeline as described by [Vihola and Franks (2020)](https://academic.oup.com/biomet/article/107/2/381/5721278) for inferring paramenters under the same model.

Further tools include

3. R scripts for post-processing ABC output and illustrating the resulting model fit,
4. a Python script for a corresponding ABC sampler under the diploid Beta-Xi-coalescent with exponential population growth, described in e.g. [Koskela and Wilke Berenguer (2019)](https://www.sciencedirect.com/science/article/pii/S0025556418303523), and
5. a Python script for simulating the Kingman coalescent with specified population size histories.

All aforementioned Python scripts use [msprime](https://tskit.dev/msprime/docs/stable/intro.html).

# Dependencies

ds-tree has been tested on Kubuntu 18.04 using the g++ 7.5.0 compiler.
To run, it requires a copy of the Gnu Scientific Library (tested on version 2.4).
Parts of the repository (specifically, items 4, and 5 above) also require msprime 1.0 or above (tested on version 1.0), along with its dependencies.
Item 3 above requires R (tested using RStudio 1.4.1106 with R version 4.1.0).

# Compilation
To compile the Durrett-Schweinsberg simulator and its corresponding ABC sampler, call `make` and `make sim` in the project root.
The other parts of the repository do not require compilation.
