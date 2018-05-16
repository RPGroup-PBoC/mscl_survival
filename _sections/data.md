---
title: Data & Analysis
order: 3
icon: fa-database
---

### `.csv` files
All single-cell measurements of the osmotic-shock experiments can be downloaded as a `.csv` here:
* [**`mscl_survival_data.csv`**](mscl_survival_data.csv)

Raw images and MCMC chains are hosted on the CaltechDATA data repository and can be downloaded through the links below.

* [**Raw Image Data** (DOI: XXX)]()
* [**MCMC Sampling Chains (`.csv`)** (DOI: XXX)]()

All Python code for processing the raw data sets can be reached through the `master` branch of this repository, [`github.com/rpgroup-pboc/mscl_survival.git`](http://www.github.com/rpgroup-pboc/mscl_survival.git). All MCMC was performed using the [Stan probabailistic programming language.](www.mc-stan.org) The `.stan` models used for the hierarchical logistic regression can be downloaded below:

* [**`hierarchical_calibration_factor.stan`**](code/stan/hierarchical_calbration_factor.stan)
* [**`hierarchical_logistic_regression.stan`**](code/stan/hierarchical_logistic_regression.stan)
* [**`complete_analysis.stan`**](code/stan/complete_analysis.stan) **\|** This model samples both hierarchical models simultaneously. 

We have also included a Jupyter Notebook which describes the image processing procedure undertaken in this work. 

* [**Appendix A:** Image Processing Pipeline](code/notebooks/appendix_A_image_processing.html) **\|** [download as `.ipynb`](code/notebooks/appendix_A_image_processing.ipynb)