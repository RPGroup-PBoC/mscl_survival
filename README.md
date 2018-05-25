# `mscl_survival` #

This repository contains all Python code, Jupyter Notebooks, small data files, and the full text for the publication "Connecting the dots between mechanosensitive channel abundance, osmotic shock, and survival at single-cell resolution."  We encourage the scientific community to fork this repository and open constructive issues that will improve and/or correct any writing in this work. The general architecture of the repository is described below. 


## `code` ## 

This repository contains all processing, analysis, and figure code used in this work. Within `code`, there five subfolders separated by document type. 

1. **`analysis`**
    Contains Python scripts (`.py`) that conduct the major analysis described in the work. This includes all Markov chain Monte Carlo analyses. 

2. **`figs`**
    Contains Python scripts (`.py`) used to generate the main text and supplementary information figures. 
3. **`notebooks`**
    Contains a Jupyter notebook that walks through the image analysis pipeline. 
4. **`processing`**
    Contains all experimental processing. Each folder within `processing` is separated by date, strain, and shock rate. Each folder contains a Python script `.py` that performs all image processing and saves a `.csv` for that experimental run. 
5. **`stan`**
    All `stan` models which were used in the MCMC. 

## `data/csv` ##

This folder contains a variety of `.csv` files generated from the MCMC analysis. The primary file, `.mscl_survival_stats.csv`, contains all measurements that were used in the formal analysis.

## `doc` ##

Herein lies the main text of the paper. It was written in MarkDown syntax (`.md`) and converted to LaTex, PDF, and HTML through [Pandoc](https://pandoc.org/). The file, `make.py` builds the entire paper using information provided in `headers/build_main.yaml`.

## `figs` ##

This folder contains a select number of the figures generated in this work. 

## `mscl` ##

This is a Python module which contains all invariant functions used for processing, analysis, and plotting. It is structured as a typical Python module that can be imported. There are several components to this module, described below. 

1. **`image.py`**
    Contains functions related to the loading, manipulation, and processing of the raw microscopy images and parsing the `.xml` files produced through manual identification of survivors and fatalities.

2. **`mcmc.py`**
    A set of functions that manipulate the output from `pystan` sampling.

3.**`plotting.py`**
    Functions for generating standard plots and setting the plotting style.

5.**`stats.py`**
    A set of functions that compute numerous frequently used statistics.


## Contact ##
Any questions regarding this repository should be directed to Griffin Chure through the issue tracker. Please do not send private messages or emails. 

## License ##
This work is licensed under a [Creative Commons Attribution CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/). All code contained herein is licenced under an [MIT license]() which is as follows

```
	(c) Copyright 2018 Griffin Chure, Heun Jin Lee, and Rob Phillips

	Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

	The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
```
