# Quantifying critical N dilution curves across G × E × M effects for potato using a partially-pooled Bayesian hierarchical method

## Overview
This repoistory contains the data, analysis, and manuscript version history of the Bohman et al. (2022) manuscript in progress: "Quantifying critical N dilution curves across G × E × M effects for potato using a partially-pooled Bayesian hierarchical method".

### Authors
* Brian J. Bohman (@bohm0072)
* Michael J. Culshaw-Maurer (@MCMaurer)
* Feriel Ben Abdallah
* Claudia Giletto
* Gilles Bélanger
* Fabián G. Fernández
* Yuxin Miao
* David J. Mulla
* Carl J. Rosen

### Abstract
Multiple critical N dilution curves [CNDCs] have been previously developed for potato; however, attempts to directly compare differences in CNDCs across genotype [G], environment [E], and management [M] interactions have been confounded by nonuniform statistical methods, biased experimental data, and lack of proper quantification of uncertainty in the critical N concentration [%N<sub>c</sub>]. This study implements a partiallypooled Bayesian hierarchical method to develop CNDCs for previously published and newly reported experimental data, systematically evaluates the difference in %N<sub>c</sub> [∆%N<sub>c</sub>] across G × E × M effects, and directly compare CNDCs from the Bayesian framework to CNDCs from conventional statistical methods. The partially-pooled Bayesian hierarchical method implemented in this study has the advantage of being less susceptible to inferential bias at the level of individual G × E × M interactions compared to alternative statistical methods that result from insufficient quantity and quality of experimental datasets (e.g., unbalanced distribution of N limiting and non-N limiting observations). This method also allows for a direct statistical comparison of differences in %N<sub>c</sub> across levels of the G × E × M interactions. Where found to be significant, ∆%N<sub>c</sub> was attributed to variation in the timing of tuber initiation (e.g., maturity class) and the relative rate of tuber bulking (e.g., planting density) across G x E × M interactions. In addition to using the median value for %N<sub>c</sub> (i.e., CNDC), the lower and upper boundary values for the credible region (i.e., CNDC<sub>lo</sub> and CNDC<sub>up</sub>) derived using the Bayesian framework should be used in calculation of N nutrition index (and other calculations) to account for uncertainty in %N<sub>c</sub>. Overall, this study provides additional evidence that %N<sub>c</sub> is dependent upon G × E × M interactions; therefore, evaluation of crop N status or N use efficiency must account for variation in %N<sub>c</sub> across G × E × M interactions.

## Organization

### `manuscript`
This directory contains all materials used to develop the compiled manuscript, including analysis scripts, draft manuscript versions, figure image files, stored model objects, and csv formatted tables

### `data`
This directory contains all source data, scripts to format source data, and analysis-ready data used for subsequent analysis in the `manuscript` directory

### `literature review`
This directory contains a small sub-set of previously published work on this topic, including published origin of source data and reference methods.

### `notes`
This directory contains various notes.

### `prelim analysis`
This directory contains a set of minimal examples demonstrating various analysis methods, including data visualization and a comparison of model fit using `brms` vs. `JAGS`

### `renv`
This directory contains the `renv` files used to reproduce this analysis.  Note: this is currently incompatible with the `manuscript/analysis/model-fit.R` script due to conflicts with `rstan` from file naming issues... In order to run `manuscript/analysis/model-fit.R`, need to ensure that system packages match `renv.lock` file and disable `renv` with `renv::deactivate()`.
