# cndc_bayesian_eval
Minimum example of Bayesian method to fit critical nitrogen dilution curve (CNDC) with intended to explore new methods to evaluate curves

## Use and Objective
This repo is intended to develop methods to allow for the comparison of critical N dilution curve fits using Bayesian methods and based on the previous work of [Makowksi et al. (2020)](http://doi.org/10.1016/j.eja.2020.126076). This repo contains a portion of a broader analysis used for working on a manuscript in progress.

## Files
- `Makowski-2020.pdf`: PDF of Makowksi et al. (2020) paper  
- `Makowski-2020-Appendix-D.R`: Manually transposed copy of R script from Appendix D of Makowksi et al. (2020)  
- `Giletto-2020.xlsx`: Data copied from Appendix A and B of [Giletto et al. (2020)](http://doi.org/10.1016/j.eja.2020.126114)  
- `get-data.R`: Script used to combine Giletto et al. (2020) with other to-be-publised data from Rosen Lab (e.g., Bohman data)  
- `data.csv`: Output from `get-data.R` with Giletto and Bohman data properly formatted for use in `fit-model.R`  
- `fit-model.R`: Script used to fit Bayesian hierarchical model to fit CNDC using framework of Makowksi et al. (2020)  
- `cndc-bayesian_eval.Rproj`: RStudio project for use in this analysis  

## Fitted Model
- run script `manuscript/analysis/model-get.R` to retrieve previously fitted model
- run script `manuscript/analysis/model-fit.R` to fit model locally

## Reproducible Environment
RStudio project now includes `renv` environment control. Note: this is currently incompatible with the `model-fit.R` script due to conflicts with `rstan` from file naming issues... In order to run `model-fit.R`, need to ensure that system packages match `renv.lock` file and disable `renv` with `renv::deactivate()`.
