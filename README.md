A PhD project titled "Mapping and predicting the risk of stillbirth and preterm birth in England"

This repository includes all code related to the project.

This repository contains:

Directory: 1_prediction_modelling

Scripts:

Directory: 2_spatial_multilevel_modelling

Scripts: 

inla_models.R - fitting Bayesian hierachical models (MSOA level) for preterm birth below 37 weeks gestation with and without spatial dependancies;
                output: 1 table with predictive performance and fit for 4 models (unadjusted and adjusted with IIDRE or SSRE); 
                output: 3 plots from unadjusted SSRE model (point estimates, significance, and 95%CI width)
