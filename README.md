# Genetic diversity impacts climate-induced species range shifts

This repository contains code and data for the study entitled *Genetic diversity impacts climate-induced species range shifts*.

In this study, we test whether species' genetic diversity helps explain the rates at which species are shifting their geographic ranges in response to climate change.

## Instructions

See the `scripts/` folder for analysis pipelines and recreation of figures. The `data/` folder contains the preprocessed dataset.

## Modeling Framework

We modeled absolute range shift rates using generalized linear mixed-effects models (GLMMs) with a Gamma error distribution and log link. Fixed effects included genetic diversity, climate change velocity, and range position (trailing edge, centroid, leading edge), along with all two-way and three-way interactions. Five methodological covariates were included to account for heterogeneity among studies, and taxonomic class was included as a random intercept. Observations were weighted to balance contributions across range positions, and continuous predictors were standardized.

Because inference based on p-values in GLMMs is debated, we used a bootstrap approach for statistical inference. For each bootstrap iteration, observations were resampled with replacement, weights were recalculated, and the GLMM was refitted using maximum likelihood. From 10,000 converged bootstrap models, we estimated mean coefficients, 95% confidence intervals, and the significance of effects based on whether bootstrap distributions overlapped zero. Model fit was assessed using marginal R² (lognormal method), and hierarchical variance partitioning was used to quantify the relative contribution of predictors. Sensitivity analyses using range-position-specific models confirmed the robustness of results.

## 
