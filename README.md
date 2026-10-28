# Biodiversity-ecosystem functioning relationships in Mediterranean fish communities

This repository contains the R code used to investigate biodiversity-ecosystem functioning (BEF) relationships in fish communities across four northern Mediterranean sub-regions: Northern Spain, the Gulf of Lion, the Adriatic Sea and the Aegean Sea. It analyses 23 years (1999–2021) of Mediterranean International Trawl Survey (MEDITS) data (10,034 hauls, 331 species) with piecewise structural equation modelling, relating environmental conditions and anthropogenic pressures to taxonomic diversity, functional diversity, community-weighted maximum body length and fish biomass. Fish biomass is treated as a biomass-related proxy for ecosystem functioning.

## Purpose
- Compute taxonomic, functional and trait-based diversity metrics for each trawl haul
- Assemble environmental and anthropogenic covariates at haul positions
- Select and fit the sub-models of a piecewise SEM with spatial correlation structures
- Summarise direct, indirect and total effects and run sensitivity analyses

## Methods and Workflow
Scripts are numbered in execution order and must be run from the project root.

1. **Diversity metrics**: taxonomic, functional (mFD) and community-weighted trait metrics per haul ([1-Diversity_metrics.r](1-Diversity_metrics.r)).
2. **Covariates**: Copernicus Marine environmental variables, fishing pressure and gravity extracted at haul positions ([2-Environmental_anthropogenic_covariatesR](2-Environmental_anthropogenic_covariatesR)).
3. **SEM dataset**: merging, filtering and collinearity screening ([3-Build_dataset.R](3-Build_dataset.R)).
4. **Sub-model selection** (spatial correlation structure, then `MuMIn::dredge()`), one script per response:
   [4a](4a-Build_submodel_taxo_richness.R), [4b](4b-Build_submodel_taxo_evenness.R), [4c](4c-Build_submodel_log_cwm_length.R), [4d](4d-Build_submodel_fric.R), [4e](4e-Build_submodel_feve.R), [4f](4f-Build_submodel_biomass.R).
5. **Residual spatial autocorrelation**: Moran's I of the final sub-models ([5-Check_spatial_autocor.R](5-Check_spatial_autocor.R)).
6. **Piecewise SEM**: final model fitted to all GSAs ([6-Build_SEM_model_ALL_GSA_clean.R](6-Build_SEM_model_ALL_GSA_clean.R)).
7. **Sensitivity analyses**: alternative biomass specifications and spatial block cross-validation ([7-Sensitivity_analyses.R](7-Sensitivity_analyses.R)).
8. **Summary statistics**: correlations and among-GSA comparisons ([8-Exploratory_analyses.R](8-Exploratory_analyses.R)).
9. **Figures**: [Figure_1.R](Figure_1.R) to [Figure_6.R](Figure_6.R) and [Figure_S1.R](Figure_S1.R) to [Figure_S5.R](Figure_S5.R).

## Inputs and Data
- [Data/](Data): MEDITS survey data (`MEDITS_data_6_7_17_18_22_23.RData`), species traits (`MEDITS_fish_traits.RData`), and the intermediate objects `community_GSA_*.RData` and `SEM_data_all_GSA.RData`.
- [Data/Environment/](Data/Environment): environmental covariates (`environmental_data_GSA_*.RData`). Script 2 additionally expects Copernicus `.nc` downloads, `FPI_*_stack.RData` and `r_gravity_Med.RData`, which are **not present** in the repository.
- [Tools/](Tools): saved helper functions (`.RData`) used by scripts 1-3 and Figure 5. Present: the three effect-computation helpers used by `Figure_5.R`, `plot_env_na` and `select_axes_broken_stick_function`; the other helpers loaded by scripts 1-3 are missing.
- Script 2 requires a Copernicus Marine account and a Python environment (via `reticulate`); credentials must be supplied by the user in place of the placeholders in the script.
- Script 3 sources an external helper (`~/Scripts R/Divers/Variance_Inflation_Factor.R`) located outside this repository.

## Outputs
- [Outputs/](Outputs): `taxonomic_diversity_*`, `functional_diversity_*`, `other_diversity_metrics_*`, `diversity_*`, `SEM_model_ALL_GSA.RData`, `SEM_results_ALL_GSA.RData`.
- `Sub_model_selection/` (not in the repository; create it before running scripts 4a–4f): `dredged_*.RData`, only if the `save()` calls in those scripts are uncommented.
- [Outputs/Moran_I_residuals.xlsx](Outputs/Moran_I_residuals.xlsx): residual Moran's I summary (script 5).
- [Outputs/Sensitivity_analyses/](Outputs/Sensitivity_analyses): cross-validation results (`cv_results.rds`, `cv_*.csv`).
- [Figures/](Figures): PNG figures (and `Figure_S3` PDF).

Not every analysis can be rerun from the repository alone: the inputs listed above as missing, are needed for a full rerun. Figures 3 and 5 can use the saved SEM outputs; most other figures need only `Data/SEM_data_all_GSA.RData`.

## Installation and Environment
- **R version**: developed with R 4.5.x
- **Packages**: listed in [0-Load libraries.R](0-Load%20libraries.R), which installs missing CRAN packages and attaches the packages required by the scripts.

```r
source("0-Load libraries.R")
```

## Basic Usage
Each analysis script sources the loader and clears the workspace on start-up, so scripts are independent and read their inputs from disk:

```r
source("./0-Load libraries.R")
rm(list = ls())
```

Run in order from the project root:

```r
source("1-Diversity_metrics.r")
source("2-Environmental_anthropogenic_covariatesR")
source("3-Build_dataset.R")
source("4a-Build_submodel_taxo_richness.R")   # several hours
source("4b-Build_submodel_taxo_evenness.R")   # several hours
source("4c-Build_submodel_log_cwm_length.R")   # several hours
source("4d-Build_submodel_fric.R")   # several hours
source("4e-Build_submodel_feve.R")   # several hours
source("4f-Build_submodel_biomass.R")   # several hours
source("5-Check_spatial_autocor.R")
source("6-Build_SEM_model_ALL_GSA_clean.R")
source("7-Sensitivity_analyses.R")
source("8-Exploratory_analyses.R")
```

Figure scripts (`Figure_*.R`) are run individually once the required inputs exist.

## Repository Structure
- [0-Load libraries.R](0-Load%20libraries.R): dependency installation and loading.
- `1-` to `3-`: diversity metrics, covariates and SEM dataset.
- `4a-` to `4f-`: sub-model selection per response variable.
- `5-` to `8-`: residual checks, SEM, sensitivity analyses and summary statistics.
- `Figure_*.R`: main and supplementary figures.
- [Data/](Data), [Tools/](Tools), [Outputs/](Outputs), [Figures/](Figures): inputs, helper functions, saved outputs and figures.

## Computational Notes
- Scripts `4a`–`4f` call `MuMIn::dredge()` on a large candidate set while re-estimating spatial correlation parameters for each candidate; each script can take **several hours** (`4b` runs in parallel on all but one core). They are not needed to reproduce the SEM (script 6 uses hard-coded final formulas) or the figures.
- `mFD::alpha.fd.multidim()` in script 1 is computationally intensive.
- Script 7 refits models repeatedly for cross-validation; its runtime has not been measured.
- The `save()` calls in scripts 4a–4f and the `dredged_*.RData` loads in script 5 are commented out; uncomment them to pass dredging results between scripts.
