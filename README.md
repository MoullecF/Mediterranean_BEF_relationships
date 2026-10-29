# Biodiversity-ecosystem functioning relationships in Mediterranean fish communities
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [![Data: CC BY 4.0](https://img.shields.io/badge/Data%20License-CC%20BY%204.0-brightgreen.svg)](https://creativecommons.org/licenses/by/4.0/)

This repository contains the R code used to investigate biodiversity-ecosystem functioning (BEF) relationships in fish communities across four northern Mediterranean sub-regions: Northern Spain, the Gulf of Lion, the Adriatic Sea and the Aegean Sea. It analyses 23 years (1999-2021) of Mediterranean International Trawl Survey (MEDITS) data (10,034 hauls, 331 species) with piecewise structural equation modelling, relating environmental conditions and anthropogenic pressures to taxonomic diversity, functional diversity, community-weighted maximum body length and fish biomass. Fish biomass is treated as a biomass-related proxy for ecosystem functioning.

All figures can be regenerated directly from the provided outputs, without rerunning scripts 1-8.

## Purpose ![target](https://img.shields.io/badge/-purpose-37474f?style=flat-square)
- Compute taxonomic, functional and trait-based diversity metrics for each trawl haul
- Assemble environmental and anthropogenic covariates at haul positions
- Select and fit the sub-models of a piecewise SEM with spatial correlation structures
- Summarise direct, indirect and total effects and run sensitivity analyses

## Methods and Workflow ![workflow](https://img.shields.io/badge/-workflow-00695c?style=flat-square)
1. **Diversity metrics**: taxonomic, functional (mFD) and community-weighted trait metrics per haul ([1-Diversity_metrics.r](R/1-Diversity_metrics.r)).
2. **Covariates**: Copernicus Marine environmental variables, fishing pressure and gravity extracted at haul positions ([2-Environmental_anthropogenic_covariates.R](R/2-Environmental_anthropogenic_covariates.R)).
3. **SEM dataset**: merging, filtering and collinearity screening ([3-Build_dataset.R](R/3-Build_dataset.R)).
4. **Sub-model selection** (spatial correlation structure, then `MuMIn::dredge()`), one script per response:
   [4a](R/4a-Build_submodel_taxo_richness.R), [4b](R/4b-Build_submodel_taxo_evenness.R), [4c](R/4c-Build_submodel_log_cwm_length.R), [4d](R/4d-Build_submodel_fric.R), [4e](R/4e-Build_submodel_feve.R), [4f](R/4f-Build_submodel_biomass.R).
5. **Residual spatial autocorrelation**: Moran's I of the final sub-models ([5-Check_spatial_autocor.R](R/5-Check_spatial_autocor.R)).
6. **Piecewise SEM**: final model fitted to all GSAs ([6-Build_SEM_model.R](R/6-Build_SEM_model.R)).
7. **Sensitivity analyses**: alternative biomass specifications and spatial block cross-validation ([7-Sensitivity_analyses.R](R/7-Sensitivity_analyses.R)).
8. **Spatial exploratory analyses**: correlations and among-GSA comparisons ([8-Exploratory_analyses.R](R/8-Exploratory_analyses.R)).
9. **Figures**: [Figure_1.R](R/Figure_1.R) to [Figure_6.R](R/Figure_6.R) (main text) and [Figure_S1.R](R/Figure_S1.R) to [Figure_S5.R](R/Figure_S5.R) (supporting information).

## Inputs and Data ![data](https://img.shields.io/badge/-data-283593?style=flat-square)
- [Data/](Data): MEDITS survey data (`MEDITS_data_6_7_17_18_22_23.RData`), species traits (`MEDITS_fish_traits.RData`), and the intermediate objects `community_GSA_*.RData` and `SEM_data_all_GSA.RData`. For the MEDITS data cleaning and consistency checks, see [COISPA/RoME](https://github.com/COISPA/RoME), an R package for checking MEDITS survey data.
- [Data/Environment/](Data/Environment): environmental covariates (`environmental_data_GSA_*.RData`). Script 2 additionally expects Copernicus `.nc` downloads, `FPI_*_stack.RData` and `r_gravity_Med.RData`; these raw files are too large for the repository and are available on request. Copernicus data can also be downloaded from the public catalogue (free user account required).
- [Tools/](Tools): saved helper functions (`.RData`) used by scripts 1-3 and Figure 5.

## Outputs ![outputs](https://img.shields.io/badge/-outputs-5d4037?style=flat-square)
- [Outputs/](Outputs): `taxonomic_diversity_*`, `functional_diversity_*`, `other_diversity_metrics_*`, `diversity_*`, `SEM_model_ALL_GSA.RData`, `SEM_results_ALL_GSA.RData`.
- `Sub_model_selection/`: `dredged_*.RData`.
- [Outputs/Moran_I_residuals.xlsx](Outputs/Moran_I_residuals.xlsx): residual Moran's I summary.
- [Outputs/Sensitivity_analyses/](Outputs/Sensitivity_analyses): cross-validation results (`cv_results.rds`, `cv_*.csv`).
- [Figures/](Figures): PNG figures.

## Installation and Environment ![setup](https://img.shields.io/badge/-setup-546e7a?style=flat-square)
- **R version**: developed with R 4.5.x
- **Core packages**: listed in [0-Load libraries.R](R/0-Load%20libraries.R), which installs missing CRAN packages and attaches the packages required by the scripts.

```r
source(here::here("R/0-Load libraries.R"))
```

## Basic Usage ![run](https://img.shields.io/badge/-run-1b5e20?style=flat-square)
Run the scripts in order. All scripts are in `R/` and use `here::here()` for paths, so they can be run from any working directory within the project. Each script sources the loader (`R/0-Load libraries.R`) and clears the workspace on start-up, so scripts are independent.

1. Compute diversity metrics from the MEDITS data:
   ```r
   source(here::here("R/1-Diversity_metrics.r"))
   ```
2. Extract environmental and anthropogenic covariates (requires the raw files listed in `Data/Environment/`):
   ```r
   source(here::here("R/2-Environmental_anthropogenic_covariates.R"))
   ```
3. Assemble the SEM dataset:
   ```r
   source(here::here("R/3-Build_dataset.R"))
   ```
4. Build the sub-models: spatial correlation structure, then predictor selection with `MuMIn::dredge()` (several hours (>20h) per script):
   ```r
   source(here::here("R/4a-Build_submodel_taxo_richness.R"))
   source(here::here("R/4b-Build_submodel_taxo_evenness.R"))
   source(here::here("R/4c-Build_submodel_log_cwm_length.R"))
   source(here::here("R/4d-Build_submodel_fric.R"))
   source(here::here("R/4e-Build_submodel_feve.R"))
   source(here::here("R/4f-Build_submodel_biomass.R"))
   ```
5. Check residual spatial autocorrelation:
   ```r
   source(here::here("R/5-Check_spatial_autocor.R"))
   ```
6. Fit the SEM:
   ```r
   source(here::here("R/6-Build_SEM_model.R"))
   ```
7. Run sensitivity and spatial exploratory analyses:
   ```r
   source(here::here("R/7-Sensitivity_analyses.R"))
   source(here::here("R/8-Exploratory_analyses.R"))
   ```
8. Generate figures as needed (`R/Figure_*.R`). All required outputs are already provided in the repository, so figures can be reproduced without rerunning steps 1 to 7.

## Repository Structure ![folders](https://img.shields.io/badge/-folders-3949ab?style=flat-square)
- [0-Load libraries.R](R/0-Load%20libraries.R): dependency installation and loading.
- `1-` to `3-`: diversity metrics, covariates and SEM dataset.
- `4a-` to `4f-`: sub-model selection per response variable.
- `5-` to `8-`: residual checks, SEM, sensitivity analyses and summary statistics.
- `Figure_*.R`: main and supplementary figures.
- [:file\_folder: **Tools**](Tools/): Helper functions
- [:file\_folder: **Data**](Data/): Input data
- [:file\_folder: **Outputs**](Outputs/): Saved outputs
- [:file\_folder: **Figures**](Figures/): Saved figures

## Computational Notes
- Scripts `4a`–`4f` call `MuMIn::dredge()` on a large candidate set while re-estimating spatial correlation parameters for each candidate; each script can take **several hours**. They are not needed to reproduce the SEM or the figures.
- `mFD::alpha.fd.multidim()` in script 1 is computationally intensive.
- Script 7 refits models repeatedly for cross-validation (can be long to run)

## Citation ![cite](https://img.shields.io/badge/-cite-4e342e?style=flat-square)
If you use this code or data, please cite:
- **Code**: This GitHub repository (XXXXX/Mediterranean_BEF_relationships, main branch).
- **Associated publication**: Manuscript under review.

## License ![license](https://img.shields.io/badge/-license-263238?style=flat-square)
This project is released under the GNU General Public License v3.0. See [LICENSE](LICENSE) for details.