# ==============================================================================
# Compute taxonomic, functional and trait-based diversity per haul
#
# Note: Functional diversity (alpha.fd.multidim) is computationally intensive
# ==============================================================================

# ==============================================================================
# 1. Load libraries and prepare the environment
# ==============================================================================

# Load libraries
source("./0-Load libraries.R")

rm(list = ls())

# ==============================================================================
# 2. Load data and define the community
# ==============================================================================

# Load data
data <- get(load("Data/MEDITS_data_6_7_17_18_22_23.RData")) ; names(data)

# Define the area of interest
GSA <- c(6, 7, 17, 18, 22, 23)

# Clean the data
community <- data %>% 
  dplyr::select(             
    "area" = "gsa", "haul_id" = "haul_id", "haul_dur" = "haul_dur", "gear" = "gear",                  
    "taxon" = "taxon", "wgt_cpua" = "biomass_density", "num_cpua" = "abundance_density",           
    "depth" = "depth", "year" = "year", "month" = "month", "latitude" = "latitude", "longitude" = "longitude") %>% 
  filter(area %in% GSA)

# Check data
length(unique(community$haul_id))      
length(unique(community$taxon)) 

# Save filtered community data
# save(community, file = paste0("Data/community_GSA_", paste0(GSA, collapse = "_"), ".RData"))

# ==============================================================================
# 3. Taxonomic diversity
# ==============================================================================

# Taxonomic diversity
community <- get(load(paste0("Data/community_GSA_", paste0(GSA, collapse = "_"), ".RData")))

### Compute richness, evenness (Pielou's), and diversity (Shannon-Wiener) ----
taxonomic_diversity <- community %>%
  group_by(haul_id) %>%
  summarize(tric = vegan::specnumber(taxon),                      # richness
            tdiv = vegan::diversity(wgt_cpua, index = "shannon"), # diversity
            teve = vegan::diversity(wgt_cpua, index = "shannon")/ log(specnumber(taxon))) # evenness

# Identify hauls with taxonomic richness = 1 and remove them from the taxonomic diversity dataset (since evenness is undefined for richness = 1)
hauls_richness_1 <- taxonomic_diversity %>%
  filter(tric == 1) %>%
  pull(haul_id)

# Remove hauls with taxonomic richness = 1 from the taxonomic diversity dataset and community dataset (since they will be excluded from the functional diversity analysis as well)
taxonomic_diversity <- taxonomic_diversity %>%
  filter(!haul_id %in% hauls_richness_1)

community <- community %>%
  filter(!haul_id %in% hauls_richness_1)

# Save
# save(taxonomic_diversity, file = paste0("Outputs/taxonomic_diversity_GSA_", paste0(GSA, collapse = "_"), ".RData"))
  
# ==============================================================================
# 4. Functional diversity
# ==============================================================================

# Functional diversity
traits <- get(load("Data/MEDITS_fish_traits.RData"))

# Merge traits with survey species
community_traits <- data.frame(taxon = unique(community$taxon)) %>%
  left_join(traits, by = c("taxon" = "Medits_code"))

# Check how many traits are NAs
missing_traits <- community_traits[rowSums(is.na(community_traits)) > 0, ] 
colSums(is.na(missing_traits)) 

# Traits to exclude
exclude_t <- c("Length_infinity", "Natural_mortality", "Weight_infinity", "Length_maturity",
               "Parental_care", "Water_column_position", "Aspect_ratio", "Age_maturity")

# Number of hauls with fewer than 3
n_hauls_fewer_than_3 <- community %>%
  group_by(haul_id) %>%
  summarize(n_species = n_distinct(taxon)) %>%
  filter(n_species < 3) %>%
  nrow()

# Species matrix
sm <- community %>% 
  dplyr::select(haul_id, taxon = taxon, count = wgt_cpua) %>%
  spread(taxon, count) %>%                        
  filter(rowSums(!is.na(.)) > 2) %>%  # Prior to analyses, we excluded hauls with less than three species because these might lead to uninformative trait spaces (Mammola et al., 2024)
  mutate(across(everything(), ~ replace_na(., 0))) %>% 
  column_to_rownames(var = "haul_id")

# Traits matrix
tm <- community_traits %>%
  dplyr::select(!(exclude_t)) %>% # exlude unwanted traits
  column_to_rownames(var = "taxon") %>% 
  mutate(across(where(is.character), as.factor)) # set to factors to compute gawdis

# Check that tm's row names match sm's column names (must be TRUE)
all(sort(rownames(tm)) == sort(colnames(sm))) 

# List of traits included in the analysis
colnames(tm)

# Check collinearity among all traits and Variance Inflation Factor (VIF)
df_cor <- tm[, which(sapply(tm, is.numeric))]
psych::pairs.panels(df_cor, stars = TRUE, method = "spearman", hist.col = "lightblue", density = TRUE, ellipses = TRUE)

# Build trait type reference table
fish_traits_cat <- data.frame(
  trait_name = colnames(tm),
  trait_type = sapply(tm, function(col) {
    if (is.numeric(col))        "Q" else   # quantitative
    if (is.ordered(col))        "O" else   # ordinal
                                "N"        # nominal
  }), row.names = NULL)

# Compute functional distance matrix using Gower's distance
sp_dist_fish <- mFD::funct.dist(sp_tr = tm, tr_cat = fish_traits_cat, metric = "gower", scale_euclid = "scale_center", ordinal_var = "classic", weight_type = "equal", stop_if_NA = TRUE)

# Check trait contribution to the dissimilarity matrix
attr(sp_dist_fish, "correls")

# Compute multidimensional functional spaces and assess their quality
fspaces_quality_fish <- mFD::quality.fspaces(sp_dist = sp_dist_fish, maxdim_pcoa = 10, deviation_weighting = "absolute", fdist_scaling = FALSE, fdendro = "average")

# Illustrating the quality of the selected functional spaces
mFD::quality.fspaces.plot(
  fspaces_quality = fspaces_quality_fish, quality_metric = "mad",
  fspaces_plot = c("tree_average", "pcoa_2d", "pcoa_3d", "pcoa_4d", "pcoa_5d", "pcoa_6d"),
  name_file = NULL, range_dist = NULL, range_dev = NULL, range_qdev = NULL,
  gradient_deviation = c(neg = "darkblue", nul = "grey80", pos = "darkred"),
  gradient_deviation_quality = c(low = "yellow", high = "red"), x_lab = "Trait-based distance")

# Test correlation between functional axes and traits
sp_faxes_coord_fish <- fspaces_quality_fish$details_fspaces$sp_pc_coord

fish_tr_faxes <- mFD::traits.faxes.cor(sp_tr = tm, sp_faxes_coord = sp_faxes_coord_fish[ , c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")], plot = FALSE)

# Print traits with significant effect:
fish_tr_faxes[which(fish_tr_faxes$p.value < 0.05), ]

# Look number of dimensions to include in the functional space and compute variance explained by the selected axes 
load("Tools/select_axes_broken_stick_function.RData")
fspace_axes <- select_axes_broken_stick(pc_eigenvalues = fspaces_quality_fish$details_fspaces$pc_eigenvalues, mad_vals = fspaces_quality_fish$quality_fspaces$mad, keep = c("lowest_mad"), plot = TRUE, verbose = TRUE)

# Functional alpha diversity indices in a multidimensional space
# number of dimensions of the functional space to include in the analysis
n_dimensions <- 4 # or fspace_axes$n_axes

s <- sm %>% filter(rowSums(. > 0, na.rm = TRUE) > n_dimensions) %>% # convex hull requirement
  as.matrix() 

# compute relative abundance of each species 
s_rel <- funrar::make_relative(as.matrix(s))

# Compute functional diversity indices in a multidimensional space
# NOTE: mFD::alpha.fd.multidim() is computationally intensive and may take a long time to run. Do not rerun it
# use the saved Outputs/functional_diversity_GSA_*.RData
alpha_fd_indices_fish <- mFD::alpha.fd.multidim(
  sp_faxes_coord = sp_faxes_coord_fish[ , paste0("PC", 1:n_dimensions)],
  asb_sp_w = s, ind_vect = c("feve", "fric", "fdis", "fmpd", "fnnd", "fdiv", "fori", "fspe", "fide"),
  scaling = TRUE, check_input = TRUE, details_returned = TRUE)

# Functional redundancy index
functio_diversity <-alpha.fd.hill(asb_sp_w = s_rel, sp_dist = sp_dist_fish, q = c(1,2), tau = "mean")
summary(functio_diversity$asb_FD_Hill)
functio_q1 <- functio_diversity$asb_FD_Hill[,1]
taxo_q1 <- exp(vegan::diversity(s_rel))
functio_redund_q1 <- 1-(functio_q1/taxo_q1)

# Extract the indices from the list
functional_diversity <- alpha_fd_indices_fish$functional_diversity_indices %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "haul_id") %>%
  dplyr::select(-sp_richn, -paste0("fide_PC", 1:n_dimensions)) %>%
  left_join(tibble::enframe(functio_redund_q1, name = "haul_id", value = "fredund"), by = "haul_id")

# Save
# save(functional_diversity, file = paste0("Outputs/functional_diversity_GSA_", paste0(GSA, collapse = "_"), ".RData"))

# ==============================================================================
# 5. Community-weighted and trophic metrics
# ==============================================================================

### Community weighted metrics and trophic metrics ----
# Community weighted mean and variance of selected continuous traits
trait_list <- c("Length_max") 

weighted_metrics <- community %>%
  left_join(community_traits, by = "taxon") %>%
  group_by(haul_id) %>%
  summarize(across(all_of(trait_list), list(
        cwm = ~weighted.mean(.x, wgt_cpua, na.rm = TRUE),
        cwv = ~modi::weighted.var(.x, wgt_cpua, na.rm = TRUE)), .names = "{.col}_{.fn}"))

# Benthivory vs planktivory
feeding <- community %>%
  left_join(community_traits, by = "taxon") %>%
  group_by(haul_id, Trophic_guild) %>%
  summarize(bm = sum(wgt_cpua, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Trophic_guild, values_from = bm, values_fill = 0) %>%
  mutate(
    bent_vs_plank = if_else(
      planktivorous > 0,
      benthivorous / planktivorous,
      NA_real_  # hauls with no planktivorous biomass → undefined ratio
    )
  ) %>%
  dplyr::select(haul_id, bent_vs_plank)

# Join
other_diversity_metrics <- left_join(weighted_metrics, feeding, by = "haul_id")

# Save
# save(other_diversity_metrics, file = paste0("Outputs/other_diversity_metrics_GSA_", paste0(GSA, collapse = "_"), ".RData"))

# Join all metrics
# Load data
load(paste0("Outputs/taxonomic_diversity_GSA_", paste0(GSA, collapse = "_"), ".RData"))
load(paste0("Outputs/functional_diversity_GSA_", paste0(GSA, collapse = "_"), ".RData"))
load(paste0("Outputs/other_diversity_metrics_GSA_", paste0(GSA, collapse = "_"), ".RData"))

# keep only hauls present in all datasets / just to be sure
common_hauls <- Reduce(intersect, list(taxonomic_diversity$haul_id, functional_diversity$haul_id, other_diversity_metrics$haul_id))

# Join
diversity <- taxonomic_diversity %>%
  filter(haul_id %in% common_hauls) %>%
  left_join(functional_diversity %>% filter(haul_id %in% common_hauls), by = "haul_id") %>%
  left_join(other_diversity_metrics %>% filter(haul_id %in% common_hauls), by = "haul_id")

# Save
# save(diversity, file = paste0("Outputs/diversity_GSA_", paste0(GSA, collapse = "_"), ".RData"))
