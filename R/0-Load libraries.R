# ==============================================================================

# ==============================================================================
# 1. Package lists
# ==============================================================================

required_packages <- c(
  "reticulate", "copernicusR", "usethis", "lubridate", "terra",
  "MASS", "mgcv", "nlme", "gstat", "ape", "car", "sp", "geoR", "spdep",
  "MuMIn", "pbapply", "piecewiseSEM", "marginaleffects", "ggeffects",
  "vegan", "mFD", "BAT", "modi", "funrar", "janitor", "GGally",
  "rstatix", "multcompView", "writexl", "ggmap", "ggcorrplot", "scico",
  "DiagrammeR", "DiagrammeRsvg", "rsvg",
  "rnaturalearth", "sf",
  "tibble", "tidyr", "dplyr", "ggplot2", "patchwork"
)

# Installed but only used through `pkg::fun()`
namespace_only_packages <- c(
  "usdm", "fastDummies", "stringr", "forcats", "viridis", "psych", "remotes"
)

# ==============================================================================
# 2. Install missing packages
# ==============================================================================

installed_packages <- rownames(installed.packages())

to_install <- setdiff(
  c(required_packages, namespace_only_packages),
  installed_packages
)

if (length(to_install) > 0) {
  message("Installing missing CRAN packages: ", paste(to_install, collapse = ", "))
  tryCatch(
    install.packages(to_install, dependencies = TRUE),
    error = function(e) {
      stop("Package installation failed: ", conditionMessage(e), call. = FALSE)
    }
  )
}

# ==============================================================================
# 3. Load packages
# ==============================================================================

missing_after_install <- setdiff(
  c(required_packages, namespace_only_packages),
  rownames(installed.packages())
)

if (length(missing_after_install) > 0) {
  stop(
    "Required packages could not be installed: ",
    paste(missing_after_install, collapse = ", "),
    call. = FALSE
  )
}

for (pkg in required_packages) {
  tryCatch(
    suppressPackageStartupMessages(
      library(pkg, character.only = TRUE)
    ),
    error = function(e) {
      stop("Package '", pkg, "' could not be loaded: ", conditionMessage(e),
           call. = FALSE)
    }
  )
}

rm(installed_packages, to_install, missing_after_install, pkg)
