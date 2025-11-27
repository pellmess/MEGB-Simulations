library(pacman)
p_load(saeSim)



source("MEM_Boosting_VersionMERT.R")
source("MEM_Boosting_VersionREEM.R")
source("gtb.R")

# devtools::install_github("krennpa/SAEforest", force = TRUE)
library(SAEforest)
p_load(tidyverse)
p_load(checkmate)
p_load(ranger)
p_load(emdi)
p_load(lme4)
p_load(tictoc)
p_load(mboost)
p_load(randomForest)
p_load(fastDummies)
p_load(xgboost)
# Load auxiliary functions for benchmark estimators
source("auxiliar/BHF_estimation.R") # Functions to obtain BHF mean estimates
source("auxiliar/ebp.R")

# Install MEGB package from local source
# remotes::install_local("MEGB_0.0.0.9000.tar.gz", dependencies = TRUE, type = "source")
library(MEGB)

set.seed(1)




# Choose which estimators you want to use
estimators <- c(
  "bhf",
  "BoostMERT",   # MERTBoosting, Salditt et al.
  "BoostREEM",   # REEMBoosting, Salditt et al.
  "ebp",
  "MBOOST",      # mboost, Salditt et al.
  "MBOOST_T",    # tree-based mboost, Salditt et al.
  "MEGB", 
  "MERF",        # MERF, Krennmair
  "RF",
  "GB",          # Gradient Boosting implementation by Salditt
  "xgboost"
)

run_mod <- function(mod = "mod1", which_estimators = "ebp") {
  load(paste0(mod, "/model_specf_objects.RData"))
  
  ######
  r <<- 1:5
  
  ######
  # Preparations for different estimators
  vars <- input_elements[[1]]$sample %>% names()
  x.names <- vars[grep("x", vars)]
  
  # Define target and predictor variable names
  y.names <- "y"
  
  formel <- as.formula(paste(y.names, paste(x.names, collapse = " + "), sep = " ~ "))
  
  group.id <- "idD"
  
  # Source estimator scripts; each automatically saves results as .csv and .RData
  
  ###### Loop over estimator names
  for (est_name in which_estimators) {
    print(paste0(mod, " ", est_name, " started"))
    source(paste0("estimators/", est_name, ".R"), local = TRUE)
    print(paste0(mod, " ", est_name, " successfully completed"))
  }
}

# mods <- paste0("mod", c(1:4))
mods <- paste0("mod", 1)

for (mod_name in mods) {
  run_mod(mod_name, which_estimators = estimators)
}
