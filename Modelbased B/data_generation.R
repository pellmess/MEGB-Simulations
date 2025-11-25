#############################################################################
# Sample Generation for the MEM Boosting Simulation I                                                            
#############################################################################


#==============================================================================
# Preparation
#==============================================================================
rm(list=ls())
library(mvtnorm)
library(tidyverse)
set.seed(1)
#-------------------------------------------------
#- set simulation parameters 
#-------------------------------------------------

nPredictors <- 9
mu <- rep(0, nPredictors)
SIG <- diag( 0.5, nPredictors )
SIG[upper.tri(SIG)] <- 0.4
SIG <- SIG + t(SIG) # predictor covariance matrix
ICC <- c( .05, .25, .5 )
Sigma2_eps <- 1
n_domains <- c(100, 250) # number of known areas
ni_domain <- c(10, 30) # mean number of observations per area
N_i <- 300

study_design <- expand.grid(
  nPredictors = nPredictors,
  ICC         = ICC,
  Sigma2_eps  = Sigma2_eps,
  n_domains         = n_domains,
  ni_domain    = ni_domain,
  N_i = N_i,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)


n_mc <- 100 # number of monte_carlo sim-runs

#- path for saving
DIR <- " "

#- for saving sample statistics
sampleStats <- NULL


#==============================================================================
# generate samples
#==============================================================================

# Iterate over study_design
for(i in 1:nrow(study_design)) {
  
  # draw parameters from study_design
  tmp_nPred     <- study_design[i, "nPredictors"]
  tmp_ICC       <- study_design[i, "ICC"]
  tmp_sigma2eps <- study_design[i, "Sigma2_eps"]
  tmp_nDom      <- study_design[i, "n_domains"]
  tmp_ni        <- study_design[i, "ni_domain"]
  tmp_Ni        <- study_design[i, "N_i"]
  
  # simulate  amount n_mc many datasets
  lapply(seq_len(n_mc), function(simrun){
    TotalObs <-  tmp_nDom * tmp_Ni
    idu <- 1:TotalObs # identifier for each observation (just for sanity checks)
    idd <- rep(1:tmp_nDom, each = tmp_Ni) # area indicator

    #- predictors
    X <- data.frame( mvtnorm::rmvnorm( n = TotalObs, mean = mu, sigma = SIG ) )
    
    #- mean function
    g <- with( X, 2*X1 + X2^2 + 4*(X3 > 0) + 2*log( abs(X1) ) * X3 )
    
    #- errors
    eps <- rnorm( n = TotalObs, 0, sqrt( tmp_sigma2eps ) )
    
    #- random effects
    Sigma2_b <- ( tmp_ICC * tmp_sigma2eps ) / ( 1 - tmp_ICC )
    b <- rnorm( n = tmp_nDom, 0, sqrt( Sigma2_b ) )
    
    #- outcome
    pop_data <- data.frame( idu = idu, idd = idd,   
                       X, g = g, eps = eps, b = rep(b, each = tmp_Ni) )
    pop_data$Y <- with( pop_data, g + b + eps )
    pop_data$Y <- ifelse(pop_data$Y < 0, 0.01, pop_data$Y)

    ni <- round( runif(tmp_nDom, min = tmp_ni - 5, max = tmp_ni + 5) )
    
    # draw per area a sample of size n_i
    smp_data <- pop_data %>%
      group_by(idd) %>%
      group_modify(~ dplyr::slice_sample(.x, n = ni[.y$idd])) %>%
      ungroup()
    
    result <- list(smp_data = smp_data, pop_data = pop_data, setting = study_design[i,], simrun = simrun)
    
    stub <- sprintf(
      "n_domains=%d_ni_domain=%d_N_i=%d_ICC=%.2f_Sigma2_eps=%g_nPredictors=%d_rep=%03d",
      tmp_nDom, tmp_ni, tmp_Ni, tmp_ICC, tmp_sigma2eps, tmp_nPred, simrun
    )
    
    saveRDS(result, file = file.path("data/", paste0("run_", stub, ".rds")), compress = "xz")
  })
  
  
}

