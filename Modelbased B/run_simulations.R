################################################################################
# MEM Boosting Simulation I                                                      
################################################################################

set.seed(1)
#==============================================================================
# Preparation
#==============================================================================

#remotes::install_local("MEGB_0.0.0.9000.tar.gz", dependencies = TRUE, type = "source")
library(MEGB)
1
# devtools::install_github("krennpa/SAEforest", force = TRUE)

rm(list = ls())


#- load packages and functions
library(lme4)
library(randomForest)
library(mboost)
library(tictoc)
library(tidyverse)
source("MEM_Boosting_VersionMERT.R")
source("MEM_Boosting_VersionREEM.R")
source("gtb.R")

#-------------------------------------------------
#- set simulation parameters 
#-------------------------------------------------

#- specify model
model <- Y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9
random <- ~ 1 | idd
modelLMM <- Y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + (1 | idd)
modelLMMcor  <- Y ~ X1 + I(X2 ^ 2) + I(X3 > 0) + I(log(abs(X1))):X3 + (1 |
                                                                         idd)

#- set parameters for MEM boosting
conv_memboost = 0.001
maxIter_memboost = 100
verbose_memboost = F

minIter_memboost = 0
cp = 0.01

n.trees = 100
interaction.depth = 6
shrinkage = .1

bag.fraction = 0.5
minsplit = 20



#- for getting the samples 
DIR <- "data"
allFiles <- dir(DIR, full.names = TRUE, recursive = TRUE)
# 

#==============================================================================
# simulation start
#==============================================================================



for(i in 1:length(allFiles)){
  
  path <- paste0("results/")
  dir.create(path, showWarnings = F)
  # gradient process ----------------------------------------------------------
  if(file.exists(paste0(path, i, ".txt"))){
    warning(paste("Run no.", i, "has already been completed, skipping"))
    next
  } else{
    # Check if this simulation run has already finished
    write.table(x = i, paste0(path, i, ".txt"))
  
    sim_obj <- readRDS(allFiles[i])
    smp_data <- sim_obj$smp_data[, c("Y", paste0("X", 1:9), "idd")] |> as.data.frame()
    smp_data$idd <- as.factor(smp_data$idd)
    pop_data <- sim_obj$pop_data[, c(paste0("X", 1:9), "idd")] |> as.data.frame()
    pop_data$idd <- as.factor(pop_data$idd)
    #-------------------------------------------------------------------------------
    #- Standard boosting
    #-------------------------------------------------------------------------------
    # 
    tictoc::tic("Standard Boosting L2")
    s_Boost_L2 <- Sys.time()
    fitBoost_L2 <- gtb(model, data = smp_data,
                        n.trees = n.trees,
                        loss = "L2",
                        interaction.depth = interaction.depth,
                        shrinkage = shrinkage,
                        bag.fraction = bag.fraction,
                        minsplit = minsplit )
    e_Boost_L2 <- Sys.time()
    tictoc::toc(log = T)
    time_Boost_L2 <- difftime(e_Boost_L2, s_Boost_L2, units = "secs")[[1]]

    pop_data_pred <- predict.gtb( fitBoost_L2, pop_data, n.trees = n.trees )
    area_est_L2 <- cbind.data.frame(pop_data_pred, idd = pop_data$idd) |> 
      group_by(idd) |> 
      summarise(mw = mean(pop_data_pred))



  #-------------------------------------------------------------------------------
  #- MEM Boosting 
  #-------------------------------------------------------------------------------
  
  tictoc::tic("MERT Boosting L2")
  s_BoostMERT_L2 <- Sys.time()
  fitBoostMERT_L2 <- mem_boost_mert( model, data = smp_data,
                                     random = random,
                                     loss = "L2",
                                     n.trees = n.trees,
                                     interaction.depth = interaction.depth,
                                     shrinkage = shrinkage,
                                     bag.fraction = bag.fraction,
                                     minsplit = minsplit,
                                     maxIter_memboost = maxIter_memboost,
                                     minIter_memboost = minIter_memboost,
                                     verbose_memboost = verbose_memboost,
                                     cp = cp )
  e_BoostMERT_L2 <- Sys.time()
  time_BoostMERT_L2 <- difftime(e_BoostMERT_L2, s_BoostMERT_L2, units = "secs")[[1]]
  tictoc::toc()

  tictoc::tic("REEM Boosting L2")
  s_BoostREEM_L2 <- Sys.time()
  fitBoostREEM_L2 <- mem_boost_reem( model, data = smp_data,
                                     random = random,
                                     loss = "L2",
                                     n.trees = n.trees,
                                     interaction.depth = interaction.depth,
                                     shrinkage = shrinkage,
                                     bag.fraction = bag.fraction,
                                     minsplit = minsplit,
                                     maxIter_memboost = maxIter_memboost,
                                     minIter_memboost = minIter_memboost,
                                     verbose_memboost = verbose_memboost,
                                     cp = cp )
  e_BoostREEM_L2 <- Sys.time()
  time_BoostREEM_L2 <- difftime(e_BoostREEM_L2, s_BoostREEM_L2, units = "secs")[[1]]
  tictoc::toc()


  fitList_MEMBoost <- list(   BoostMERT_L2 = fitBoostMERT_L2,
                              BoostREEM_L2 = fitBoostREEM_L2 )

  #- Calculate predicted values
  Yhat_MEMBoost <- lapply(fitList_MEMBoost, function(X){

    # 1) random effects as named vector (names = idd of training data)
    ids_re <- as.character(unique(smp_data$idd))
    re_vec <- setNames(as.numeric(X$raneffs), ids_re)

    # 2) Align with the rows of pop_data (no sorting, no merging)
    re_row <- re_vec[ as.character(pop_data$idd) ]

    # Safety: incase unknown levels exist -> report error
    if (anyNA(re_row)) {
      stop("MEM-Boost: unknown idd found in pop_data (not present in training).
         Ensure all pop_data$idd values appear in smp_data$idd.")
    }
    

    # 3) fixed effects prediction + adding random effects
    fhat <- predict.gtb(X$boosting_ensemble, pop_data)  # same row order as pop_data
    pop_data_pred <- fhat + re_row

    # 4) aggregate to area-level + stable order and Row names
    out <- cbind.data.frame(pop_data_pred, idd = pop_data$idd) |>
      dplyr::mutate(idd = as.integer(as.character(idd))) |>
      dplyr::group_by(idd) |>
      dplyr::summarise(mw = mean(pop_data_pred), .groups = "drop") |>
      dplyr::arrange(idd)

    rownames(out) <- out$idd
    out
  })

  area_est_BoostMERT_L2 <- Yhat_MEMBoost[[1]]
  area_est_BoostREEM_L2 <- Yhat_MEMBoost[[2]]

  #-------------------------------------------------------------------------------
  #- Model-based boosting using simple linear models as base-learners
  #-------------------------------------------------------------------------------

  dfs_mboost <- list( smp_data = smp_data,
                      pop_data = pop_data)
  dfs_mboost <- lapply( dfs_mboost, function(Y){
    # add intercept variable
    Y$int <- rep(1, nrow(Y))
    # mean-center the predictors
    predX <- grepl("X[1-9]", colnames(Y))
    Y[, predX] <- apply(Y[, predX], 2, function(x) scale(x, center = T, scale = F))
    return(Y)  } )


  tictoc::tic("MBOOST_L2")
  s_MBOOST_L2 <- Sys.time()
  fitMBOOST_L2 <- mboost::gamboost( Y ~ bols(int, intercept = F) + # intercept
                                      bols(X1, intercept = F) + bols(X2, intercept = F) +
                                      bols(X3, intercept = F) +
                                      bols(X4, intercept = F) + bols(X5, intercept = F) +
                                      bols(X6, intercept = F) +
                                      bols(X7, intercept = F) + bols(X8, intercept = F) +
                                      bols(X9, intercept = F) +
                                      brandom(idd),
                                    data = dfs_mboost$smp_data)
  cvm <- cvrisk( fitMBOOST_L2 ) # find the optimal number of boosting iterations
  fitMBOOST_L2[ mstop(cvm) ] # choose the model based on optimal number of boosting iterations
  e_MBOOST_L2 <- Sys.time()
  time_MBOOST_L2 <- difftime(e_MBOOST_L2, s_MBOOST_L2, units = "secs")[[1]]
  tictoc::toc(log = TRUE)

  pop_data_pred <- predict(fitMBOOST_L2, dfs_mboost$pop_data)
  area_est_MBOOST_L2 <- cbind.data.frame(pop_data_pred, idd = pop_data$idd) |> 
    group_by(idd) |> 
    summarise(mw = mean(pop_data_pred))
  
  
  
  #-------------------------------------------------------------------------------
  #- Model-based boosting using a tree as base learners
  #-------------------------------------------------------------------------------
  
  # using trees as base learner and additionally a random intercept
  tictoc::tic("MBOOST_L2 trees")
  s_MBOOST_T_L2 <- Sys.time()
  fitMBOOST_T_L2 <- mboost::gamboost( Y ~ btree(X1, X2, X3, X4, X5, X6, X7, X8, X9) +
                                        brandom(idd),
                                      data = dfs_mboost$smp_data)
  cvm <- cvrisk( fitMBOOST_T_L2 )
  fitMBOOST_T_L2[ mstop(cvm) ]
  e_MBOOST_T_L2 <- Sys.time()
  time_MBOOST_T_L2 <- difftime(e_MBOOST_T_L2, s_MBOOST_T_L2, units = "secs")[[1]]
  tictoc::toc(log = TRUE)

  pop_data_pred <- predict(fitMBOOST_T_L2, newdata = dfs_mboost$pop_data)
  area_est_MBOOST_T_L2 <- cbind.data.frame(pop_data_pred, idd = pop_data$idd) |> 
    group_by(idd) |> 
    summarise(mw = mean(pop_data_pred))

  # 
  # 
  #-------------------------------------------------------------------------------
  #- Standard Random Forest 
  #-------------------------------------------------------------------------------
  tictoc::tic("Standard RF")
  s_RF <- Sys.time()
  fitRF <- randomForest::randomForest(formula = model, data = smp_data)
  e_RF <- Sys.time()
  tictoc::toc(log = T)
  time_RF <- difftime(e_RF, s_RF, units = "secs")[[1]]

  pop_data_pred <- predict(fitRF, pop_data)
  area_est_RF <- cbind.data.frame(pop_data_pred, idd = pop_data$idd) |> 
    group_by(idd) |> 
    summarise(mw = mean(pop_data_pred))

  #-------------------------------------------------------------------------------
  #- MERF 
  #-------------------------------------------------------------------------------
    library(SAEforest)
    ts.MERF = Sys.time()
    MERF <- SAEforest_model(
      Y = smp_data$Y,
      X = smp_data[,paste0("X", 1:9)],
      dName = "idd",
      smp_data = smp_data,
      pop_data = pop_data,
      seed = 123, meanOnly = T
    )
    te.MERF = Sys.time()

    time_MERF <- difftime(te.MERF, ts.MERF, units = "secs")[[1]]
    area_est_MERF <- MERF$Indicators
    names(area_est_MERF) <- c("idd", "mw")
    # sort in numeric order
    area_est_MERF$idd <- as.integer(as.character(area_est_MERF$idd))
    area_est_MERF <- area_est_MERF[order(area_est_MERF$idd), , drop = FALSE]
    row.names(area_est_MERF) <- NULL
  #-------------------------------------------------------------------------------
  #- EBP 
  #-------------------------------------------------------------------------------
   formul_ebp <-
     as.formula(paste("Y", paste(paste0("X", 1:9), collapse = " + "), sep = " ~ "))
   ts.ebp = Sys.time()
   ebp_results <- emdi::ebp(
     fixed = formul_ebp,
     pop_data = pop_data,
     smp_data =  smp_data,
     pop_domains = "idd",
     smp_domains = "idd",
     transformation = "box.cox",
   )
   te.ebp= Sys.time()
   time_ebp <- difftime(te.ebp, ts.ebp, units = "secs")[[1]]
   area_est_ebp <- ebp_results$ind %>% dplyr::select(domain = Domain, Mean)
   names(area_est_ebp) <- c("idd", "mw")


    
  
  #-------------------------------------------------------------------------------
  #- MEGB 
  #-------------------------------------------------------------------------------
  params <- list(
    eta = 0.1,              # learning rate
    nrounds = 10000,        # number of trees
    max_depth = 3,          # maximal depth of trees
    min_child_weight = 3,   # minimal child weight
    subsample = 0.5,        # proportion of the sample used to train each tree
    lambda = 1,             # L2-regularisation
    alpha = 0,              # L1-regularisation
    gamma = 0.9
  )
  tictoc::tic("MEGB")
  s_MEGB <- Sys.time()
  megb_obj <- MEGB::megb(Y = smp_data$Y,
                         X = smp_data[,paste0("X", 1:9)],
                         dom_name = "idd",
                         smp_data = smp_data,
                         pop_data = pop_data,
                         gradient_params = params,
                         seed = 1,
                         mse = F,
                         bootstrap_cores = 1,
                         B = 100)
  e_MEGB <- Sys.time()
  tictoc::toc(log = T)
  time_MEGB <- difftime(e_MEGB, s_MEGB, units = "secs")[[1]]
  area_est_MEGB <- megb_obj$Indicators
  names(area_est_MEGB) <- c("idd", "mw")
  area_est_MEGB
  
  
    # -------------------------------------------------------------------------------
    # - BHF
    # -------------------------------------------------------------------------------

    tictoc::tic("BHF")
    s_bhf <- Sys.time()
    bhf <- saeTrafo::NER_Trafo(fixed = model,
                                    pop_data = pop_data,
                                    smp_data =  smp_data,
                                    pop_domains = "idd",
                                    smp_domains = "idd",
                                    transformation = "no")
    e_bhf <- Sys.time()
    tictoc::toc(log = T)
    time_bhf <- difftime(e_bhf, s_bhf, units = "secs")[[1]]
    area_est_bhf <- bhf$ind
    names(area_est_bhf) <- c("idd", "mw")
    area_est_bhf

    results_save <- list(
      area_est_L2 = area_est_L2,
      area_est_BoostMERT_L2 = area_est_BoostMERT_L2,
      area_est_BoostREEM_L2 = area_est_BoostREEM_L2,
      area_est_MBOOST_L2 = area_est_MBOOST_L2,
      area_est_MBOOST_T_L2 = area_est_MBOOST_T_L2,
      area_est_RF = area_est_RF,
      area_est_MEGB = area_est_MEGB,
      area_est_ebp = area_est_ebp,
      area_est_MERF = area_est_MERF,
      area_est_bhf = area_est_bhf,
      time_MEGB = time_MEGB,
      time_ebp = time_ebp,
      time_MERF = time_MERF,
      time_bhf = time_bhf,
      time_Boost_L2 = time_Boost_L2,
      time_BoostMERT_L2 = time_BoostMERT_L2,
      time_BoostREEM_L2 = time_BoostREEM_L2,
      time_MBOOST_L2 = time_MBOOST_L2,
      time_MBOOST_T_L2 = time_MBOOST_T_L2,
      time_RF = time_RF,
      sim_obj = sim_obj

    )
    saveRDS(results_save, file = sprintf("results/results_%04d.rds", i))
    
  }
}




