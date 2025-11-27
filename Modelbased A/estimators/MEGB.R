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


for (i in r) {
  print(i)
  no <- i

  path <- paste0(mod, "/results/",
                 "MEGB/")
  dir.create(path, showWarnings = F)
  # gradient process ----------------------------------------------------------
  if(file.exists(paste0(path, i, ".txt"))){
    warning(paste("Run no.", i, "has already been completed, skipping"))
    next
  } else{
    write.table(x = i, paste0(path, i, ".txt"))
    
    input.elements.no <- input_elements[[no]]
    input.elements.no$pop[,group.id] <- as.factor(input.elements.no$pop[,group.id])
	input.elements.no$sample[,group.id] <- as.factor(input.elements.no$sample[,group.id])

    megb_obj <- MEGB::megb(Y = input.elements.no$sample[, y.names],
                               X =  input.elements.no$sample[, x.names, drop = F],
                               dom_name = group.id,
                               smp_data = input.elements.no$sample,
                               pop_data = input.elements.no$pop,
                               gradient_params = params,
                               seed = 1,
                               mse = F,
                               bootstrap_cores = 1,
                               B = 100)

    megb_results <- megb_obj$Indicators
    names(megb_results) <- c("domain", "Mean")
    megb_results$run <- no
    
    csv_path <- paste0(path, "/", sprintf("%03d", no), "_gradientEstimates.csv")
    message(paste("csv_path:", csv_path))
    
    write.table(
      megb_results,
      file = csv_path,
      row.names = F,
      col.names = T
    )
    
    mse_results <- megb_obj$MSE_Estimates$MSE_estimates
    
    csv_path <- paste0(path, "/", sprintf("%03d", no), "_mse_results.csv")
    message(paste("csv_path:", csv_path))
    
    write.table(
      mse_results,
      file = csv_path,
      row.names = F,
      col.names = T
    )


sigma_e <- megb_obj$error_sd
    
    csv_path <- paste0(path, "/", sprintf("%03d", no), "_sigma_e.csv")
    message(paste("csv_path:", csv_path))
    
    write.table(
      sigma_e,
      file = csv_path,
      row.names = F,
      col.names = T
    )

    
  }
}

